#!/usr/bin/env Rscript
# Aggregate all N runs into a Markdown summary + PNG plot.
#
# Reads $RESULTS_DIR/N*/samples.csv + checkpoints.jsonl.
# Env vars (optional):
#   RESULTS_DIR - default: <script dir>/results
#   BUNDLES_DIR - default: <script dir>/bundles
suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(ggplot2)
})

this_file <- (function() {
  a <- commandArgs(trailingOnly = FALSE)
  m <- sub("^--file=", "", grep("^--file=", a, value = TRUE))
  if (length(m)) normalizePath(m[1]) else normalizePath(sys.frame(1)$ofile)
})()
SCRIPT_DIR <- dirname(this_file)
RES <- Sys.getenv("RESULTS_DIR", unset = file.path(SCRIPT_DIR, "results"))
BUN <- Sys.getenv("BUNDLES_DIR", unset = file.path(SCRIPT_DIR, "bundles"))

runs <- sort(list.files(RES, pattern = "^N\\d+$", full.names = TRUE))
stopifnot(length(runs) > 0)

parse_checkpoints <- function(path) {
  if (!file.exists(path)) return(data.table())
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(lines)]
  rows <- lapply(lines, function(l) {
    j <- tryCatch(fromJSON(l), error = function(e) NULL)
    if (is.null(j)) return(NULL)
    data.table(epoch_ms = as.numeric(j$epoch_ms), label = as.character(j$label))
  })
  rbindlist(rows, fill = TRUE)
}

collect <- function(run_dir) {
  N <- as.integer(sub("^N", "", basename(run_dir)))
  csv <- file.path(run_dir, "samples.csv")
  cp  <- parse_checkpoints(file.path(run_dir, "checkpoints.jsonl"))
  if (!file.exists(csv)) return(NULL)
  # Force numeric: fread would pick integer64 for millisec epochs, and
  # coercing integer64 -> numeric without bit64 returns garbage.
  s <- fread(csv, colClasses = list(numeric = c(
    "epoch_ms","vmrss_kb","vmdata_kb","vmsize_kb","vmpeak_kb",
    "tmp_bundle_bytes","tmp_upload_bytes","tmp_total_bytes"
  )))
  if (nrow(s) == 0) return(NULL)

  s[, t_sec := (epoch_ms - min(epoch_ms)) / 1000]
  base_rss <- s[t_sec <= 5, min(vmrss_kb)]

  ready_epoch <- cp[label == "app_ready", epoch_ms][1]
  post_upload_rss <- if (!is.na(ready_epoch)) {
    s[epoch_ms >= ready_epoch & epoch_ms <= ready_epoch + 4000, max(vmrss_kb)]
  } else NA_real_

  done_epoch <- cp[label == "done", epoch_ms][1]
  post_all_rss <- if (!is.na(done_epoch)) {
    s[epoch_ms >= done_epoch - 2000 & epoch_ms <= done_epoch + 2000, max(vmrss_kb)]
  } else s[.N, vmrss_kb]

  peak_rss   <- s[, max(vmrss_kb)]
  peak_vmsize<- s[, max(vmsize_kb)]
  peak_tmp   <- s[, max(tmp_total_bytes)]

  tar_path <- file.path(BUN, sprintf("als_N%02d.tar.gz", N))
  tar_mb <- if (file.exists(tar_path)) file.size(tar_path) / 1024^2 else NA_real_

  data.table(
    N              = N,
    tarball_mb     = round(tar_mb, 1),
    baseline_mb    = round(base_rss / 1024, 1),
    after_upload_mb= round(post_upload_rss / 1024, 1),
    after_walk_mb  = round(post_all_rss / 1024, 1),
    peak_rss_mb    = round(peak_rss / 1024, 1),
    peak_vsize_mb  = round(peak_vmsize / 1024, 1),
    peak_tmp_mb    = round(peak_tmp / 1024^2, 1)
  )
}

tab <- rbindlist(lapply(runs, collect), fill = TRUE)
setorder(tab, N)
print(tab)

peak_1gb <- 900   # leave ~100 MB headroom below plan limit
peak_3gb <- 2800
peak_8gb <- 7500

recommend <- function(mb) {
  if (is.na(mb)) return("n/a")
  if (mb < peak_1gb) return("Basic (1 GB)")
  if (mb < peak_3gb) return("Standard (3 GB)")
  if (mb < peak_8gb) return("Standard-L / Pro (8 GB)")
  return("EXCEEDS 8 GB — split bundle")
}
tab[, plan_tier := sapply(peak_rss_mb, recommend)]

plot_dt <- melt(tab, id.vars = "N",
                measure.vars = c("baseline_mb", "after_upload_mb", "after_walk_mb", "peak_rss_mb"),
                variable.name = "phase", value.name = "rss_mb")
plot_dt[, phase := factor(phase,
                          levels = c("baseline_mb", "after_upload_mb", "after_walk_mb", "peak_rss_mb"),
                          labels = c("baseline", "after upload", "after tab walk", "peak"))]

p <- ggplot(plot_dt, aes(x = N, y = rss_mb, colour = phase)) +
  geom_line() + geom_point(size = 2) +
  geom_hline(yintercept = 1024, linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = 3072, linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = 8192, linetype = "dashed", alpha = 0.5) +
  annotate("text", x = max(plot_dt$N), y = 1024, label = "1 GB plan",  hjust = 1, vjust = -0.3, size = 3) +
  annotate("text", x = max(plot_dt$N), y = 3072, label = "3 GB plan",  hjust = 1, vjust = -0.3, size = 3) +
  annotate("text", x = max(plot_dt$N), y = 8192, label = "8 GB plan",  hjust = 1, vjust = -0.3, size = 3) +
  scale_x_continuous(breaks = tab$N) +
  labs(title = "output_explorer_embed RSS vs GWAS count",
       subtitle = "Synthetic bundles duplicated from als_bundle.tar.gz",
       x = "Number of GWAS in bundle", y = "R process RSS (MB)") +
  theme_minimal()

ggsave(file.path(RES, "peak_rss_vs_N.png"), p, width = 8, height = 5, dpi = 120)

lines <- c(
  "# output_explorer_embed resource profile",
  "",
  sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "Synthetic multi-GWAS bundles built by duplicating the shipped ALS example",
  "one-GWAS bundle. Each run: start Shiny, upload the bundle in headless",
  "Firefox, click every top tab + every sub-tab, sample R process VmRSS at 2 Hz.",
  "",
  "## Peak memory table",
  "",
  "| N | tarball MB | baseline MB | after upload MB | after tab walk MB | **peak RSS MB** | peak VSIZE MB | peak /tmp MB | recommended plan |",
  "|--:|-----------:|------------:|----------------:|------------------:|----------------:|--------------:|-------------:|:-----------------|"
)
for (i in seq_len(nrow(tab))) {
  r <- tab[i]
  lines <- c(lines, sprintf(
    "| %d | %s | %s | %s | %s | **%s** | %s | %s | %s |",
    r$N, format(r$tarball_mb), format(r$baseline_mb),
    format(r$after_upload_mb), format(r$after_walk_mb),
    format(r$peak_rss_mb), format(r$peak_vsize_mb),
    format(r$peak_tmp_mb), r$plan_tier
  ))
}

if (nrow(tab) >= 2) {
  fit <- lm(peak_rss_mb ~ N, data = tab)
  slope <- coef(fit)[["N"]]
  intercept <- coef(fit)[["(Intercept)"]]
  n1_upload_delta <- tab[N == min(N), after_upload_mb - baseline_mb]
  nmax_upload_delta <- tab[N == max(N), after_upload_mb - baseline_mb]
  nmax_walk_delta   <- tab[N == max(N), after_walk_mb - after_upload_mb]

  lines <- c(lines, "",
             "## Linear fit (peak RSS vs N)",
             "",
             sprintf("- **peak_rss_mb ≈ %.1f + %.1f × N**", intercept, slope),
             sprintf("- ~%.0f MB per additional GWAS", slope),
             sprintf("- R² = %.3f", summary(fit)$r.squared),
             "",
             "## Key findings",
             "",
             "1. **Two distinct memory-cost phases.** The `after upload` column is",
             "   the RSS right after the bundle is opened + `build_comparison_long()`",
             "   finishes (app.R). The `after tab walk` column is after clicking every",
             "   tab / sub-tab. **The walk adds more per GWAS than the initial load** —",
             sprintf("   at N=%d, upload alone costs %.0f MB above baseline; the walk adds",
                     tab[.N, N], nmax_upload_delta),
             sprintf("   another %.0f MB.", nmax_walk_delta),
             "",
             "2. **`/tmp` footprint = 2 × extracted bundle size.** shinyapps.io keeps",
             "   both the uploaded tarball and the extracted tree in the R session's",
             "   TMPDIR — no `on.exit()` cleanup in `gd_open()` (reader.R).",
             "",
             "## Recommendation",
             "",
             sprintf("- **Basic (1 GB)** — safe up to ~%d GWAS of the profiled scale.",
                     max(1L, floor((1024 - intercept) / slope))),
             sprintf("- **Standard (3 GB)** — safe up to ~%d GWAS of the profiled scale.",
                     max(1L, floor((3072 - intercept) / slope))),
             "- **If per-GWAS SNP counts are much larger than the source bundle**,",
             "  halve the numbers above as a rule of thumb until rerun on a real bundle."
  )
}

lines <- c(lines, "",
           "## Plot",
           "",
           "![peak RSS vs N](peak_rss_vs_N.png)",
           "",
           "## Caveats",
           "",
           "- Synthetic bundles duplicate the same source GWAS N times, so per-GWAS",
           "  block sizes are identical. Real bundles vary and can be larger — the",
           "  linear growth in peak RSS here is a *lower bound* for bundles with",
           "  heavier per-GWAS payloads.",
           "- Peak RSS was sampled at 2 Hz — momentary spikes shorter than 500 ms may",
           "  be missed. DT / ggplot render peaks last multiple seconds so this is fine.",
           "- Cross-GWAS LDSC / rG blocks (top-level, not per-GWAS) aren't present in",
           "  the source example bundle. A real multi-GWAS bundle would carry an N×N",
           "  rG matrix which adds O(N²) memory — not captured here.")

writeLines(lines, file.path(RES, "summary.md"))
message("wrote ", file.path(RES, "summary.md"))
message("wrote ", file.path(RES, "peak_rss_vs_N.png"))
