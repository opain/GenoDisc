#!/usr/bin/Rscript
# ============================================================================
# compare_cojo.R  (INTERNAL - not part of the packaged pipeline)
#
# Compare GCTA-COJO and Manc-COJO selections from the SAME pipeline run, per
# GWAS and per chromosome. Both methods must have been run (config: cojo: T and
# manc_cojo: T), writing to results/<gwas>/cojo/ and results/<gwas>/manc_cojo/.
#
# Produces, under <out>/:
#   cojo_comparison_by_chr.tsv     one row per (gwas, chr): counts, Jaccard,
#                                  shared-bJ correlation, per-method status.
#   cojo_comparison_nonshared.tsv  one row per non-shared selected SNP: distance
#                                  to the nearest SNP selected by the other method
#                                  and the LD r2 with it (same reference panel).
#   cojo_comparison_summary.md     human-readable summary incl. timing/memory.
#
# LD r2 is computed by shelling out to PLINK 1.9 (bare `plink`, as elsewhere in
# the pipeline - see scripts/finemap.R) over the same 1KG+HGDP reference the
# COJO rules use: <resdir>/data/1kg_hgdp/1KG.Phase3.<pop>.MAF_001.chr<chr>.
#
# Usage:
#   Rscript misc/cojo_comparison/compare_cojo.R \
#     --outdir tests/output_full --gwas COAD01 --resdir <resdir> [--population EUR]
# Env fallbacks (mirrors misc/internal_results_review/_paths.R):
#   GENODISC_OUTDIR, GENODISC_RESDIR, GENODISC_GWAS, GENODISC_OUT
# ============================================================================

suppressMessages({
  library(optparse)
  library(data.table)
})

option_list <- list(
  make_option("--outdir", type = "character",
              default = Sys.getenv("GENODISC_OUTDIR", unset = NA),
              help = "Pipeline output dir containing results/<gwas>/{cojo,manc_cojo}/ [required]"),
  make_option("--gwas", type = "character",
              default = Sys.getenv("GENODISC_GWAS", unset = NA),
              help = "GWAS name(s), comma-separated [required]"),
  make_option("--resdir", type = "character",
              default = Sys.getenv("GENODISC_RESDIR", unset = NA),
              help = "Resource dir holding data/1kg_hgdp/ reference bfiles [required]"),
  make_option("--population", type = "character", default = "EUR",
              help = "Reference population for the LD panel [default EUR]"),
  make_option("--chromosomes", type = "character", default = "1:22",
              help = "Chromosomes, as an R-style range '1:22' or comma list [default 1:22]"),
  make_option("--plink", type = "character", default = "plink",
              help = "PLINK 1.9 executable (bare name assumes it is on PATH) [default plink]"),
  make_option("--out", type = "character",
              default = Sys.getenv("GENODISC_OUT", unset = NA),
              help = "Output dir [default <outdir>/cojo_comparison]")
)
opt <- parse_args(OptionParser(option_list = option_list))

stopifnot_arg <- function(x, nm) if (is.na(x) || x == "") stop(paste0("--", nm, " is required"), call. = FALSE)
stopifnot_arg(opt$outdir, "outdir"); stopifnot_arg(opt$gwas, "gwas"); stopifnot_arg(opt$resdir, "resdir")

gwas_names <- trimws(strsplit(opt$gwas, ",")[[1]])
chroms <- if (grepl(":", opt$chromosomes)) {
  rr <- as.integer(strsplit(opt$chromosomes, ":")[[1]]); rr[1]:rr[2]
} else as.integer(trimws(strsplit(opt$chromosomes, ",")[[1]]))
out_dir <- if (is.na(opt$out) || opt$out == "") file.path(opt$outdir, "cojo_comparison") else opt$out
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

jma_path    <- function(g, method, chr) file.path(opt$outdir, "results", g, method, paste0(g, "_chr", chr, ".jma.cojo"))
status_path <- function(g, method) {
  infix <- if (method == "manc_cojo") "manc_cojo" else "cojo"
  file.path(opt$outdir, "results", g, method, paste0(g, ".GW.", infix, ".status.csv"))
}
bench_path  <- function(nm) file.path(opt$outdir, "benchmarks", nm)
ref_prefix  <- function(chr) file.path(opt$resdir, "data", "1kg_hgdp",
                                       paste0("1KG.Phase3.", opt$population, ".MAF_001.chr", chr))

read_jma <- function(p) {
  if (!file.exists(p) || file.info(p)$size == 0) return(NULL)
  dt <- tryCatch(fread(p), error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0) return(NULL)
  # GCTA effect allele = refA; Manc-COJO effect allele = A1. Normalise to EA.
  if ("refA" %in% names(dt)) dt[, EA := refA] else if ("A1" %in% names(dt)) dt[, EA := A1] else dt[, EA := NA_character_]
  dt[, .(SNP, bp, EA, bJ)]
}

# Pairwise LD r2 over the reference panel for all `snps` on chromosome `chr`.
# One PLINK call per chr; returns a data.table (SNP_A, SNP_B, R2) or empty.
ld_r2_matrix <- function(chr, snps) {
  snps <- unique(snps)
  if (length(snps) < 2) return(data.table(SNP_A = character(), SNP_B = character(), R2 = numeric()))
  ref <- ref_prefix(chr)
  if (!file.exists(paste0(ref, ".bed"))) {
    warning("Reference bfile missing for chr", chr, " (", ref, "); r2 set NA")
    return(NULL)
  }
  tmp <- tempfile(paste0("ld_chr", chr, "_")); on.exit(unlink(paste0(tmp, c(".snps", ".ld", ".log", ".nosex")), force = TRUE), add = TRUE)
  writeLines(snps, paste0(tmp, ".snps"))
  args <- c("--bfile", ref, "--extract", paste0(tmp, ".snps"),
            "--r2", "--ld-window-kb", "1000000", "--ld-window", "999999",
            "--ld-window-r2", "0", "--out", tmp)
  res <- suppressWarnings(system2(opt$plink, args, stdout = FALSE, stderr = FALSE))
  if (!file.exists(paste0(tmp, ".ld"))) return(data.table(SNP_A = character(), SNP_B = character(), R2 = numeric()))
  ld <- tryCatch(fread(paste0(tmp, ".ld")), error = function(e) NULL)
  if (is.null(ld) || nrow(ld) == 0) return(data.table(SNP_A = character(), SNP_B = character(), R2 = numeric()))
  ld[, .(SNP_A, SNP_B, R2)]
}

# Look up r2 between two SNPs in an (undirected) LD table.
lookup_r2 <- function(ld, a, b) {
  if (is.null(ld) || nrow(ld) == 0) return(NA_real_)
  hit <- ld[(SNP_A == a & SNP_B == b) | (SNP_A == b & SNP_B == a), R2]
  if (length(hit) == 0) NA_real_ else hit[1]
}

read_bench <- function(nm) {
  p <- bench_path(nm)
  if (!file.exists(p)) return(NULL)
  tryCatch(fread(p), error = function(e) NULL)
}
# Sum wall-seconds / take max RSS (MB) across a set of benchmark files.
bench_agg <- function(names) {
  rows <- rbindlist(lapply(names, read_bench), fill = TRUE)
  if (nrow(rows) == 0) return(list(s = NA_real_, max_rss = NA_real_, n = 0L))
  list(s = sum(rows$s, na.rm = TRUE),
       max_rss = if ("max_rss" %in% names(rows)) max(rows$max_rss, na.rm = TRUE) else NA_real_,
       n = nrow(rows))
}

by_chr <- list(); nonshared <- list(); summaries <- list()

for (g in gwas_names) {
  gcta_status <- if (file.exists(status_path(g, "cojo"))) fread(status_path(g, "cojo")) else NULL
  manc_status <- if (file.exists(status_path(g, "manc_cojo"))) fread(status_path(g, "manc_cojo")) else NULL

  all_shared_bJ <- data.table(bJ_gcta = numeric(), bJ_manc_aligned = numeric())

  for (chr in chroms) {
    gcta <- read_jma(jma_path(g, "cojo", chr))
    manc <- read_jma(jma_path(g, "manc_cojo", chr))
    n_g <- if (is.null(gcta)) 0L else nrow(gcta)
    n_m <- if (is.null(manc)) 0L else nrow(manc)
    g_snps <- if (is.null(gcta)) character(0) else gcta$SNP
    m_snps <- if (is.null(manc)) character(0) else manc$SNP
    shared <- intersect(g_snps, m_snps)
    uni <- union(g_snps, m_snps)
    jacc <- if (length(uni) == 0) NA_real_ else length(shared) / length(uni)

    # Shared-SNP bJ, aligned to GCTA's effect allele.
    bJ_cor <- NA_real_
    if (length(shared) >= 1) {
      mg <- merge(gcta[SNP %in% shared], manc[SNP %in% shared], by = "SNP", suffixes = c("_g", "_m"))
      mg[, bJ_m_aligned := ifelse(EA_g == EA_m, bJ_m, -bJ_m)]
      all_shared_bJ <- rbind(all_shared_bJ, data.table(bJ_gcta = mg$bJ_g, bJ_manc_aligned = mg$bJ_m_aligned))
      if (nrow(mg) >= 3) bJ_cor <- suppressWarnings(cor(mg$bJ_g, mg$bJ_m_aligned))
    }

    # Non-shared SNPs: nearest other-method SNP + LD r2.
    only_g <- setdiff(g_snps, m_snps); only_m <- setdiff(m_snps, g_snps)
    if (length(only_g) + length(only_m) > 0) {
      ld <- ld_r2_matrix(chr, uni)
      add_nonshared <- function(src, src_dt, other_dt, method_lbl) {
        for (s in src) {
          bp_s <- src_dt[SNP == s, bp][1]
          if (is.null(other_dt) || nrow(other_dt) == 0) {
            nearest <- NA_character_; dist <- NA_real_; r2 <- NA_real_
          } else {
            d <- abs(other_dt$bp - bp_s); k <- which.min(d)
            nearest <- other_dt$SNP[k]; dist <- d[k]; r2 <- lookup_r2(ld, s, nearest)
          }
          nonshared[[length(nonshared) + 1L]] <<- data.table(
            gwas = g, chr = chr, selected_by = method_lbl, snp = s, bp = bp_s,
            nearest_other_snp = nearest, dist_bp = dist, r2_with_nearest = r2)
        }
      }
      add_nonshared(only_g, gcta, manc, "gcta")
      add_nonshared(only_m, manc, gcta, "manc")
    }

    gst <- if (!is.null(gcta_status)) gcta_status$status[gcta_status$chr == chr][1] else NA_character_
    mst <- if (!is.null(manc_status)) manc_status$status[manc_status$chr == chr][1] else NA_character_
    by_chr[[length(by_chr) + 1L]] <- data.table(
      gwas = g, chr = chr, n_gcta = n_g, n_manc = n_m, n_shared = length(shared),
      jaccard = round(jacc, 4), bJ_cor = round(bJ_cor, 4),
      gcta_status = gst, manc_status = mst)
  }

  # Genome-wide roll-up + timing/memory for this GWAS.
  gcta_bench <- bench_agg(paste0("cojo_", g, "_chr", chroms, ".tsv"))
  manc_bench <- bench_agg(paste0("manc_cojo_", g, "_chr", chroms, ".tsv"))
  gcta_proc <- bench_agg(paste0("process_cojo_", g, ".tsv"))
  manc_proc <- bench_agg(paste0("process_manc_cojo_", g, ".tsv"))
  # Manc per-chr peak memory from the run logs.
  manc_peak <- NA_real_
  peaks <- c()
  for (chr in chroms) {
    lg <- file.path(opt$outdir, "logs", paste0("manc_cojo-", g, "-chr", chr, ".log"))
    if (file.exists(lg)) {
      ln <- grep("Peak memory", readLines(lg, warn = FALSE), value = TRUE)
      if (length(ln)) peaks <- c(peaks, as.numeric(sub(".*Peak memory:\\s*([0-9.]+)\\s*GB.*", "\\1", ln[1])))
    }
  }
  if (length(peaks)) manc_peak <- max(peaks, na.rm = TRUE)

  bc <- rbindlist(by_chr)[gwas == g]
  gw_bJ_cor <- if (nrow(all_shared_bJ) >= 3) round(suppressWarnings(cor(all_shared_bJ$bJ_gcta, all_shared_bJ$bJ_manc_aligned)), 4) else NA_real_
  tot_shared <- sum(bc$n_shared); tot_union <- sum(bc$n_gcta) + sum(bc$n_manc) - tot_shared
  summaries[[g]] <- list(
    n_gcta = sum(bc$n_gcta), n_manc = sum(bc$n_manc), n_shared = tot_shared,
    jaccard = if (tot_union > 0) round(tot_shared / tot_union, 4) else NA_real_,
    bJ_cor = gw_bJ_cor,
    gcta_ref_too_small = if (!is.null(gcta_status)) gcta_status[status == "reference_too_small", paste(chr, collapse = ",")] else "",
    manc_exceeds = if (!is.null(manc_status)) manc_status[status == "exceeds_reference_n", paste(chr, collapse = ",")] else "",
    gcta_s = gcta_bench$s, manc_s = manc_bench$s,
    gcta_rss = gcta_bench$max_rss, manc_rss = manc_bench$max_rss,
    gcta_proc_s = gcta_proc$s, manc_proc_s = manc_proc$s,
    manc_peak_gb = manc_peak)
}

by_chr_dt <- rbindlist(by_chr)
nonshared_dt <- if (length(nonshared)) rbindlist(nonshared) else
  data.table(gwas = character(), chr = integer(), selected_by = character(), snp = character(),
             bp = numeric(), nearest_other_snp = character(), dist_bp = numeric(), r2_with_nearest = numeric())

fwrite(by_chr_dt,   file.path(out_dir, "cojo_comparison_by_chr.tsv"),   sep = "\t")
fwrite(nonshared_dt, file.path(out_dir, "cojo_comparison_nonshared.tsv"), sep = "\t")

# ---- Markdown summary -------------------------------------------------------
md <- c("# GCTA-COJO vs Manc-COJO comparison", "",
        paste0("Pipeline outdir: `", opt$outdir, "`  ",
               "Reference: `", opt$population, "` (1KG+HGDP)  GWAS: ", paste(gwas_names, collapse = ", ")),
        "")
for (g in gwas_names) {
  s <- summaries[[g]]
  spdup <- if (!is.na(s$gcta_s) && !is.na(s$manc_s) && s$manc_s > 0) round(s$gcta_s / s$manc_s, 1) else NA
  md <- c(md,
    paste0("## ", g), "",
    paste0("- Selected SNPs: GCTA ", s$n_gcta, ", Manc ", s$n_manc, ", shared ", s$n_shared,
           " (genome-wide Jaccard ", s$jaccard, ")"),
    paste0("- Shared-SNP joint-effect (bJ) correlation, allele-aligned: ", s$bJ_cor),
    paste0("- GCTA chromosomes flagged `reference_too_small`: ", ifelse(s$gcta_ref_too_small == "", "none", s$gcta_ref_too_small)),
    paste0("- Manc chromosomes flagged `exceeds_reference_n` (kept): ", ifelse(s$manc_exceeds == "", "none", s$manc_exceeds)),
    paste0("- Wall time (sum over chr): GCTA ", round(s$gcta_s, 1), "s vs Manc ", round(s$manc_s, 1),
           "s", ifelse(is.na(spdup), "", paste0("  (~", spdup, "x faster)"))),
    paste0("- process step wall time: GCTA ", round(s$gcta_proc_s, 1), "s vs Manc ", round(s$manc_proc_s, 1), "s"),
    paste0("- Max RSS (benchmarks): GCTA ", round(s$gcta_rss, 1), " MB vs Manc ", round(s$manc_rss, 1), " MB"),
    paste0("- Manc peak memory (from .log): ", ifelse(is.na(s$manc_peak_gb), "n/a", paste0(s$manc_peak_gb, " GB"))),
    "")
  # per-chr status differences
  bc <- by_chr_dt[gwas == g & (n_gcta != n_manc | gcta_status != manc_status)]
  if (nrow(bc)) {
    md <- c(md, "Chromosomes where counts or status differ:", "",
            "| chr | n_gcta | n_manc | shared | jaccard | bJ_cor | gcta_status | manc_status |",
            "|----:|-------:|-------:|-------:|--------:|-------:|-------------|-------------|")
    for (i in seq_len(nrow(bc))) md <- c(md, paste0("| ", paste(bc[i, .(chr, n_gcta, n_manc, n_shared, jaccard, bJ_cor, gcta_status, manc_status)], collapse = " | "), " |"))
    md <- c(md, "")
  }
  nd <- nonshared_dt[gwas == g]
  if (nrow(nd)) {
    md <- c(md, paste0("Non-shared selected SNPs: ", nrow(nd),
                       " (median LD r2 with nearest other-method SNP: ",
                       round(median(nd$r2_with_nearest, na.rm = TRUE), 3),
                       "; median distance ", round(median(nd$dist_bp, na.rm = TRUE) / 1000, 1), " kb)."), "")
  }
}
writeLines(md, file.path(out_dir, "cojo_comparison_summary.md"))

cat("Wrote:\n  ", file.path(out_dir, "cojo_comparison_by_chr.tsv"),
    "\n  ", file.path(out_dir, "cojo_comparison_nonshared.tsv"),
    "\n  ", file.path(out_dir, "cojo_comparison_summary.md"), "\n", sep = "")
