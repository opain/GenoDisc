#!/usr/bin/Rscript
suppressMessages(library("optparse"))

option_list = list(
  make_option("--gwas", action="store", default=NA, type='character',
              help="Name of GWAS [required]"),
  make_option("--config_file", action="store", default=NA, type='character',
              help="Path to config file [required]")
)

option_list <- c(option_list, list(
  make_option("--pipeline_dir", action="store", default=NA, type="character",
              help="Path to the pipeline directory [required]")
))

opt = parse_args(OptionParser(option_list=option_list))
options(pipeline_dir = opt$pipeline_dir)

suppressMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scattermore)
})
source(file.path(opt$pipeline_dir, 'scripts', 'functions', 'utils_functions.R'))

outdir  <- read_param(config = opt$config_file, param = 'outdir',  return_obj = F)
resdir  <- read_param(config = opt$config_file, param = 'resdir',  return_obj = F)
gwas_list <- read_param(config = opt$config_file, param = 'gwas_list', return_obj = T)
population <- gwas_list$population[gwas_list$name == opt$gwas]

locus_dir <- paste0(outdir, '/results/', opt$gwas, '/locus_plots/')
dir.create(locus_dir, recursive = TRUE, showWarnings = FALSE)
sentinel  <- paste0(locus_dir, opt$gwas, '.locus_plots.done')

clump_csv <- paste0(outdir, '/results/', opt$gwas, '/clump/', opt$gwas, '.GW.clump.clean.csv')
idx <- fread(clump_csv)
idx <- idx[!is.na(P) & P > 0 & P < 1e-5]
setorder(idx, P)
if (nrow(idx) > 50) idx <- idx[1:50]

if (nrow(idx) == 0) {
  message("No index variants with P < 1e-5; writing sentinel and exiting.")
  file.create(sentinel)
  quit(save = "no", status = 0)
}

# Parse .clumped per-chromosome files to build index SNP -> member SNP map
clumped_files <- Sys.glob(paste0(outdir, '/results/', opt$gwas, '/clump/', opt$gwas, '_chr*.clumped'))
clump_members_by_index <- list()
for (f in clumped_files) {
  tmp <- tryCatch(fread(f), error = function(e) NULL)
  if (is.null(tmp) || !all(c("SNP", "SP2") %in% names(tmp))) next
  for (i in seq_len(nrow(tmp))) {
    sp2 <- tmp$SP2[i]
    if (is.na(sp2) || sp2 == "NONE") {
      clump_members_by_index[[tmp$SNP[i]]] <- character(0)
    } else {
      members <- unlist(strsplit(sp2, ","))
      members <- gsub("\\(.*$", "", members)
      clump_members_by_index[[tmp$SNP[i]]] <- members
    }
  }
}

ss <- fread(paste0(outdir, '/results/', opt$gwas, '/gwas_sumstat/', opt$gwas, '.cleaned.gz'),
            select = c("CHR", "BP", "SNP", "P"))
ss <- ss[!is.na(P) & P > 0 & P <= 1 & !is.na(BP) & !is.na(CHR)]
ss[, CHR := as.integer(CHR)]

gene_locs <- fread(paste0(resdir, '/data/biomart/gene_locations.tsv'))
gene_locs <- gene_locs[chromosome_name %in% as.character(1:22)]
gene_locs[, chromosome_name := as.integer(chromosome_name)]

r2_breaks <- c(-Inf, 0.2, 0.4, 0.6, 0.8, Inf)
r2_labels <- c("[0,0.2)", "[0.2,0.4)", "[0.4,0.6)", "[0.6,0.8)", "[0.8,1.0]")
r2_colours <- c(
  "NA"        = "grey70",
  "[0,0.2)"   = "#357EBD",
  "[0.2,0.4)" = "#46B8DA",
  "[0.4,0.6)" = "#5CB85C",
  "[0.6,0.8)" = "#F0AD4E",
  "[0.8,1.0]" = "#D43F3A"
)
index_colour <- "#9632B8"

stack_genes <- function(g) {
  if (nrow(g) == 0) return(g)
  g <- g[order(eff_start)]
  row_ends <- numeric(0)
  rows <- integer(nrow(g))
  for (i in seq_len(nrow(g))) {
    placed <- FALSE
    for (r in seq_along(row_ends)) {
      if (g$eff_start[i] > row_ends[r]) {
        rows[i] <- r
        row_ends[r] <- g$eff_end[i]
        placed <- TRUE
        break
      }
    }
    if (!placed) {
      row_ends <- c(row_ends, g$eff_end[i])
      rows[i] <- length(row_ends)
    }
  }
  g[, gene_row := rows]
  g
}

sanitise_id <- function(x) gsub("[/:]", "_", x)

make_plot <- function(idx_row) {
  snp_id  <- idx_row$SNP
  chr_i   <- as.integer(idx_row$CHR)
  bp_i    <- as.integer(idx_row$BP)
  p_i     <- as.numeric(idx_row$P)
  nearest <- sub(",.*$", "", idx_row$NearestGene)

  members <- clump_members_by_index[[snp_id]]
  if (is.null(members)) members <- character(0)
  member_bps <- ss[CHR == chr_i & SNP %in% members, BP]
  all_bps <- c(bp_i, member_bps)
  min_clump_bp <- min(all_bps)
  max_clump_bp <- max(all_bps)
  win_lo <- min_clump_bp - 500000L
  win_hi <- max_clump_bp + 500000L

  win <- ss[CHR == chr_i & BP >= win_lo & BP <= win_hi]
  if (nrow(win) == 0) {
    warning(sprintf("Locus %s: no variants in window; skipping.", snp_id))
    return(invisible(NULL))
  }

  bfile <- sprintf("%s/data/1kg/1KG.Phase3.%s.MAF_001.chr%d", resdir, population, chr_i)
  tmp_prefix <- tempfile(pattern = paste0("ld_", sanitise_id(snp_id), "_"))
  window_kb_half <- ceiling(max(bp_i - win_lo, win_hi - bp_i) / 1000)

  cmd <- sprintf(
    "plink --bfile %s --r2 --ld-snp %s --ld-window-kb %d --ld-window 99999 --ld-window-r2 0 --out %s",
    bfile, shQuote(snp_id), window_kb_half, tmp_prefix
  )
  rc <- system(paste(cmd, "> /dev/null 2>&1"))
  ld_path <- paste0(tmp_prefix, ".ld")
  if (rc != 0 || !file.exists(ld_path)) {
    warning(sprintf("Locus %s: PLINK failed or produced no .ld output; skipping.", snp_id))
    on.exit(unlink(paste0(tmp_prefix, "*")), add = TRUE)
    return(invisible(NULL))
  }
  ld <- tryCatch(fread(ld_path, select = c("SNP_B", "R2")),
                 error = function(e) NULL)
  unlink(paste0(tmp_prefix, "*"))
  if (is.null(ld) || nrow(ld) == 0) {
    warning(sprintf("Locus %s: empty LD output; skipping.", snp_id))
    return(invisible(NULL))
  }
  setnames(ld, c("SNP", "R2"))

  win <- merge(win, ld, by = "SNP", all.x = TRUE, sort = FALSE)
  win[, neglogP := -log10(P)]
  win[, ld_bin := cut(R2, breaks = r2_breaks, labels = r2_labels, right = FALSE, include.lowest = TRUE)]
  win[, ld_bin := as.character(ld_bin)]
  win[is.na(ld_bin), ld_bin := "NA"]
  win[, ld_bin := factor(ld_bin, levels = names(r2_colours))]
  win[, is_index := SNP == snp_id]
  win <- win[order(ld_bin)]

  bulk <- win[is_index == FALSE]
  scatter_geom <- if (nrow(bulk) > 5000) {
    geom_scattermore(data = bulk,
                     aes(x = BP, y = neglogP, colour = ld_bin),
                     pointsize = 2.5)
  } else {
    geom_point(data = bulk,
               aes(x = BP, y = neglogP, colour = ld_bin),
               size = 1.5)
  }

  title_text <- sprintf("%s — %s (chr%d:%s)", snp_id, nearest, chr_i, format(bp_i, big.mark = ","))

  scatter <- ggplot() +
    scatter_geom +
    geom_point(data = win[is_index == TRUE],
               aes(x = BP, y = neglogP),
               shape = 23, size = 4, fill = index_colour, colour = "black") +
    geom_hline(yintercept = -log10(5e-8), colour = "red", linetype = "solid") +
    geom_hline(yintercept = -log10(1e-5), colour = "blue", linetype = "dashed") +
    scale_colour_manual(values = r2_colours, drop = FALSE, name = expression(r^2)) +
    scale_x_continuous(labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
    coord_cartesian(xlim = c(win_lo, win_hi)) +
    expand_limits(y = 0) +
    labs(x = NULL,
         y = expression(-log[10](italic(P))),
         title = title_text) +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(size = 10, face = "bold"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())

  genes_win <- gene_locs[chromosome_name == chr_i &
                          start_position <= win_hi &
                          end_position   >= win_lo]
  if (nrow(genes_win) > 0) {
    win_width_bp <- win_hi - win_lo
    chars_per_window <- 110
    label_pad <- 1.35
    genes_win[, vis_start := pmax(start_position, win_lo)]
    genes_win[, vis_end   := pmin(end_position,   win_hi)]
    genes_win[, label_x   := (vis_start + vis_end) / 2]
    genes_win[, label_half_bp := nchar(external_gene_name) *
                                  (win_width_bp / chars_per_window) / 2 * label_pad]
    genes_win[, eff_start := pmin(start_position, label_x - label_half_bp)]
    genes_win[, eff_end   := pmax(end_position,   label_x + label_half_bp)]
    genes_win <- stack_genes(genes_win)
    genes_win[, arrow_x_start := ifelse(strand == 1, end_position, start_position)]
    genes_win[, arrow_x_end   := ifelse(strand == 1,
                                         pmin(end_position   + (win_hi - win_lo) * 0.01, win_hi),
                                         pmax(start_position - (win_hi - win_lo) * 0.01, win_lo))]
    n_rows <- max(genes_win$gene_row)
    gene_plot <- ggplot(genes_win) +
      geom_rect(aes(xmin = start_position, xmax = end_position,
                    ymin = gene_row - 0.15, ymax = gene_row + 0.15),
                fill = "grey40", colour = NA) +
      geom_segment(aes(x = arrow_x_start, xend = arrow_x_end,
                       y = gene_row, yend = gene_row),
                   arrow = arrow(length = unit(0.08, "cm"), type = "closed"),
                   colour = "grey20") +
      geom_text(aes(x = label_x, y = gene_row + 0.45, label = external_gene_name),
                size = 2.5, colour = "grey20") +
      scale_x_continuous(labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
      scale_y_continuous(limits = c(0.4, n_rows + 0.9)) +
      coord_cartesian(xlim = c(win_lo, win_hi)) +
      labs(x = sprintf("Chromosome %d position (bp)", chr_i), y = NULL) +
      theme_bw(base_size = 10) +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            panel.grid = element_blank())
  } else {
    gene_plot <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "No genes in window", colour = "grey50") +
      theme_void()
  }

  n_gene_rows <- if (exists("n_rows")) n_rows else 1
  scatter_h <- 4.5
  gene_h    <- max(1.2, 0.35 * n_gene_rows + 0.6)
  combined  <- scatter / gene_plot + plot_layout(heights = c(scatter_h, gene_h))
  out_png   <- paste0(locus_dir, opt$gwas, '.locus_plot.', sanitise_id(snp_id), '.png')
  ggsave(out_png, plot = combined, width = 8, height = scatter_h + gene_h + 0.5, dpi = 120)
  invisible(out_png)
}

for (i in seq_len(nrow(idx))) {
  tryCatch(make_plot(idx[i]),
           error = function(e) warning(sprintf("Locus %s: %s", idx$SNP[i], conditionMessage(e))))
}

file.create(sentinel)
