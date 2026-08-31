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
  library(ggrepel)
})
source(file.path(opt$pipeline_dir, 'scripts', 'functions', 'utils_functions.R'))

outdir <- read_param(config = opt$config_file, param = 'outdir', return_obj = F)

ss_full <- fread(paste0(outdir, '/results/', opt$gwas, '/gwas_sumstat/', opt$gwas, '.cleaned.gz'),
                 select = c("CHR", "BP", "SNP", "P"))
ss_full <- ss_full[!is.na(P) & P > 0 & P <= 1 & !is.na(BP) & !is.na(CHR)]
ss_full[, CHR := as.integer(CHR)]
ss_full <- ss_full[CHR >= 1 & CHR <= 22]

chr_lens <- ss_full[, .(chr_max = max(BP)), by = CHR][order(CHR)]
chr_lens[, cum_offset := cumsum(as.numeric(chr_max)) - as.numeric(chr_max)]
chr_lens[, label_pos := cum_offset + as.numeric(chr_max) / 2]
chr_lens[, chr_colour := ifelse(CHR %% 2 == 1, "grey60", "grey30")]

ss <- ss_full[P < 1e-3]
if (nrow(ss) < 10) ss <- ss_full
ss <- merge(ss, chr_lens[, .(CHR, cum_offset, chr_colour)], by = "CHR")
ss[, cum_bp := BP + cum_offset]
ss[, neglogP := -log10(P)]

clump_csv <- paste0(outdir, '/results/', opt$gwas, '/clump/', opt$gwas, '.GW.clump.clean.csv')
has_clump <- file.exists(clump_csv)

idx_plot <- NULL
clump_members <- character(0)

if (has_clump) {
  idx <- fread(clump_csv, select = c("SNP", "CHR", "BP", "P", "NearestGene"))
  idx <- idx[!is.na(P) & P > 0 & CHR >= 1 & CHR <= 22]
  idx_plot <- merge(idx, chr_lens[, .(CHR, cum_offset)], by = "CHR")
  idx_plot[, cum_bp := BP + cum_offset]
  idx_plot[, neglogP := -log10(P)]
  idx_plot[, label := sub(",.*$", "", NearestGene)]
  idx_plot[label == "None", label := NA_character_]

  clumped_files <- Sys.glob(paste0(outdir, '/results/', opt$gwas, '/clump/', opt$gwas, '_chr*.clumped'))
  for (f in clumped_files) {
    tmp <- tryCatch(fread(f), error = function(e) NULL)
    if (is.null(tmp) || !("SP2" %in% names(tmp))) next
    sp2_strings <- tmp$SP2[!is.na(tmp$SP2) & tmp$SP2 != "NONE"]
    if (length(sp2_strings) == 0) next
    members <- unlist(strsplit(paste(sp2_strings, collapse = ","), ","))
    members <- gsub("\\(.*$", "", members)
    clump_members <- c(clump_members, members)
  }
  clump_members <- unique(clump_members)
}

ss[, category := "other"]
ss[SNP %in% clump_members, category := "clump_member"]
if (!is.null(idx_plot)) ss[SNP %in% idx_plot$SNP, category := "index"]

base_plot <- ggplot() +
  geom_point(data = ss[category == "other"],
             aes(x = cum_bp, y = neglogP, colour = chr_colour),
             size = 0.6) +
  scale_colour_identity() +
  geom_point(data = ss[category == "clump_member"],
             aes(x = cum_bp, y = neglogP),
             colour = "green3", size = 0.6) +
  geom_hline(yintercept = -log10(5e-8), colour = "red", linetype = "solid") +
  geom_hline(yintercept = -log10(1e-5), colour = "blue", linetype = "dashed") +
  scale_x_continuous(breaks = chr_lens$label_pos, labels = chr_lens$CHR,
                     expand = expansion(mult = c(0.01, 0.01))) +
  expand_limits(y = 0) +
  labs(x = "Chromosome", y = expression(-log[10](italic(P)))) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank())

if (!is.null(idx_plot) && nrow(idx_plot) > 0) {
  base_plot <- base_plot +
    geom_point(data = idx_plot, aes(x = cum_bp, y = neglogP),
               shape = 23, size = 3, fill = "darkgreen", colour = "black")
}

unlab_path <- paste0(outdir, '/results/', opt$gwas, '/gwas_sumstat/', opt$gwas, '.manhattan_plot.unlabelled.png')
lab_path   <- paste0(outdir, '/results/', opt$gwas, '/gwas_sumstat/', opt$gwas, '.manhattan_plot.labelled.png')

ggsave(unlab_path, plot = base_plot, width = 12, height = 5, dpi = 200)

if (!is.null(idx_plot) && nrow(idx_plot) > 0) {
  idx_labelled <- idx_plot[!is.na(label) & P < 5e-8]
  labelled_plot <- base_plot +
    geom_text_repel(data = idx_labelled,
                    aes(x = cum_bp, y = neglogP, label = label),
                    ylim = c(-log10(5e-8), NA),
                    direction = "both",
                    force = 2,
                    box.padding = 0.5,
                    point.padding = 0.3,
                    min.segment.length = 0,
                    max.overlaps = Inf,
                    segment.size = 0.3,
                    segment.colour = "grey50",
                    size = 3)
  ggsave(lab_path, plot = labelled_plot, width = 12, height = 5, dpi = 200)
} else {
  ggsave(lab_path, plot = base_plot, width = 12, height = 5, dpi = 200)
}

# Also write the raw filtered variant data (with clumping status + nearest
# gene for indexes) so the Shiny app can rebuild the plot on the fly and
# expose user-facing knobs (threshold, gene labels, theme, etc.) via
# build_manhattan_plot(). Much smaller than the two PNG blobs.
manhattan_data <- list(
  data = ss[, .(CHR, BP, SNP, P, chr_colour, cum_bp, neglogP, category)],
  offsets = chr_lens[, .(CHR, cum_offset, label_pos, chr_colour)],
  indexes = if (!is.null(idx_plot))
              idx_plot[, .(SNP, CHR, BP, P, cum_bp, neglogP, label)]
            else data.table(SNP = character(), CHR = integer(), BP = numeric(),
                            P = numeric(), cum_bp = numeric(),
                            neglogP = numeric(), label = character()),
  p_threshold = if (nrow(ss) < 10) 1 else 1e-3
)
data_path <- paste0(outdir, '/results/', opt$gwas, '/gwas_sumstat/',
                    opt$gwas, '.manhattan_data.rds')
saveRDS(manhattan_data, data_path)
