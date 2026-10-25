#!/usr/bin/Rscript
# Extract PER-VARIANT credible-set membership + PIP from the SuSiE fit objects.
#
# process_finemap.R writes only a per-CREDIBLE-SET summary (one row per cs). The
# coding-annotation (Phase 1) and CALDERA (Phase 3) steps need one row PER VARIANT,
# with its PIP and genomic coordinates/alleles. Per-variant PIP lives only in the
# saved susie_rss fit objects ({gwas}.chr{chr}.{lead}.rds [L=10] and .L1.rds [L=1]),
# where $pip is named by rsID. Coordinates + alleles come from the SAME per-chr
# reference .bim finemap.R used (the cleaned sumstats carry no CHR/BP), matched by
# rsID: V1=CHR, V2=SNP, V4=BP, V5=A1, V6=A2.
#
# Output: {gwas}.GW.credset_variants.csv with columns
#   locus, cs, L, SNP, CHR, BP, A1, A2, PIP
# A header-only file is written when no credible set exists.

suppressMessages(library("optparse"))
option_list <- list(
  make_option("--gwas", default = NA, type = "character"),
  make_option("--config_file", default = NA, type = "character"),
  make_option("--pipeline_dir", default = NA, type = "character"),
  make_option("--population", default = "EUR", type = "character",
              help = "GWAS population code; selects the ancestry-matched reference bim [default EUR]")
)
opt <- parse_args(OptionParser(option_list = option_list))
options(pipeline_dir = opt$pipeline_dir)
suppressMessages({ library(data.table); library(susieR); library(stringr) })
source(file.path(opt$pipeline_dir, "scripts", "functions", "utils_functions.R"))

outdir <- read_param(config = opt$config_file, param = "outdir", return_obj = FALSE)
resdir <- read_param(config = opt$config_file, param = "resdir", return_obj = FALSE)

out_cols <- c("locus", "cs", "L", "SNP", "CHR", "BP", "A1", "A2", "PIP")
empty_out <- function() setNames(data.table(matrix(character(0), ncol = length(out_cols))), out_cols)

fm_dir <- paste0(outdir, "/results/", opt$gwas, "/finemap/")
out_path <- paste0(fm_dir, opt$gwas, ".GW.credset_variants.csv")

rds_files <- list.files(path = fm_dir, pattern = "\\.rds$", full.names = TRUE)
if (length(rds_files) == 0) { fwrite(empty_out(), out_path); quit(save = "no") }

# Parse chr + lead + L from the filename ({gwas}.chr{chr}.{lead}[.L1].rds). The gwas
# name itself may contain 'chr' (e.g. COAD01chr22), so strip the gwas prefix first.
parse_name <- function(f) {
  b <- sub("\\.rds$", "", basename(f))
  L <- if (grepl("\\.L1$", b)) 1L else 10L
  b <- sub("\\.L1$", "", b)
  rest <- sub(paste0("^", opt$gwas, "\\.chr"), "", b)   # -> "{chr}.{lead}"
  chr <- sub("\\..*$", "", rest)
  lead <- sub("^[^.]*\\.", "", rest)
  list(chr = chr, lead = lead, L = L)
}

rows <- list()
for (f in rds_files) {
  meta <- parse_name(f)
  fit <- readRDS(f)
  cs <- summary(fit)$cs
  if (is.null(cs) || nrow(cs) == 0) next
  for (j in seq_len(nrow(cs))) {
    idx <- as.numeric(unlist(str_split(cs$variable[j], ",")))
    idx <- idx[!is.na(idx)]
    if (length(idx) == 0) next
    rows[[length(rows) + 1L]] <- data.table(
      locus = meta$lead, cs = cs$cs[j], L = meta$L, CHR = meta$chr,
      SNP = names(fit$pip)[idx], PIP = as.numeric(fit$pip[idx]))
  }
}
per_var <- rbindlist(rows, use.names = TRUE, fill = TRUE)
if (is.null(per_var) || nrow(per_var) == 0) { fwrite(empty_out(), out_path); quit(save = "no") }

# Join BP/A1/A2 from the per-chr reference bim (per chromosome present).
per_var[, BP := NA_integer_][, A1 := NA_character_][, A2 := NA_character_]
for (chr in unique(per_var$CHR)) {
  bim_path <- paste0(resdir, "/data/1kg_hgdp/1KG.Phase3.", opt$population, ".MAF_001.chr", chr, ".bim")
  if (!file.exists(bim_path)) next
  bim <- fread(bim_path, header = FALSE, select = c(2, 4, 5, 6))
  setnames(bim, c("SNP", "bBP", "bA1", "bA2"))
  sel <- per_var$CHR == chr
  m <- match(per_var$SNP[sel], bim$SNP)
  per_var$BP[sel] <- bim$bBP[m]
  per_var$A1[sel] <- bim$bA1[m]
  per_var$A2[sel] <- bim$bA2[m]
}

per_var <- per_var[, ..out_cols]
setorder(per_var, locus, cs, -PIP)
fwrite(per_var, out_path)
