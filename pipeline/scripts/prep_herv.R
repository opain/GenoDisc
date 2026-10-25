#!/usr/bin/Rscript
# Prepare the HERV retrotranscriptome FUSION panel (Duarte, Pain et al. 2024,
# figshare 22179655) after download/extraction into the panel directory. The
# archive nests the weights under users/.../EUR_all_1/CMC_wrapped/ and contains
# BOTH gene (HGNC symbol) and HERV (FAMILY_cytoband) features. Per decision D2 we
# keep HERV features ONLY. Writes the canonical <panel>.pos (PANEL WGT ID CHR P0
# P1 N) with CMC_wrapped/ moved directly under the panel dir so WGT resolves with
# --weights_dir <panel_dir>. HERV IDs are kept verbatim (NO Ensembl remap).
suppressMessages(library(optparse)); suppressMessages(library(data.table))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--panel_dir", type = "character"),
  make_option("--panel", type = "character", default = "CMC.BRAIN.HERV"),
  make_option("--n", type = "integer", default = 563L, help = "EUR sample size")
)))
pd <- opt$panel_dir
out_pos <- file.path(pd, paste0(opt$panel, ".pos"))
HERV_RE <- "_[0-9XY]+[pq][0-9]"   # FAMILY_cytoband, e.g. ERVLE_1p36.32a

# Pick the plain .pos (header exactly WGT ID CHR P0 P1), not the CMC_wrapped_2.pos
# variant whose DIR column hard-codes an absolute author path.
cands <- list.files(pd, pattern = "\\.pos$", recursive = TRUE, full.names = TRUE)
cands <- cands[basename(cands) != basename(out_pos)]
pick <- NULL
for (f in cands) {
  h <- strsplit(trimws(readLines(f, 1)), "[ \t]+")[[1]]
  if (identical(h, c("WGT", "ID", "CHR", "P0", "P1"))) { pick <- f; break }
}
if (is.null(pick)) stop("No plain 'WGT ID CHR P0 P1' .pos found under ", pd)

src_wdir <- file.path(dirname(pick), "CMC_wrapped")
if (!dir.exists(src_wdir)) stop("Weights dir not found: ", src_wdir)

# Move the weights into a subdir named exactly <panel> and set WGT =
# '<panel>/<file>'. Both consumers then resolve: FUSION.assoc_test reads
# <weights_dir>/<WGT>, FeaturePred.V2.0.R reads <weights_dir>/<PANEL>/<basename>
# (FeaturePred.V2.0.R:229) - both require the files under <panel_dir>/<panel>/.
dest_wdir <- file.path(pd, opt$panel)
if (normalizePath(src_wdir) != normalizePath(dest_wdir, mustWork = FALSE)) {
  if (dir.exists(dest_wdir)) unlink(dest_wdir, recursive = TRUE)
  if (!file.rename(src_wdir, dest_wdir)) stop("Failed to move ", src_wdir, " -> ", dest_wdir)
}

pos <- fread(pick)
is_herv <- grepl(HERV_RE, pos$ID)
cat("prep_herv: ", sum(is_herv), " HERV of ", nrow(pos), " features (dropping ",
    sum(!is_herv), " gene features)\n", sep = "")
pos <- pos[is_herv]
if (nrow(pos) == 0L) stop("No HERV features matched ", HERV_RE)
pos[, WGT := file.path(opt$panel, basename(WGT))]
if (!file.exists(file.path(pd, pos$WGT[1])))
  stop("WGT path does not resolve under panel_dir after move: ", file.path(pd, pos$WGT[1]))

# Drop the non-HERV weight files to save space.
keep_bn <- basename(pos$WGT)
all_w <- list.files(dest_wdir, pattern = "\\.wgt\\.RDat$", full.names = TRUE)
unlink(all_w[!(basename(all_w) %in% keep_bn)])

pos[, PANEL := opt$panel]
pos[, N := opt$n]
out <- pos[, .(PANEL, WGT, ID, CHR, P0, P1, N)]
fwrite(out, out_pos, sep = " ", quote = FALSE)

# Tidy the now-empty nested extraction tree.
if (dir.exists(file.path(pd, "users"))) unlink(file.path(pd, "users"), recursive = TRUE)
cat("prep_herv: wrote", nrow(out), "HERV features (N=", opt$n, ") to ", out_pos, "\n")
