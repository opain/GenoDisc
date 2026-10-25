#!/usr/bin/Rscript
# Prepare the fetal-brain gene-level FUSION panel (Bray & Pain 2020, figshare
# 11637036) after download/extraction into the panel directory. Writes the
# canonical <panel>.pos (PANEL WGT ID CHR P0 P1 N) that run_twas reads with
# --weights_dir <panel_dir>. IDs are already unversioned Ensembl gene IDs, so
# (unlike insert_n_nongtex.R) there is NO biomaRt symbol->Ensembl remap here.
suppressMessages(library(optparse)); suppressMessages(library(data.table))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--panel_dir", type = "character", help = "Panel directory (contains the extracted archive)"),
  make_option("--panel", type = "character", default = "FETAL.BRAIN.RNASEQ"),
  make_option("--n", type = "integer", default = 120L, help = "Sample size (fetal brains)")
)))
pd <- opt$panel_dir
out_pos <- file.path(pd, paste0(opt$panel, ".pos"))

# Locate the single source .pos produced by extraction (ignore our own output).
cands <- list.files(pd, pattern = "\\.pos$", recursive = TRUE, full.names = TRUE)
cands <- cands[basename(cands) != basename(out_pos)]
if (length(cands) != 1L) stop("Expected exactly one source .pos under ", pd, ", found ", length(cands))
src <- cands[1]
pos <- fread(src)
stopifnot(all(c("WGT", "ID", "CHR", "P0", "P1") %in% names(pos)))

# WGT is relative to the directory containing the source .pos. The archive nests
# the weights under a subdir whose name collides with the outer dir, so rather
# than flatten, move the referenced weight files into a clean wgt/ subdir under
# the panel dir and rewrite WGT to 'wgt/<file>' (resolves with --weights_dir
# <panel_dir>). Basenames are unique (Fetal_ENSG...).
wdir <- dirname(src)
if (!file.exists(file.path(wdir, pos$WGT[1])))
  stop("WGT path does not resolve under ", wdir, ": ", pos$WGT[1])
wgt_dir <- file.path(pd, "wgt")
dir.create(wgt_dir, showWarnings = FALSE)
new_wgt <- file.path("wgt", basename(pos$WGT))
moved <- file.rename(file.path(wdir, pos$WGT), file.path(pd, new_wgt))
if (!all(moved)) stop("Failed to move ", sum(!moved), " weight files into ", wgt_dir)
pos[, WGT := new_wgt]
if (!file.exists(file.path(pd, pos$WGT[1])))
  stop("WGT path does not resolve under panel_dir after move: ", file.path(pd, pos$WGT[1]))
# Drop the leftover extracted archive dir(s), keeping only wgt/.
for (d in list.dirs(pd, recursive = FALSE))
  if (normalizePath(d) != normalizePath(wgt_dir)) unlink(d, recursive = TRUE)

# Features must be Ensembl gene IDs (combine_twas.R merges on ensembl_gene_id).
bad <- !grepl("^ENSG[0-9]+", pos$ID)
if (any(bad)) stop(sum(bad), " fetal feature IDs are not Ensembl gene IDs (e.g. ",
                   paste(head(pos$ID[bad], 3), collapse = ", "), ")")
pos[, ID := sub("\\..*$", "", ID)]   # defensive: strip any version suffix

pos[, PANEL := opt$panel]
pos[, N := opt$n]
out <- pos[, .(PANEL, WGT, ID, CHR, P0, P1, N)]
fwrite(out, out_pos, sep = " ", quote = FALSE)
cat("prep_fetal_brain: wrote", nrow(out), "gene features (N=", opt$n, ") to ", out_pos, "\n")
