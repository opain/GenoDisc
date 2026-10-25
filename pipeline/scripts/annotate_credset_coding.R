#!/usr/bin/Rscript
# Annotate SuSiE credible-set variants (GRCh37/hg19) with coding consequence.
#  - AlphaMissense (missense pathogenicity + class) via a tabix join on the hg19 table
#  - bcftools csq (stop/frameshift/splice/inframe/start-lost/synonymous ...) via an
#    Ensembl GRCh37 GFF3 + FASTA
# REF/ALT are resolved from the reference FASTA (our A1/A2 are GWAS effect/other
# alleles, not genomic ref/alt). Produces a per-variant long table and a per-cs
# summary that flags credible sets carrying a meaningful share of PIP on
# protein-altering variants.
#
# Inputs come from extract_credset_variants.R; shells out to samtools/tabix/bcftools
# (all in envs/coding_annot.yaml). Writes header-only outputs when there are no
# credible-set variants, so packaging always finds the files.

suppressMessages(library("optparse"))
option_list <- list(
  make_option("--gwas", default = NA, type = "character"),
  make_option("--config_file", default = NA, type = "character"),
  make_option("--pipeline_dir", default = NA, type = "character"),
  make_option("--alphamissense", default = NA, type = "character", help = "bgzipped+tabixed AlphaMissense_hg19"),
  make_option("--gff3", default = NA, type = "character", help = "Ensembl GRCh37 GFF3 (bgzipped) for bcftools csq"),
  make_option("--fasta", default = NA, type = "character", help = "GRCh37 reference FASTA (faidx-indexed)"),
  make_option("--pip_thresh", default = 0.1, type = "double")
)
opt <- parse_args(OptionParser(option_list = option_list))
options(pipeline_dir = opt$pipeline_dir)
suppressMessages({ library(data.table); library(stringr) })
source(file.path(opt$pipeline_dir, "scripts", "functions", "utils_functions.R"))

outdir <- read_param(config = opt$config_file, param = "outdir", return_obj = FALSE)
fm_dir <- paste0(outdir, "/results/", opt$gwas, "/finemap/")
long_path <- paste0(fm_dir, opt$gwas, ".GW.credset_coding.csv")
summ_path <- paste0(fm_dir, opt$gwas, ".GW.credset_coding_summary.csv")

long_cols <- c("locus", "cs", "L", "SNP", "CHR", "BP", "A1", "A2", "REF", "ALT",
               "PIP", "consequence", "gene", "am_pathogenicity", "am_class", "protein_altering")
summ_cols <- c("locus", "cs", "L", "n_variants", "pip_total", "pip_protein_altering",
               "frac_protein_altering", "top_coding_snp", "top_coding_gene",
               "top_coding_consequence", "top_coding_am_class", "top_coding_pip", "flag_protein_altering")
write_empty <- function() {
  fwrite(setNames(data.table(matrix(character(0), ncol = length(long_cols))), long_cols), long_path)
  fwrite(setNames(data.table(matrix(character(0), ncol = length(summ_cols))), summ_cols), summ_path)
}

cv <- tryCatch(fread(paste0(fm_dir, opt$gwas, ".GW.credset_variants.csv")), error = function(e) NULL)
if (is.null(cv) || nrow(cv) == 0) { write_empty(); quit(save = "no") }
# CHR is read as integer by fread; keep it character so joins against the
# (character) AlphaMissense/FASTA contig keys match.
cv[, CHR := as.character(CHR)]

# consequences counted as protein-altering (coding-sequence or splice-site changing)
PA_CSQ <- c("missense", "stop_gained", "stop_lost", "start_lost", "frameshift",
            "inframe_deletion", "inframe_insertion", "splice_acceptor", "splice_donor")

run <- function(cmd, args) system2(cmd, args, stdout = TRUE, stderr = FALSE)

# --- resolve REF/ALT from the reference FASTA (SNVs only) --------------------
cv[, is_snv := nchar(A1) == 1 & nchar(A2) == 1 & A1 %in% c("A","C","G","T") & A2 %in% c("A","C","G","T")]
uniq <- unique(cv[is_snv == TRUE, .(CHR, BP)])
uniq[, region := paste0(CHR, ":", BP, "-", BP)]
reg_file <- tempfile()
writeLines(uniq$region, reg_file)
# samtools faidx returns one FASTA record per region: ">CHR:BP-BP\n<base>"
fa <- run("samtools", c("faidx", shQuote(opt$fasta), "-r", shQuote(reg_file)))
ref_base <- character(0); hdr <- NULL
for (ln in fa) {
  if (startsWith(ln, ">")) hdr <- sub("^>", "", ln)
  else if (!is.null(hdr)) { ref_base[hdr] <- toupper(ln); hdr <- NULL }
}
uniq[, REF := unname(ref_base[region])]
cv <- merge(cv, uniq[, .(CHR, BP, REF)], by = c("CHR", "BP"), all.x = TRUE, sort = FALSE)
# ALT = the GWAS allele that is not the reference base; ambiguous if neither/both match
cv[, ALT := NA_character_]
cv[is_snv == TRUE & A1 == REF & A2 != REF, ALT := A2]
cv[is_snv == TRUE & A2 == REF & A1 != REF, ALT := A1]

# --- AlphaMissense join (missense pathogenicity) ----------------------------
cv[, am_pathogenicity := NA_real_][, am_class := NA_character_]
am_targets <- unique(cv[!is.na(ALT), .(CHR, BP, REF, ALT)])
if (nrow(am_targets) > 0 && !is.na(opt$alphamissense)) {
  # tabix -R takes a BED file (0-based, half-open), not region strings.
  amreg <- tempfile(fileext = ".bed")
  writeLines(paste0("chr", am_targets$CHR, "\t", am_targets$BP - 1L, "\t", am_targets$BP), amreg)
  am <- run("tabix", c("-R", shQuote(amreg), shQuote(opt$alphamissense)))
  if (length(am) > 0) {
    amdt <- fread(text = paste(am, collapse = "\n"), header = FALSE)
    # cols: CHROM POS REF ALT genome uniprot transcript protein_variant am_pathogenicity am_class
    setnames(amdt, 1:10, c("CHROM","POS","aREF","aALT","genome","uniprot","transcript","pv","am_path","am_cls"))
    amdt[, CHR := sub("^chr", "", CHROM)][, BP := as.integer(POS)]
    amdt <- amdt[, .(am_path = am_path[1], am_cls = am_cls[1]),
                 by = .(CHR, BP, aREF, aALT)]
    cv <- merge(cv, amdt, by.x = c("CHR","BP","REF","ALT"), by.y = c("CHR","BP","aREF","aALT"),
                all.x = TRUE, sort = FALSE)
    cv[!is.na(am_path), am_pathogenicity := as.numeric(am_path)]
    cv[!is.na(am_cls),  am_class := am_cls]
    cv[, c("am_path","am_cls") := NULL]
  }
}

# --- bcftools csq (all protein-altering consequences) -----------------------
cv[, consequence := NA_character_][, gene := NA_character_]
vcf_dt <- cv[!is.na(ALT), .(CHR, BP, SNP, REF, ALT)]
if (nrow(vcf_dt) > 0 && !is.na(opt$gff3) && !is.na(opt$fasta)) {
  setorder(vcf_dt, CHR, BP)
  vcf <- tempfile(fileext = ".vcf")
  con <- file(vcf, "w")
  writeLines(c("##fileformat=VCFv4.2", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"), con)
  writeLines(sprintf("%s\t%d\t%s\t%s\t%s\t.\t.\t.", vcf_dt$CHR, vcf_dt$BP, vcf_dt$SNP, vcf_dt$REF, vcf_dt$ALT), con)
  close(con)
  csq <- run("bcftools", c("csq", "-f", shQuote(opt$fasta), "-g", shQuote(opt$gff3),
                           "-p", "a", "-l", "-Ov", shQuote(vcf)))
  csq <- csq[!startsWith(csq, "#")]
  if (length(csq) > 0) {
    cd <- fread(text = paste(csq, collapse = "\n"), header = FALSE, sep = "\t")
    # INFO is the last column; BCSQ=consequence|gene|transcript|biotype|strand|aa|...
    cd[, info := cd[[ncol(cd)]]]
    cd[, bcsq := sub(".*BCSQ=", "", info)]
    cd[, bcsq := sub(";.*", "", bcsq)]
    cd[, first := sub(",.*", "", bcsq)]         # most-severe consequence is listed first
    parts <- tstrsplit(cd$first, "|", fixed = TRUE)
    cd[, csq_cons := sub("[@*].*", "", parts[[1]])]  # strip compound/downstream markers
    cd[, csq_gene := if (length(parts) >= 2) parts[[2]] else NA_character_]
    cd[, SNP := cd[[3]]]
    # csq emits one record per transcript/consequence; keep the first (most-severe)
    # per variant so the merge doesn't duplicate credible-set rows.
    cd <- cd[!duplicated(SNP)]
    cv <- merge(cv, cd[, .(SNP, csq_cons, csq_gene)], by = "SNP", all.x = TRUE, sort = FALSE)
    cv[!is.na(csq_cons), consequence := csq_cons]
    cv[!is.na(csq_gene) & nzchar(csq_gene), gene := csq_gene]
    cv[, c("csq_cons","csq_gene") := NULL]
  }
}

# AlphaMissense is authoritative for the missense label/gene where present.
cv[!is.na(am_class) & (is.na(consequence) | consequence == "synonymous"), consequence := "missense"]
cv[is_snv == FALSE, consequence := fifelse(is.na(consequence), "not_annotated_non_snv", consequence)]
cv[, protein_altering := !is.na(consequence) & consequence %in% PA_CSQ]

# --- write long table -------------------------------------------------------
for (cc in setdiff(long_cols, names(cv))) cv[, (cc) := NA]
setorder(cv, locus, cs, -PIP)
fwrite(cv[, ..long_cols], long_path)

# --- per-cs summary ---------------------------------------------------------
cv[, pa_pip := fifelse(protein_altering, PIP, 0)]
summ <- cv[, {
  top <- .SD[protein_altering == TRUE][which.max(PIP)]
  .(L = L[1],
    n_variants = .N,
    pip_total = sum(PIP, na.rm = TRUE),
    pip_protein_altering = sum(pa_pip, na.rm = TRUE),
    frac_protein_altering = ifelse(sum(PIP, na.rm = TRUE) > 0, sum(pa_pip, na.rm = TRUE) / sum(PIP, na.rm = TRUE), 0),
    top_coding_snp = if (nrow(top)) top$SNP else NA_character_,
    top_coding_gene = if (nrow(top)) top$gene else NA_character_,
    top_coding_consequence = if (nrow(top)) top$consequence else NA_character_,
    top_coding_am_class = if (nrow(top)) top$am_class else NA_character_,
    top_coding_pip = if (nrow(top)) top$PIP else NA_real_)
}, by = .(locus, cs)]
summ[, flag_protein_altering := pip_protein_altering >= opt$pip_thresh]
setcolorder(summ, summ_cols)
setorder(summ, -pip_protein_altering)
fwrite(summ, summ_path)
