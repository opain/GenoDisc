#!/usr/bin/Rscript

# Per-gene evidence behind the DrugTargetor drug- and ATC-class enrichments, for the
# Shiny drill-down. For one (gwas, panel) it emits, per drug and per ATC class
# (L2/L3/L4), the member genes that carry the signal: signed membership, the gene's
# TWAS Z/P, the contribution (membership x TWAS.Z), and the evidence provenance
# (which source DB(s) and activity_type(s) support the drug-gene link). This is what
# lets a user see, e.g., that the lipid class "matches" CAD because ezetimibe's
# DSIGDB signature UP-regulates PCSK9 (a strong CAD risk gene), while the true
# therapeutic target NPC1L1 isn't TWAS-modelled.
#
# Built sparsely from the RAW DrugTargetor database (not the dense .prop), so the
# same pass yields both membership (identical convention to
# format_drug_targetor_atc_for_twas_gsea.R) and provenance.

suppressMessages(library(optparse))
option_list <- list(
  make_option("--gwas", type="character", default=NA, help="GWAS ID [required]"),
  make_option("--panel", type="character", default=NA, help="TWAS panel/weight [required]"),
  make_option("--resdir", type="character", default="resources", help="Resources dir"),
  make_option("--outdir", type="character", default=NA, help="Pipeline outdir [required]"),
  make_option("--pipeline_dir", type="character", default=NA, help="Pipeline dir"),
  make_option("--topk", type="integer", default=50L, help="Max member genes kept per drug/class (by |contribution|)"),
  make_option("--levels", type="character", default="l2,l3,l4", help="ATC levels to emit")
)
opt <- parse_args(OptionParser(option_list=option_list))
suppressMessages(library(data.table))

level_width <- c(l2=3L, l3=4L, l4=5L)
levels <- trimws(strsplit(opt$levels, ",", fixed=TRUE)[[1]])
if(!all(levels %in% names(level_width))) stop("--levels must be a subset of l2,l3,l4")

## --- inputs -----------------------------------------------------------------
universe <- unique(fread(file.path(opt$resdir,"data/magma/NCBI37.3.gene.loc"), header=FALSE)$V6)

# Modelled genes for this panel: per-gene TWAS Z/P (this is the only signal source)
tw <- fread(file.path(opt$outdir,"results",opt$gwas,"twas",
                      paste0(opt$gwas,"_twas_",opt$panel,"_GW_clean.txt.gz")),
            select=c("external_gene_name","TWAS.Z","TWAS.P"))
setnames(tw, c("gene","TWAS.Z","TWAS.P"))
tw <- unique(tw[!is.na(TWAS.Z)], by="gene")

# Raw assertions + sign, same mapping as format_drug_targetor_atc_for_twas_gsea.R
db <- fread(file.path(opt$resdir,"data/drug_targetor/wholedatabase_for_targetor"),
            sep="\t", header=TRUE)
db <- db[gene %in% universe]
db[, sgn := fifelse(activity_type %in% c("DECREASED_EXPRESSION","NEGATIVE_RESPONSE","OPPOSITE_RESPONSE"), -1L,
             fifelse(activity_type %in% c("INCREASED_EXPRESSION","POSITIVE_RESPONSE"), 1L, 0L))]
db[, drug_name := sub(".*NAME:([^|]+)\\|.*","\\1", atc)]
db[, atc7 := sub("^ATC:([^|]+)\\|.*","\\1", atc)]   # composite (may be comma-joined)

topk <- function(dt, by) dt[order(-abs(contribution))][, head(.SD, opt$topk), by=by]

## --- DRUG-level evidence ----------------------------------------------------
# one signed membership per (drug, gene): -1 precedence; keep provenance
dg <- db[sgn!=0L, .(membership = if(any(sgn==-1L)) -1L else 1L,
                    activity = paste(sort(unique(activity_type)), collapse="; "),
                    source   = paste(sort(unique(source)),        collapse="; ")),
         by=.(atc, drug_name, gene)]
drug_summary <- dg[, .(n_target=.N), by=.(atc, drug_name)]
dge <- merge(dg, tw, by="gene")                       # modelled members only
dge[, contribution := membership * TWAS.Z]
drug_summary <- merge(drug_summary,
                      dge[, .(n_modelled=.N), by=.(atc,drug_name)], by=c("atc","drug_name"), all.x=TRUE)
drug_summary[is.na(n_modelled), n_modelled := 0L]
drug_evidence <- topk(dge, c("atc","drug_name"))[
  , .(atc, drug_name, gene, membership, activity, source,
      TWAS.Z=round(TWAS.Z,3), TWAS.P, contribution=round(contribution,3))]

## --- CLASS-level evidence (L2/L3/L4) ----------------------------------------
# explode composite ATC codes -> per-drug code list
udrug <- unique(db[, .(atc)])
codes <- strsplit(sub("^ATC:([^|]+)\\|.*","\\1", udrug$atc), ",", fixed=TRUE)
drug_codes <- data.table(atc=rep(udrug$atc, lengths(codes)), code=unlist(codes))

class_evidence <- list(); class_summary <- list()
for(lv in levels){
  k <- level_width[[lv]]
  dc <- unique(drug_codes[, .(atc, class=substr(code,1,k))])
  # (class, gene): mean signed membership over drugs in class targeting the gene (matches .prop),
  # with aggregated provenance and supporting-drug count
  m <- merge(dc, dg, by="atc", allow.cartesian=TRUE)         # atc, class, drug_name, gene, membership, activity, source
  cg <- m[, .(membership = mean(membership),
              n_drugs    = uniqueN(atc),
              activity   = paste(sort(unique(unlist(strsplit(activity,"; ")))), collapse="; "),
              source     = paste(sort(unique(unlist(strsplit(source,  "; ")))), collapse="; ")),
          by=.(class, gene)]
  csum <- cg[, .(level=lv, n_target=.N), by=.(code=class)]
  cge <- merge(cg, tw, by="gene")
  cge[, contribution := membership * TWAS.Z]
  csum <- merge(csum, cge[, .(n_modelled=.N), by=.(code=class)], by="code", all.x=TRUE)
  csum[is.na(n_modelled), n_modelled := 0L]
  ce <- cge[order(-abs(contribution))][, head(.SD, opt$topk), by=class][
    , .(level=lv, code=class, gene, membership=round(membership,3), n_drugs, activity, source,
        TWAS.Z=round(TWAS.Z,3), TWAS.P, contribution=round(contribution,3))]
  class_evidence[[lv]] <- ce
  class_summary[[lv]]  <- csum
}

out <- list(
  panel = opt$panel,
  drug_evidence  = drug_evidence,
  drug_summary   = drug_summary,
  class_evidence = rbindlist(class_evidence),
  class_summary  = rbindlist(class_summary)
)
outfile <- file.path(opt$outdir,"results",opt$gwas,"twas","drugtargetor",
                     paste0("twas_gsea_drugtargetor_evidence_",opt$panel,".rds"))
saveRDS(out, outfile)
cat(sprintf("evidence[%s]: %d drug rows, %d class rows -> %s\n",
            opt$panel, nrow(drug_evidence), nrow(out$class_evidence), outfile))
