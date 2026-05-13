#!/usr/bin/env Rscript
# Extract per-trait counts of upstream signals and per-method gene lists:
#   - SNP-level: GW-sig clumped lead variants; COJO conditionally independent signals
#   - SuSiE finemapping: credible sets and unique mapped genes (L1, L10)
#   - TWAS-FUSION expression: FDR<0.05 AND COLOC_logical=TRUE, unique genes (any panel)
#   - SMR expression: FDR<0.05 AND HEIDI>0.01, unique genes (any panel)
#   - PWAS-FUSION protein: FDR<0.05 AND COLOC_logical=TRUE, unique genes (any panel)
#   - SMR protein: FDR<0.05 AND HEIDI>0.01, unique genes
#   - MAGMA gene-level: FDR<0.05
#   - Tissue enrichment: FDR<0.05 and retained after conditional analysis
# Also produce per-trait gene lists (one per method) for the cross-trait overlap
# analysis in cross_trait_overlap.R.

suppressMessages({ library(data.table) })

.script_dir <- local({
  args <- commandArgs(trailingOnly = FALSE)
  fa <- grep("^--file=", args, value = TRUE)
  if(length(fa) > 0) dirname(normalizePath(sub("^--file=", "", fa[1]))) else getwd()
})
source(file.path(.script_dir, "_paths.R"))
outdir <- tables_dir

x <- readRDS(rds_path)
traits <- setdiff(names(x), "configuration")

# ---- Per-trait scalar counts ----
fdr_th   <- 0.05
heidi_th <- 0.01  # SMR convention: p_HEIDI > 0.01 supports causal vs pleiotropy

# GW-significance threshold for clumping and COJO. The pipeline's
# .GW.clump.clean.csv / .GW.cojo.clean.csv files actually retain signals down
# to P<1e-5 (a broader clumping threshold). For a count of "GW-significant
# independent loci/signals" we apply the canonical P<5e-8 filter here.
gw_p_threshold <- 5e-8
count_clump <- function(d){
  if(is.null(d) || nrow(d) == 0) return(0L)
  sum(as.data.table(d)$P < gw_p_threshold, na.rm = TRUE)
}
count_cojo <- function(d){
  if(is.null(d) || nrow(d) == 0) return(0L)
  sum(as.data.table(d)$P < gw_p_threshold, na.rm = TRUE)
}
count_susie <- function(d){ if(is.null(d)) NA_integer_ else nrow(d) }
hc_fusion_exp <- function(d){
  if(is.null(d) || nrow(d) == 0) return(character(0))
  d <- as.data.table(d)
  unique(d[TWAS.P.FDR < fdr_th & COLOC_logical == TRUE, `Gene Symbol`])
}
hc_smr_exp <- function(d){
  if(is.null(d) || nrow(d) == 0) return(character(0))
  d <- as.data.table(d)
  unique(d[p_SMR.FDR < fdr_th & p_HEIDI > heidi_th, `Gene Symbol`])
}
hc_fusion_prot <- function(d){
  if(is.null(d) || nrow(d) == 0) return(character(0))
  d <- as.data.table(d)
  unique(d[pwas_all.P.FDR < fdr_th & COLOC_logical == TRUE, `Gene Symbol`])
}
hc_smr_prot <- function(d){
  if(is.null(d) || nrow(d) == 0) return(character(0))
  d <- as.data.table(d)
  unique(d[p_SMR.FDR < fdr_th & p_HEIDI > heidi_th, `Gene Symbol`])
}
hc_magma <- function(d){
  if(is.null(d) || nrow(d) == 0) return(character(0))
  d <- as.data.table(d)
  unique(d[P.FDR < fdr_th, ID])
}
hc_tissues <- function(spec){
  if(is.null(spec) || is.null(spec$res)) return(list(fdr=character(0), retained=character(0)))
  d <- as.data.table(spec$res)
  list(fdr      = d[P.FDR < fdr_th, Tissue],
       retained = if(is.null(spec$keep)) character(0) else spec$keep)
}

summary_rows <- list()
gene_lists   <- list()  # trait -> method -> character vector

for(tr in traits){
  tr_d <- x[[tr]]

  # SNP-level
  n_clump <- count_clump(tr_d$snp_assoc$clump)
  n_cojo  <- count_cojo(tr_d$snp_assoc$cojo)
  n_susie_L1  <- count_susie(tr_d$snp_assoc$susie$L1)
  n_susie_L10 <- count_susie(tr_d$snp_assoc$susie$L10)

  # Genes from finemapping
  genes_L1  <- tr_d$mol_assoc$finemap$L1
  genes_L10 <- tr_d$mol_assoc$finemap$L10
  if(is.null(genes_L1))  genes_L1  <- character(0)
  if(is.null(genes_L10)) genes_L10 <- character(0)

  # Molecular associations
  twas_genes      <- hc_fusion_exp(tr_d$mol_assoc$exp$fusion$res)
  smr_exp_genes   <- hc_smr_exp(tr_d$mol_assoc$exp$smr$results)
  pwas_genes      <- hc_fusion_prot(tr_d$mol_assoc$protein$fusion$results)
  smr_prot_genes  <- hc_smr_prot(tr_d$mol_assoc$protein$smr$results)
  magma_genes     <- hc_magma(tr_d$mol_assoc$magma)
  # Nearest-gene loci restricted to GW-significant clumps/COJO signals
  if(is.null(tr_d$snp_assoc$clump) || nrow(tr_d$snp_assoc$clump) == 0){
    nearest_clump <- character(0)
  } else {
    c_dt <- as.data.table(tr_d$snp_assoc$clump)
    nearest_clump <- unique(unlist(strsplit(c_dt[P < gw_p_threshold, NearestGene], ',', fixed = TRUE)))
    nearest_clump <- nearest_clump[nzchar(nearest_clump) & nearest_clump != "None"]
  }
  if(is.null(tr_d$snp_assoc$cojo) || nrow(tr_d$snp_assoc$cojo) == 0){
    nearest_cojo <- character(0)
  } else {
    c_dt <- as.data.table(tr_d$snp_assoc$cojo)
    nearest_cojo <- unique(unlist(strsplit(c_dt[P < gw_p_threshold, NearestGene], ',', fixed = TRUE)))
    nearest_cojo <- nearest_cojo[nzchar(nearest_cojo) & nearest_cojo != "None"]
  }

  # Tissue
  tlist <- hc_tissues(tr_d$tissue$specific)

  summary_rows[[tr]] <- data.table(
    Trait = tr,
    Clump_GWsig_loci             = n_clump,
    COJO_independent_signals     = n_cojo,
    SuSiE_L1_credible_sets       = n_susie_L1,
    SuSiE_L10_credible_sets      = n_susie_L10,
    SuSiE_L1_genes               = length(genes_L1),
    SuSiE_L10_genes              = length(genes_L10),
    TWAS_FUSION_hc_genes         = length(twas_genes),
    SMR_expression_hc_genes      = length(smr_exp_genes),
    PWAS_FUSION_hc_genes         = length(pwas_genes),
    SMR_protein_hc_genes         = length(smr_prot_genes),
    MAGMA_FDR_genes              = length(magma_genes),
    NearestGene_clump_loci       = length(nearest_clump),
    NearestGene_cojo_loci        = length(nearest_cojo),
    Tissue_FDR_enriched          = length(tlist$fdr),
    Tissue_retained_conditional  = length(tlist$retained)
  )

  gene_lists[[tr]] <- list(
    SuSiE_L1       = genes_L1,
    SuSiE_L10      = genes_L10,
    TWAS_FUSION    = twas_genes,
    SMR_expression = smr_exp_genes,
    PWAS_FUSION    = pwas_genes,
    SMR_protein    = smr_prot_genes,
    MAGMA          = magma_genes,
    Nearest_clump  = nearest_clump,
    Nearest_cojo   = nearest_cojo,
    Tissue_FDR        = tlist$fdr,
    Tissue_retained   = tlist$retained
  )
}

per_trait_summary <- rbindlist(summary_rows, use.names = TRUE)
fwrite(per_trait_summary, file.path(outdir, "upstream_summary_per_trait.tsv"), sep = "\t")

# ---- Long gene-table for the overlap analysis ----
long_gene <- rbindlist(lapply(traits, function(tr){
  L <- gene_lists[[tr]]
  rbindlist(lapply(names(L), function(method){
    if(length(L[[method]]) == 0) return(NULL)
    data.table(Trait = tr, Method = method, Gene = L[[method]])
  }))
}))
fwrite(long_gene, file.path(outdir, "upstream_high_confidence_genes_long.tsv"), sep = "\t")

# ---- Per-method gene-trait wide matrix (gene x trait, value = 1 if in list) ----
for(method in c("SuSiE_L10","TWAS_FUSION","SMR_expression","PWAS_FUSION","SMR_protein","MAGMA")){
  sub <- long_gene[Method == method]
  if(nrow(sub) == 0) next
  wide <- dcast(sub, Gene ~ Trait, value.var = "Method", fun.aggregate = length)
  wide[, N_traits := rowSums(.SD), .SDcols = traits[traits %in% names(wide)]]
  setorder(wide, -N_traits)
  fwrite(wide, file.path(outdir, paste0("genes_x_trait_", method, ".tsv")), sep = "\t")
}

# ---- Cross-trait overlap counts: how many high-confidence genes are shared per trait pair, per method ----
overlap_pairs <- list()
for(method in c("SuSiE_L10","TWAS_FUSION","SMR_expression","PWAS_FUSION","SMR_protein","MAGMA")){
  for(i in seq_along(traits)){
    for(j in seq_along(traits)){
      g_i <- gene_lists[[traits[i]]][[method]]
      g_j <- gene_lists[[traits[j]]][[method]]
      overlap <- length(intersect(g_i, g_j))
      overlap_pairs[[length(overlap_pairs)+1]] <- data.table(
        Method = method, TraitA = traits[i], TraitB = traits[j],
        N_A = length(g_i), N_B = length(g_j),
        N_shared = overlap,
        Jaccard = if(length(union(g_i, g_j)) > 0) overlap / length(union(g_i, g_j)) else 0
      )
    }
  }
}
overlap_dt <- rbindlist(overlap_pairs)
fwrite(overlap_dt, file.path(outdir, "cross_trait_overlap_genes.tsv"), sep = "\t")

# ---- Multi-trait high-confidence genes (>=2 traits in any one method) ----
recur <- long_gene[, .(N_traits = uniqueN(Trait), Traits = paste(sort(unique(Trait)), collapse=",")), by = .(Method, Gene)]
recur <- recur[N_traits >= 2][order(Method, -N_traits, Gene)]
fwrite(recur, file.path(outdir, "cross_trait_recurring_genes.tsv"), sep = "\t")

# ---- Tissue enrichment cross-trait ----
tissue_long <- rbindlist(lapply(traits, function(tr){
  s <- x[[tr]]$tissue$specific$res
  if(is.null(s) || nrow(s) == 0) return(NULL)
  d <- as.data.table(copy(s))
  d[, Trait := tr]
  keep_vec <- if(is.null(x[[tr]]$tissue$specific$keep)) character(0) else x[[tr]]$tissue$specific$keep
  d[, Retained := Tissue %in% keep_vec]
  d[, FDR_sig := P.FDR < fdr_th]
  d
}))
fwrite(tissue_long, file.path(outdir, "tissue_enrichment_all_traits.tsv"), sep = "\t")

# Tissue x trait wide matrix
tissue_wide <- dcast(tissue_long[, .(Trait, Tissue, P.FDR)],
                     Tissue ~ Trait, value.var = "P.FDR", fill = NA)
fwrite(tissue_wide, file.path(outdir, "tissue_fdr_by_trait_matrix.tsv"), sep = "\t")

cat("=== Per-trait upstream summary ===\n")
print(per_trait_summary)
cat("\n=== Cross-trait recurring genes by method (top 20 per method) ===\n")
for(m in unique(recur$Method)){
  cat("\n--- ", m, " ---\n", sep="")
  print(head(recur[Method == m], 15))
}
