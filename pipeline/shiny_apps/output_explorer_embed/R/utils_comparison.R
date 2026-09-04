# Cross-GWAS comparison data layer
#
# All views in comparison mode read a single long data.table produced by
# build_comparison_long(). Filtering and pivoting is done on that tibble,
# never on the raw package. Column semantics:
#
#   gwas          character   GWAS name (matches gd_gwas())
#   entity_type   character   "tissue" | "atc" | (later) "gene" | "locus" | "drug"
#   entity_id     character   canonical identifier (Tissue name, ATC code, gene symbol, ...)
#   entity_label  character   human-readable label (may equal entity_id)
#   method        character   scoring method (e.g. "MAGMA-tissue", "TWAS-GSEA-ATC")
#   panel         character   panel name where applicable, NA otherwise
#   statistic     numeric     effect / z-like statistic where available
#   se            numeric     standard error where available
#   p             numeric     nominal p-value
#   fdr           numeric     FDR-adjusted p-value
#   evidence      logical     module-specific formal-evidence flag (e.g. tissue "retained")
#   direction     character   "matches" | "opposes" | other (from Direction column)
#   reversal_z    numeric     drug-signature reversal z (TWAS-GSEA only)
#   n_units       integer     count within the entity (N Gene, N Drugs, ...)

.gd_long_cols <- c(
  "gwas", "entity_type", "entity_id", "entity_label",
  "method", "panel", "statistic", "se", "p", "fdr",
  "evidence", "direction", "reversal_z", "n_units"
)

.gd_empty_long <- function() {
  data.table::data.table(
    gwas         = character(0),
    entity_type  = character(0),
    entity_id    = character(0),
    entity_label = character(0),
    method       = character(0),
    panel        = character(0),
    statistic    = numeric(0),
    se           = numeric(0),
    p            = numeric(0),
    fdr          = numeric(0),
    evidence     = logical(0),
    direction    = character(0),
    reversal_z   = numeric(0),
    n_units      = integer(0)
  )
}

.gd_long_row <- function(gwas, entity_type, entity_id, method,
                        entity_label = entity_id, panel = NA_character_,
                        statistic = NA_real_, se = NA_real_,
                        p = NA_real_, fdr = NA_real_,
                        evidence = NA, direction = NA_character_,
                        reversal_z = NA_real_, n_units = NA_integer_) {
  data.table::data.table(
    gwas         = as.character(gwas),
    entity_type  = as.character(entity_type),
    entity_id    = as.character(entity_id),
    entity_label = as.character(entity_label),
    method       = as.character(method),
    panel        = as.character(panel),
    statistic    = as.numeric(statistic),
    se           = as.numeric(se),
    p            = as.numeric(p),
    fdr          = as.numeric(fdr),
    evidence     = as.logical(evidence),
    direction    = as.character(direction),
    reversal_z   = as.numeric(reversal_z),
    n_units      = as.integer(n_units)
  )
}

.gd_long_tissue <- function(gd, gwas) {
  spec <- safe_access(gd_read(gd, gwas, "tissue"), "specific")
  if (is.null(spec) || is.null(spec$res) || nrow(spec$res) == 0) return(NULL)
  d <- as.data.frame(spec$res)
  keep <- if (is.null(spec$keep)) character(0) else spec$keep
  n_units <- suppressWarnings(as.integer(d[["N Gene"]]))
  .gd_long_row(
    gwas         = gwas,
    entity_type  = "tissue",
    entity_id    = d$Tissue,
    method       = "MAGMA-tissue",
    entity_label = d$Tissue,
    panel        = NA_character_,
    statistic    = suppressWarnings(as.numeric(d$BETA)),
    se           = suppressWarnings(as.numeric(d$SE)),
    p            = suppressWarnings(as.numeric(d$P)),
    fdr          = suppressWarnings(as.numeric(d$P.FDR)),
    evidence     = d$Tissue %in% keep,
    n_units      = n_units
  )
}

.gd_long_atc_magma <- function(gd, gwas) {
  m <- safe_access(gd_read(gd, gwas, "tx/atc"), "magma")
  if (is.null(m) || nrow(m) == 0) return(NULL)
  m <- as.data.frame(m)
  .gd_long_row(
    gwas         = gwas,
    entity_type  = "atc",
    entity_id    = m[["ATC Code"]],
    method       = "MAGMA-ATC",
    entity_label = paste0(m[["ATC Code"]], ": ", m[["ATC Description"]]),
    panel        = NA_character_,
    p            = suppressWarnings(as.numeric(m$P)),
    fdr          = suppressWarnings(as.numeric(m$P.FDR)),
    n_units      = suppressWarnings(as.integer(m[["N Drugs"]]))
  )
}

# The NearestGene column in snp_assoc$clump / $cojo is a comma-separated list
# of neighbouring gene annotations, e.g.:
#   "C1orf222, TMEM52 (+24.56kb), CALML6 (+26.53kb)"
# The first entry has no distance suffix (it is inside / directly beside the
# lead SNP). Extract just that first gene as the canonical locus label so
# rows match across GWAS.
.gd_first_nearest_gene <- function(x) {
  first <- sub(",.*$", "", x)
  first <- sub("\\s*\\([-+][0-9.]+\\s*kb\\)\\s*$", "", first)
  trimws(first)
}

.gd_long_locus_clump <- function(gd, gwas) {
  d <- safe_access(gd_read(gd, gwas, "snp_assoc"), "clump")
  if (is.null(d) || nrow(d) == 0) return(NULL)
  d <- as.data.frame(d)
  gene <- .gd_first_nearest_gene(d$NearestGene)
  missing <- is.na(gene) | !nzchar(gene) | tolower(gene) == "none"
  gene[missing] <- d$SNP[missing]
  # Drop duplicate (gene) hits per GWAS: keep the row with the smallest p.
  ord <- order(suppressWarnings(as.numeric(d$P)), na.last = TRUE)
  d <- d[ord, , drop = FALSE]
  gene <- gene[ord]
  dup <- duplicated(gene)
  d <- d[!dup, , drop = FALSE]; gene <- gene[!dup]
  .gd_long_row(
    gwas         = gwas,
    entity_type  = "locus",
    entity_id    = gene,
    method       = "clump",
    entity_label = gene,
    panel        = NA_character_,
    statistic    = suppressWarnings(as.numeric(d$BETA)),
    se           = suppressWarnings(as.numeric(d$SE)),
    p            = suppressWarnings(as.numeric(d$P)),
    fdr          = NA_real_,
    n_units      = 1L
  )
}

.gd_long_locus_cojo <- function(gd, gwas) {
  d <- safe_access(gd_read(gd, gwas, "snp_assoc"), "cojo")
  if (is.null(d) || nrow(d) == 0) return(NULL)
  d <- as.data.frame(d)
  gene <- .gd_first_nearest_gene(d$NearestGene)
  missing <- is.na(gene) | !nzchar(gene) | tolower(gene) == "none"
  gene[missing] <- d$SNP[missing]
  # COJO reports both marginal (P) and conditional-joint (pJ) p-values.
  # The pJ is what identifies signals that stay significant after
  # conditioning; that is the value users want to compare across traits.
  pj <- suppressWarnings(as.numeric(d$pJ))
  ord <- order(pj, na.last = TRUE)
  d <- d[ord, , drop = FALSE]
  gene <- gene[ord]; pj <- pj[ord]
  dup <- duplicated(gene)
  d <- d[!dup, , drop = FALSE]; gene <- gene[!dup]; pj <- pj[!dup]
  .gd_long_row(
    gwas         = gwas,
    entity_type  = "locus",
    entity_id    = gene,
    method       = "COJO",
    entity_label = gene,
    panel        = NA_character_,
    statistic    = suppressWarnings(as.numeric(d$BETA)),
    se           = suppressWarnings(as.numeric(d$SE)),
    p            = pj,
    fdr          = NA_real_,
    n_units      = 1L
  )
}

.gd_long_gene_magma <- function(gd, gwas) {
  m <- gd_read(gd, gwas, "mol_assoc/magma")
  if (is.null(m) || nrow(m) == 0) return(NULL)
  m <- as.data.frame(m)
  .gd_long_row(
    gwas         = gwas,
    entity_type  = "gene",
    entity_id    = m$ID,
    method       = "MAGMA-gene",
    entity_label = m$ID,
    panel        = NA_character_,
    p            = suppressWarnings(as.numeric(m$P)),
    fdr          = suppressWarnings(as.numeric(m$P.FDR))
  )
}

.gd_long_gene_twas_fusion <- function(gd, gwas) {
  r <- safe_access(gd_read(gd, gwas, "mol_assoc/exp/fusion"), "res")
  if (is.null(r) || nrow(r) == 0) return(NULL)
  r <- as.data.frame(r)
  .gd_long_row(
    gwas         = gwas,
    entity_type  = "gene",
    entity_id    = r[["Gene Symbol"]],
    method       = "TWAS-FUSION",
    entity_label = r[["Gene Symbol"]],
    panel        = as.character(r$PANEL),
    statistic    = suppressWarnings(as.numeric(r$TWAS.Z)),
    p            = suppressWarnings(as.numeric(r$TWAS.P)),
    fdr          = suppressWarnings(as.numeric(r$TWAS.P.FDR)),
    evidence     = as.logical(r$COLOC_logical)
  )
}

.gd_long_gene_smr_exp <- function(gd, gwas) {
  r <- safe_access(gd_read(gd, gwas, "mol_assoc/exp/smr"), "results")
  if (is.null(r) || nrow(r) == 0) return(NULL)
  r <- as.data.frame(r)
  # HEIDI evidence: p_HEIDI > 0.05 supports colocalisation (no evidence of
  # heterogeneity). NA HEIDI cannot be interpreted as supportive.
  ev <- !is.na(r$p_HEIDI) & suppressWarnings(as.numeric(r$p_HEIDI)) > 0.05
  ev[is.na(r$p_HEIDI)] <- NA
  gene_id <- r[["Gene Symbol"]]
  gene_id[is.na(gene_id)] <- r[["Ensembl ID"]][is.na(gene_id)]
  .gd_long_row(
    gwas         = gwas,
    entity_type  = "gene",
    entity_id    = gene_id,
    method       = "SMR-expression",
    entity_label = gene_id,
    panel        = as.character(r$PANEL),
    statistic    = suppressWarnings(as.numeric(r$b_SMR)),
    se           = suppressWarnings(as.numeric(r$se_SMR)),
    p            = suppressWarnings(as.numeric(r$p_SMR)),
    fdr          = suppressWarnings(as.numeric(r$p_SMR.FDR)),
    evidence     = ev
  )
}

.gd_long_gene_pwas_fusion <- function(gd, gwas) {
  r <- safe_access(gd_read(gd, gwas, "mol_assoc/protein/fusion"), "results")
  if (is.null(r) || nrow(r) == 0) return(NULL)
  r <- as.data.frame(r)
  .gd_long_row(
    gwas         = gwas,
    entity_type  = "gene",
    entity_id    = r[["Gene Symbol"]],
    method       = "PWAS-FUSION",
    entity_label = r[["Gene Symbol"]],
    panel        = as.character(r$PANEL),
    statistic    = suppressWarnings(as.numeric(r$pwas_all.Z)),
    p            = suppressWarnings(as.numeric(r$pwas_all.P)),
    fdr          = suppressWarnings(as.numeric(r$pwas_all.P.FDR)),
    evidence     = as.logical(r$COLOC_logical)
  )
}

.gd_long_gene_smr_prot <- function(gd, gwas) {
  r <- safe_access(gd_read(gd, gwas, "mol_assoc/protein/smr"), "results")
  if (is.null(r) || nrow(r) == 0) return(NULL)
  r <- as.data.frame(r)
  ev <- !is.na(r$p_HEIDI) & suppressWarnings(as.numeric(r$p_HEIDI)) > 0.05
  ev[is.na(r$p_HEIDI)] <- NA
  gene_id <- r[["Gene Symbol"]]
  gene_id[is.na(gene_id)] <- r[["Ensembl ID"]][is.na(gene_id)]
  .gd_long_row(
    gwas         = gwas,
    entity_type  = "gene",
    entity_id    = gene_id,
    method       = "SMR-protein",
    entity_label = gene_id,
    panel        = as.character(r$PANEL),
    statistic    = suppressWarnings(as.numeric(r$b_SMR)),
    se           = suppressWarnings(as.numeric(r$se_SMR)),
    p            = suppressWarnings(as.numeric(r$p_SMR)),
    fdr          = suppressWarnings(as.numeric(r$p_SMR.FDR)),
    evidence     = ev
  )
}

.gd_long_drug_magma <- function(gd, gwas) {
  m <- safe_access(gd_read(gd, gwas, "tx/drug"), "magma")
  if (is.null(m) || nrow(m) == 0) return(NULL)
  m <- as.data.frame(m)
  .gd_long_row(
    gwas         = gwas,
    entity_type  = "drug",
    entity_id    = m$Name,
    method       = "MAGMA-drug",
    entity_label = m$Name,
    panel        = NA_character_,
    statistic    = suppressWarnings(as.numeric(m$BETA)),
    se           = suppressWarnings(as.numeric(m$SE)),
    p            = suppressWarnings(as.numeric(m$P)),
    fdr          = suppressWarnings(as.numeric(m$P.FDR)),
    n_units      = suppressWarnings(as.integer(m[["N Genes"]]))
  )
}

.gd_long_drug_gsea <- function(gd, gwas) {
  g <- safe_access(gd_read(gd, gwas, "tx/drug"), "twas_gsea")
  if (is.null(g) || nrow(g) == 0) return(NULL)
  g <- as.data.frame(g)
  # Direction is authoritative for sign; do not derive from Estimate.
  .gd_long_row(
    gwas         = gwas,
    entity_type  = "drug",
    entity_id    = g$Name,
    method       = "TWAS-GSEA-drug",
    entity_label = g$Name,
    panel        = as.character(g$Panel),
    statistic    = suppressWarnings(as.numeric(g$Estimate)),
    se           = suppressWarnings(as.numeric(g$SE)),
    p            = suppressWarnings(as.numeric(g$P)),
    fdr          = suppressWarnings(as.numeric(g$P.FDR)),
    direction    = as.character(g$Direction),
    reversal_z   = suppressWarnings(as.numeric(g$Reversal_Z)),
    n_units      = suppressWarnings(as.integer(g[["N Genes"]]))
  )
}

.gd_long_cmap_pert <- function(gd, gwas) {
  d <- safe_access(gd_read(gd, gwas, "tx/cmap"), "drug")
  if (is.null(d) || nrow(d) == 0) return(NULL)
  d <- as.data.frame(d)
  # Reduction dims collapsed into one panel key so pick_best_per_cell can
  # take the min-P representative across cell / dose / time / panel per
  # (gwas, cmap_name).
  panel_key <- paste(d$Panel, d$cell_iname, d$pert_itime, d$pert_idose,
                      sep = " / ")
  .gd_long_row(
    gwas         = gwas,
    entity_type  = "cmap",
    entity_id    = d$cmap_name,
    method       = "CMAP-perturbation",
    entity_label = d$cmap_name,
    panel        = panel_key,
    statistic    = suppressWarnings(as.numeric(d$Estimate)),
    se           = suppressWarnings(as.numeric(d$SE)),
    p            = suppressWarnings(as.numeric(d$P)),
    fdr          = suppressWarnings(as.numeric(d$P.FDR)),
    direction    = as.character(d$Direction),
    reversal_z   = suppressWarnings(as.numeric(d$Reversal_Z))
  )
}

.gd_long_cmap_moa <- function(gd, gwas) {
  d <- safe_access(gd_read(gd, gwas, "tx/cmap"), "moa")
  if (is.null(d) || nrow(d) == 0) return(NULL)
  d <- as.data.frame(d)
  panel_key <- paste(d$Panel, d$Cell_Line, sep = " / ")
  .gd_long_row(
    gwas         = gwas,
    entity_type  = "cmap",
    entity_id    = d$MOA,
    method       = "CMAP-MOA",
    entity_label = d$MOA,
    panel        = panel_key,
    statistic    = suppressWarnings(as.numeric(d$Estimate)),
    p            = suppressWarnings(as.numeric(d$P)),
    fdr          = suppressWarnings(as.numeric(d$P.FDR)),
    direction    = as.character(d$Direction),
    reversal_z   = suppressWarnings(as.numeric(d$Reversal_Z)),
    n_units      = suppressWarnings(as.integer(d[["N Drugs"]]))
  )
}

.gd_long_atc_gsea <- function(gd, gwas) {
  g <- safe_access(gd_read(gd, gwas, "tx/atc"), "twas_gsea")
  if (is.null(g) || nrow(g) == 0) return(NULL)
  g <- as.data.frame(g)
  # Direction is authoritative for sign; do not derive from Estimate.
  .gd_long_row(
    gwas         = gwas,
    entity_type  = "atc",
    entity_id    = g[["ATC Code"]],
    method       = "TWAS-GSEA-ATC",
    entity_label = paste0(g[["ATC Code"]], ": ", g[["ATC Description"]]),
    panel        = as.character(g$Panel),
    statistic    = suppressWarnings(as.numeric(g$Estimate)),
    p            = suppressWarnings(as.numeric(g$P)),
    fdr          = suppressWarnings(as.numeric(g$P.FDR)),
    direction    = as.character(g$Direction),
    reversal_z   = suppressWarnings(as.numeric(g$Reversal_Z)),
    n_units      = suppressWarnings(as.integer(g[["N Drugs"]]))
  )
}

#' Build the cross-GWAS long tibble
#'
#' Assembles one row per (gwas, entity_type, entity_id, method, panel) from the
#' per-GWAS blocks in `gd`. Views filter and pivot this once-built tibble; they
#' must not reach into the raw package. Missing blocks are skipped silently
#' (the block simply contributes zero rows).
#'
#' @param gd A gd_result opened with gd_open()
#' @param gwas_vec Character vector of GWAS names to include
#' @param entity_types Which entity types to populate. Phase 1 supports
#'   "tissue" and "atc".
#' @return data.table with the canonical long-format columns (see .gd_long_cols)
build_comparison_long <- function(gd, gwas_vec,
                                  entity_types = c("tissue", "atc", "gene",
                                                    "locus", "drug", "cmap")) {
  if (length(gwas_vec) == 0) return(.gd_empty_long())
  parts <- list()
  for (g in gwas_vec) {
    if ("tissue" %in% entity_types) parts[[length(parts) + 1L]] <- .gd_long_tissue(gd, g)
    if ("atc"    %in% entity_types) {
      parts[[length(parts) + 1L]] <- .gd_long_atc_magma(gd, g)
      parts[[length(parts) + 1L]] <- .gd_long_atc_gsea(gd, g)
    }
    if ("gene"   %in% entity_types) {
      parts[[length(parts) + 1L]] <- .gd_long_gene_magma(gd, g)
      parts[[length(parts) + 1L]] <- .gd_long_gene_twas_fusion(gd, g)
      parts[[length(parts) + 1L]] <- .gd_long_gene_smr_exp(gd, g)
      parts[[length(parts) + 1L]] <- .gd_long_gene_pwas_fusion(gd, g)
      parts[[length(parts) + 1L]] <- .gd_long_gene_smr_prot(gd, g)
    }
    if ("locus"  %in% entity_types) {
      parts[[length(parts) + 1L]] <- .gd_long_locus_clump(gd, g)
      parts[[length(parts) + 1L]] <- .gd_long_locus_cojo(gd, g)
    }
    if ("drug"   %in% entity_types) {
      parts[[length(parts) + 1L]] <- .gd_long_drug_magma(gd, g)
      parts[[length(parts) + 1L]] <- .gd_long_drug_gsea(gd, g)
    }
    if ("cmap"   %in% entity_types) {
      parts[[length(parts) + 1L]] <- .gd_long_cmap_pert(gd, g)
      parts[[length(parts) + 1L]] <- .gd_long_cmap_moa(gd, g)
    }
  }
  parts <- Filter(function(x) !is.null(x) && nrow(x) > 0, parts)
  if (length(parts) == 0) return(.gd_empty_long())
  data.table::rbindlist(parts, use.names = TRUE, fill = TRUE)
}

#' Per-GWAS QC + power summary for the Overview tab
#'
#' Reads the lightweight metrics needed for the QC row. Missing values (e.g.
#' LDSC not run, focus block absent) become NA so the UI can render em-dashes.
#'
#' @return data.table with one row per GWAS: gwas, N, lambda_gc, h2, h2_se,
#'   n_sig_snp, n_snp_final
build_overview_qc <- function(gd, gwas_vec) {
  if (length(gwas_vec) == 0) {
    return(data.table::data.table(
      gwas = character(0), N = integer(0), lambda_gc = numeric(0),
      h2 = numeric(0), h2_se = numeric(0),
      n_sig_snp = integer(0), n_snp_final = integer(0)
    ))
  }
  rows <- lapply(gwas_vec, function(g) {
    qc <- gd_read(gd, g, "gwas_qc")
    focus_val <- safe_access(qc, "focus_dat", "val")
    cleaner_val <- safe_access(qc, "cleaner_dat", "val")
    ldsc_val <- safe_access(qc, "ldsc_dat", "val")

    # N: take a per-GWAS sample size from the clumping block if present,
    # otherwise fall back to NA. Clump rows all carry the same N.
    n_val <- NA_real_
    clump <- safe_access(gd_read(gd, g, "snp_assoc"), "clump")
    if (!is.null(clump) && nrow(clump) > 0 && "N" %in% names(clump)) {
      n_val <- suppressWarnings(max(as.numeric(clump$N), na.rm = TRUE))
      if (!is.finite(n_val)) n_val <- NA_real_
    }

    data.table::data.table(
      gwas        = g,
      N           = suppressWarnings(as.integer(n_val)),
      lambda_gc   = gd_qc_stat(qc, "lambda_gc"),
      h2          = if (!is.null(ldsc_val$obs_h2_est)) as.numeric(ldsc_val$obs_h2_est) else NA_real_,
      h2_se       = if (!is.null(ldsc_val$obs_h2_se)) as.numeric(ldsc_val$obs_h2_se) else NA_real_,
      n_sig_snp   = suppressWarnings(as.integer(gd_qc_stat(qc, "n_sig_snp"))),
      n_snp_final = if (!is.null(cleaner_val$n_snp_final))
                      suppressWarnings(as.integer(cleaner_val$n_snp_final))
                    else NA_integer_
    )
  })
  data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
}

#' Per-GWAS yield row for the Overview tab
#'
#' Counts significant entities per GWAS given the sidebar's sig basis +
#' threshold. Loci and gene counts read from raw blocks (Phase 1 does not
#' populate those in the long tibble). Tissue and ATC counts read from the
#' long tibble (best-per-cell for panelled methods).
#'
#' @param long data.table from build_comparison_long()
#' @param gd,gwas_vec As for build_comparison_long()
#' @param sig_basis "fdr" or "p"
#' @param sig_threshold numeric (default 0.05)
build_overview_yield <- function(long, gd, gwas_vec,
                                  sig_basis = "fdr", sig_threshold = 0.05) {
  # sig_basis is retained in the signature for API stability but is
  # ignored — the Overview always counts FDR-significant entities.
  is_sig <- function(v) !is.na(v) & v < sig_threshold

  # High-confidence gene count for a per-panel-per-feature results table.
  # `fdr_col` is the FDR column name; `coloc_expr` is a quoted expression
  # evaluated inside the data.table with `.SD` — either COLOC_logical
  # (FUSION) or p_HEIDI > 0.05 (SMR). Returns count of UNIQUE gene
  # symbols meeting both thresholds across all panels.
  hc_gene_count <- function(dt, fdr_col, coloc_expr, id_col) {
    if (is.null(dt) || nrow(dt) == 0 || !(fdr_col %in% names(dt))) return(NA_integer_)
    fdr_ok   <- is_sig(suppressWarnings(as.numeric(dt[[fdr_col]])))
    coloc_ok <- eval(coloc_expr, envir = dt)
    ids <- dt[[id_col]][fdr_ok & !is.na(coloc_ok) & coloc_ok]
    length(unique(ids[!is.na(ids) & nzchar(as.character(ids))]))
  }

  rows <- lapply(gwas_vec, function(g) {
    # Loci: use clumping as the canonical locus count.
    clump <- safe_access(gd_read(gd, g, "snp_assoc"), "clump")
    n_loci <- if (!is.null(clump)) nrow(clump) else NA_integer_

    # TWAS-FUSION HC genes: TWAS.P.FDR < threshold AND COLOC_logical.
    twas <- safe_access(gd_read(gd, g, "mol_assoc/exp/fusion"), "res")
    n_twas_hc <- hc_gene_count(twas, "TWAS.P.FDR", quote(COLOC_logical), "Gene Symbol")

    # PWAS-FUSION HC genes: pwas_all.P.FDR < threshold AND COLOC_logical.
    # Block layout differs from TWAS-FUSION: SMR and protein-FUSION expose
    # their results table under "results" (not "res"). See gd_read layout.
    pwas <- safe_access(gd_read(gd, g, "mol_assoc/protein/fusion"), "results")
    n_pwas_hc <- hc_gene_count(pwas, "pwas_all.P.FDR", quote(COLOC_logical), "Gene Symbol")

    # SMR-eQTL HC genes: p_SMR.FDR < threshold AND p_HEIDI > 0.05.
    smr_e <- safe_access(gd_read(gd, g, "mol_assoc/exp/smr"), "results")
    n_smr_expr_hc <- hc_gene_count(smr_e, "p_SMR.FDR",
                                     quote(suppressWarnings(as.numeric(p_HEIDI)) > 0.05),
                                     "Gene Symbol")

    # SMR-pQTL HC genes: p_SMR.FDR < threshold AND p_HEIDI > 0.05.
    smr_p <- safe_access(gd_read(gd, g, "mol_assoc/protein/smr"), "results")
    n_smr_prot_hc <- hc_gene_count(smr_p, "p_SMR.FDR",
                                     quote(suppressWarnings(as.numeric(p_HEIDI)) > 0.05),
                                     "Gene Symbol")

    # SuSiE fine-mapping: count of unique genes containing a full 95%
    # credible set (L1 output).
    finemap <- gd_read(gd, g, "mol_assoc/finemap")
    finemap_ids <- safe_access(finemap, "L1")
    n_susie <- if (is.null(finemap_ids)) NA_integer_
                else length(unique(finemap_ids[!is.na(finemap_ids) & nzchar(as.character(finemap_ids))]))

    # Sig tissues: from long (MAGMA-tissue).
    t_slice <- long[gwas == g & method == "MAGMA-tissue"]
    n_tissues <- if (nrow(t_slice) > 0) sum(is_sig(t_slice$fdr)) else NA_integer_

    # Sig ATC classes: best-per-cell then count sig.
    a_slice <- long[gwas == g & entity_type == "atc"]
    n_atc <- if (nrow(a_slice) > 0) {
      best <- pick_best_per_cell(a_slice, c("gwas", "method", "entity_id"))
      sum(is_sig(best$fdr))
    } else NA_integer_

    # Sig drugs: best-per-cell across MAGMA-drug + TWAS-GSEA-drug, then
    # count unique drug entities significant in ANY method.
    d_slice <- long[gwas == g & entity_type == "drug"]
    n_drugs <- if (nrow(d_slice) > 0) {
      best <- pick_best_per_cell(d_slice, c("gwas", "method", "entity_id"))
      length(unique(best$entity_id[is_sig(best$fdr)]))
    } else NA_integer_

    data.table::data.table(
      gwas          = g,
      n_loci        = suppressWarnings(as.integer(n_loci)),
      n_twas_hc     = suppressWarnings(as.integer(n_twas_hc)),
      n_pwas_hc     = suppressWarnings(as.integer(n_pwas_hc)),
      n_smr_expr_hc = suppressWarnings(as.integer(n_smr_expr_hc)),
      n_smr_prot_hc = suppressWarnings(as.integer(n_smr_prot_hc)),
      n_susie_hc    = suppressWarnings(as.integer(n_susie)),
      n_tissues_sig = suppressWarnings(as.integer(n_tissues)),
      n_atc_sig     = suppressWarnings(as.integer(n_atc)),
      n_drugs_sig   = suppressWarnings(as.integer(n_drugs))
    )
  })
  data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
}

#' Reduce panel dimension by taking the smallest p per (group) cell
#'
#' For TWAS/PWAS/SMR-style entities that appear under multiple panels, pick the
#' row with the smallest p-value within each group. Ties on p keep the first
#' row (stable). Rows with NA p are dropped from the pick but reappear if the
#' whole group is NA.
pick_best_per_cell <- function(long, group_cols = c("gwas", "entity_id")) {
  if (nrow(long) == 0) return(long)
  dt <- data.table::copy(long)
  data.table::setorderv(dt, c(group_cols, "p"), na.last = TRUE)
  dt[, .SD[1L], by = group_cols]
}

#' Recurrence count: number of GWAS in which an entity is significant
#'
#' @param long already-filtered long tibble (typically after
#'   apply_comparison_filters()); one row per (gwas, entity_id) is assumed.
#' @return data.table with entity_id + recurrence + min_p, ordered by
#'   recurrence desc then min_p asc.
count_recurrence <- function(long) {
  if (nrow(long) == 0) {
    return(data.table::data.table(entity_id = character(0),
                                    recurrence = integer(0),
                                    min_p = numeric(0)))
  }
  long[, .(recurrence = .N, min_p = min(p, na.rm = TRUE)),
       by = entity_id][order(-recurrence, min_p)]
}

#' Pivot a long slice into a wide entity x GWAS matrix
#'
#' Enforces that the returned matrix has exactly `length(gwas_vec)` columns,
#' one per selected GWAS in the requested order, filled with NA where a GWAS
#' contributed no row. This guards against the "empty column silently
#' disappears" bug that has bitten the extraction scripts elsewhere.
#'
#' @param long data.table (already filtered/reduced to one row per
#'   (gwas, entity_id) — call pick_best_per_cell() first when needed)
#' @param values_col column name to place in cells (e.g. "fdr", "p", "statistic")
#' @param gwas_vec column order for the resulting matrix
#' @param entity_order optional character vector to force row order
#' @return numeric matrix with rownames = entity_id, colnames = gwas_vec
pivot_matrix <- function(long, values_col, gwas_vec, entity_order = NULL) {
  if (length(gwas_vec) == 0) {
    m <- matrix(numeric(0), nrow = 0, ncol = 0)
    return(m)
  }
  if (nrow(long) == 0) {
    m <- matrix(NA_real_, nrow = 0, ncol = length(gwas_vec))
    colnames(m) <- gwas_vec
    return(m)
  }
  ent <- if (!is.null(entity_order)) unique(entity_order) else unique(long$entity_id)
  m <- matrix(NA_real_, nrow = length(ent), ncol = length(gwas_vec),
              dimnames = list(ent, gwas_vec))
  ri <- match(long$entity_id, ent)
  ci <- match(long$gwas, gwas_vec)
  ok <- !is.na(ri) & !is.na(ci)
  m[cbind(ri[ok], ci[ok])] <- suppressWarnings(as.numeric(long[[values_col]][ok]))
  m
}

#' Same as pivot_matrix but for a character-valued column (direction, etc.)
pivot_matrix_chr <- function(long, values_col, gwas_vec, entity_order = NULL) {
  if (length(gwas_vec) == 0) return(matrix(character(0), nrow = 0, ncol = 0))
  if (nrow(long) == 0) {
    m <- matrix(NA_character_, nrow = 0, ncol = length(gwas_vec))
    colnames(m) <- gwas_vec
    return(m)
  }
  ent <- if (!is.null(entity_order)) unique(entity_order) else unique(long$entity_id)
  m <- matrix(NA_character_, nrow = length(ent), ncol = length(gwas_vec),
              dimnames = list(ent, gwas_vec))
  ri <- match(long$entity_id, ent)
  ci <- match(long$gwas, gwas_vec)
  ok <- !is.na(ri) & !is.na(ci)
  m[cbind(ri[ok], ci[ok])] <- as.character(long[[values_col]][ok])
  m
}

#' Build a long tibble of genetic correlations (LDSC gencor)
#'
#' The GenoDisc pipeline computes bivariate LDSC between each in-bundle GWAS
#' and a fixed set of external reference traits (the `gencor_gwas_list` in
#' the config). This helper reads each selected bundle GWAS's
#' `gwas_qc$ldsc_gencor_dat$table` and returns the rows in a long tibble
#' keyed on (bundle_gwas, reference_trait).
#'
#' Returns an empty data.table (with the expected columns) when the bundle
#' has no gencor data.
#'
#' @param gd A gd_result opened with gd_open()
#' @param gwas_vec Character vector of selected bundle GWAS names
#' @return data.table with columns: gwas, ref_name, ref_label, rg, rg_se,
#'   rg_p, rg_p_fdr, n_snps, gcov_int, trait_category
build_gencor_long <- function(gd, gwas_vec) {
  empty <- data.table::data.table(
    gwas = character(0),
    ref_name = character(0), ref_label = character(0),
    rg = numeric(0), rg_se = numeric(0),
    rg_p = numeric(0), rg_p_fdr = numeric(0),
    n_snps = integer(0), gcov_int = numeric(0),
    trait_category = character(0)
  )
  if (length(gwas_vec) == 0) return(empty)

  parts <- lapply(gwas_vec, function(g) {
    tab <- safe_access(gd_read(gd, g, "gwas_qc"), "ldsc_gencor_dat", "table")
    if (is.null(tab) || nrow(tab) == 0) return(NULL)
    tab <- as.data.frame(tab)
    data.table::data.table(
      gwas          = g,
      ref_name      = as.character(tab$name),
      ref_label     = as.character(tab$label),
      rg            = suppressWarnings(as.numeric(tab$rg)),
      rg_se         = suppressWarnings(as.numeric(tab$rg_se)),
      rg_p          = suppressWarnings(as.numeric(tab$rg_p)),
      rg_p_fdr      = suppressWarnings(as.numeric(tab$rg_p_fdr)),
      n_snps        = suppressWarnings(as.integer(tab$n_snps)),
      gcov_int      = if ("gcov_int" %in% names(tab))
                        suppressWarnings(as.numeric(tab$gcov_int))
                      else NA_real_,
      trait_category = if ("trait_category" %in% names(tab))
                         as.character(tab$trait_category)
                       else NA_character_
    )
  })
  parts <- Filter(function(x) !is.null(x) && nrow(x) > 0, parts)
  if (length(parts) == 0) return(empty)
  data.table::rbindlist(parts, use.names = TRUE, fill = TRUE)
}

#' Per-GWAS QC statistics for the GWAS QC → QC Summary table
#'
#' Combines the fields shown in the single-GWAS QC Summary (N variants pre-QC,
#' identified genome build, N variants post-QC, Lambda GC, Max chi², N GWS
#' variants) plus the GWAS label if provided. Powers both single-GWAS mode
#' (one row) and the wide compare-mode table (one column per GWAS).
#'
#' @param gd A gd_result opened with gd_open()
#' @param gwas_vec Character vector of GWAS names
#' @param gwas_list Optional data.frame from gd_config(gd)$gwas_list (uses
#'   `name` / `label` columns to supply human-readable labels).
#' @return data.table with columns: gwas, label, n_var_orig, build,
#'   n_snp_final, lambda_gc, max_chi2, n_sig_snp
build_qc_summary_long <- function(gd, gwas_vec, gwas_list = NULL) {
  if (length(gwas_vec) == 0) {
    return(data.table::data.table(
      gwas = character(0), label = character(0),
      n_var_orig = integer(0), build = character(0),
      n_snp_final = integer(0),
      lambda_gc = numeric(0), max_chi2 = numeric(0),
      n_sig_snp = integer(0),
      int_est = numeric(0), int_se = numeric(0)
    ))
  }
  rows <- lapply(gwas_vec, function(g) {
    q <- gd_read(gd, g, "gwas_qc")
    cv <- safe_access(q, "cleaner_dat", "val")
    lv <- safe_access(q, "ldsc_dat", "val")
    build_val <- if (!is.null(cv$build$build)) as.character(cv$build$build)
                 else NA_character_
    if (is.null(build_val) || length(build_val) == 0 || is.na(build_val))
      build_val <- "Unknown"
    label <- if (!is.null(gwas_list) && "name" %in% names(gwas_list)) {
      lbl <- gwas_list$label[gwas_list$name == g]
      if (length(lbl) == 0 || is.na(lbl[1L])) g else as.character(lbl[1L])
    } else g
    data.table::data.table(
      gwas        = g,
      label       = label,
      n_var_orig  = suppressWarnings(as.integer(cv$n_var_orig)),
      build       = build_val,
      n_snp_final = suppressWarnings(as.integer(cv$n_snp_final)),
      lambda_gc   = gd_qc_stat(q, "lambda_gc"),
      max_chi2    = gd_qc_stat(q, "max_chi2"),
      n_sig_snp   = suppressWarnings(as.integer(gd_qc_stat(q, "n_sig_snp"))),
      int_est     = if (!is.null(lv$int_est)) as.numeric(lv$int_est) else NA_real_,
      int_se      = if (!is.null(lv$int_se))  as.numeric(lv$int_se)  else NA_real_
    )
  })
  data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
}

#' Per-GWAS LDSC SNP-heritability + LDSC intercept for the SNP-h² tab
#'
#' Reads `gwas_qc$ldsc_dat$val` for each GWAS and returns a one-row-per-GWAS
#' data.table. Missing values (LDSC not run, block absent) become NA and the
#' UI renders them as em-dashes.
#'
#' @param gd A gd_result opened with gd_open()
#' @param gwas_vec Character vector of GWAS names
#' @return data.table with columns: gwas, h2, h2_se, int_est, int_se
build_heritability_long <- function(gd, gwas_vec) {
  if (length(gwas_vec) == 0) {
    return(data.table::data.table(
      gwas = character(0),
      h2 = numeric(0), h2_se = numeric(0),
      int_est = numeric(0), int_se = numeric(0)
    ))
  }
  rows <- lapply(gwas_vec, function(g) {
    v <- safe_access(gd_read(gd, g, "gwas_qc"), "ldsc_dat", "val")
    data.table::data.table(
      gwas    = g,
      h2      = if (!is.null(v$obs_h2_est)) as.numeric(v$obs_h2_est) else NA_real_,
      h2_se   = if (!is.null(v$obs_h2_se)) as.numeric(v$obs_h2_se) else NA_real_,
      int_est = if (!is.null(v$int_est))   as.numeric(v$int_est)   else NA_real_,
      int_se  = if (!is.null(v$int_se))    as.numeric(v$int_se)    else NA_real_
    )
  })
  data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
}

#' Order the selected GWAS for display
#'
#' Applied to the column order of comparison matrices. Requires the QC table
#' from build_overview_qc() when ordering by N / lambda / h2 / n_sig_snp.
#'
#' @param gwas_vec character vector as returned by selected_gwas_multi()
#' @param mode one of "as_selected", "alphabetical", "n", "n_sig_snp", "h2"
#' @param qc optional data.table from build_overview_qc(); required for
#'   metric-based modes
order_gwas <- function(gwas_vec, mode = "as_selected", qc = NULL) {
  if (length(gwas_vec) <= 1L) return(gwas_vec)
  mode <- if (is.null(mode)) "as_selected" else as.character(mode)
  if (mode == "alphabetical") return(sort(gwas_vec))
  metric_col <- switch(mode,
    n           = "N",
    n_sig_snp   = "n_sig_snp",
    h2          = "h2",
    NA_character_
  )
  if (!is.na(metric_col) && !is.null(qc) && metric_col %in% names(qc)) {
    q <- qc[gwas %in% gwas_vec]
    v <- suppressWarnings(as.numeric(q[[metric_col]]))
    ord <- order(-v, q$gwas, na.last = TRUE)
    return(as.character(q$gwas[ord]))
  }
  gwas_vec
}

#' Apply the shared comparison-mode filters
#'
#' Filter order (each step preserves the "column present but empty" invariant
#' downstream since we never drop GWAS from the reference vector):
#'   1. panel restriction (if `cf$panel` is a non-empty character)
#'   2. evidence-required (if `cf$evidence_required` is TRUE, drop rows where
#'      evidence is FALSE / NA)
#'   3. significance threshold on the chosen basis (drop rows above threshold)
#'   4. "significant in >= k GWAS" (drop entities failing the recurrence test)
#'
#' @param long data.table from build_comparison_long()
#' @param cf named list with:
#'   sig_basis ("fdr" | "p"), sig_threshold (numeric),
#'   evidence_required (logical), panel (character or NULL/""), k_min (integer)
#' @return filtered data.table (still in long form)
apply_comparison_filters <- function(long, cf = list()) {
  if (nrow(long) == 0) return(long)
  dt <- data.table::copy(long)

  panel_pick <- cf$panel
  if (!is.null(panel_pick) && length(panel_pick) > 0 && nzchar(panel_pick[1L])) {
    dt <- dt[is.na(panel) | panel %in% panel_pick]
  }

  if (isTRUE(cf$evidence_required)) {
    dt <- dt[!is.na(evidence) & evidence == TRUE]
  }

  basis <- if (identical(cf$sig_basis, "p")) "p" else "fdr"
  thr <- if (!is.null(cf$sig_threshold)) as.numeric(cf$sig_threshold) else 0.05

  # Recurrence filter: keep entities significant in >= k_min GWAS on the basis.
  # k_min = 1 means "at least one GWAS", which drops entities that are not
  # significant anywhere. Only skip the filter when k_min is not provided or
  # is explicitly < 1.
  k_min <- if (!is.null(cf$k_min)) as.integer(cf$k_min) else 1L
  if (!is.na(k_min) && k_min >= 1L) {
    dt[, .gd_sig := !is.na(.SD[[1L]]) & .SD[[1L]] < thr, .SDcols = basis]
    counts <- dt[, .(k = sum(.gd_sig)), by = entity_id]
    keep_ids <- counts[k >= k_min, entity_id]
    dt <- dt[entity_id %in% keep_ids]
    dt[, .gd_sig := NULL]
  }

  dt
}
