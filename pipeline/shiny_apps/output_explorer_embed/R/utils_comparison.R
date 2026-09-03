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
                                                    "locus", "drug")) {
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
  key_col <- if (identical(sig_basis, "p")) "p" else "fdr"
  is_sig <- function(v) !is.na(v) & v < sig_threshold
  rows <- lapply(gwas_vec, function(g) {
    # Loci: use clumping as the canonical locus count.
    clump <- safe_access(gd_read(gd, g, "snp_assoc"), "clump")
    n_loci <- if (!is.null(clump)) nrow(clump) else NA_integer_

    # Sig genes: MAGMA gene-level FDR / P.
    magma <- gd_read(gd, g, "mol_assoc/magma")
    n_genes <- if (!is.null(magma) && "P.FDR" %in% names(magma)) {
      v <- if (identical(key_col, "p")) magma$P else magma$P.FDR
      sum(is_sig(suppressWarnings(as.numeric(v))))
    } else NA_integer_

    # Sig tissues: from long (MAGMA-tissue).
    t_slice <- long[gwas == g & method == "MAGMA-tissue"]
    n_tissues <- if (nrow(t_slice) > 0) sum(is_sig(t_slice[[key_col]])) else NA_integer_

    # Sig ATC classes: best-per-cell then count sig. MAGMA-ATC has one panel;
    # TWAS-GSEA-ATC gets best-per-cell across Panel.
    a_slice <- long[gwas == g & entity_type == "atc"]
    n_atc <- if (nrow(a_slice) > 0) {
      best <- pick_best_per_cell(a_slice, c("gwas", "method", "entity_id"))
      sum(is_sig(best[[key_col]]))
    } else NA_integer_

    data.table::data.table(
      gwas = g,
      n_loci = suppressWarnings(as.integer(n_loci)),
      n_genes_sig = suppressWarnings(as.integer(n_genes)),
      n_tissues_sig = suppressWarnings(as.integer(n_tissues)),
      n_atc_sig = suppressWarnings(as.integer(n_atc))
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

#' Build a long tibble of pairwise genetic correlations (LDSC gencor)
#'
#' Reads each selected GWAS's `gwas_qc$ldsc_gencor_dat` block and keeps only
#' pairs (i, j) where BOTH i and j are in `gwas_vec`. Where both directions
#' are present, the two rows are averaged (LDSC rg is symmetric; identical
#' point estimates in practice, but bundles may report slightly different
#' SEs). Returns an empty data.table (with the expected columns) when the
#' bundle carries no gencor data.
#'
#' @param gd A gd_result opened with gd_open()
#' @param gwas_vec Character vector of GWAS names to include
#' @return data.table with columns: gwas_i, gwas_j, rg, rg_se, rg_p,
#'   rg_p_fdr, n_snps
build_gencor_long <- function(gd, gwas_vec) {
  empty <- data.table::data.table(
    gwas_i = character(0), gwas_j = character(0),
    rg = numeric(0), rg_se = numeric(0),
    rg_p = numeric(0), rg_p_fdr = numeric(0),
    n_snps = integer(0)
  )
  if (length(gwas_vec) == 0) return(empty)

  parts <- lapply(gwas_vec, function(g) {
    tab <- safe_access(gd_read(gd, g, "gwas_qc"), "ldsc_gencor_dat")
    if (is.null(tab) || nrow(tab) == 0) return(NULL)
    tab <- as.data.frame(tab)
    # `name` in the gencor block is the OTHER trait's canonical name. Keep
    # only pairs whose other trait is also in the selection.
    keep <- tab$name %in% gwas_vec
    if (!any(keep)) return(NULL)
    tab <- tab[keep, , drop = FALSE]
    data.table::data.table(
      gwas_i   = g,
      gwas_j   = as.character(tab$name),
      rg       = suppressWarnings(as.numeric(tab$rg)),
      rg_se    = suppressWarnings(as.numeric(tab$rg_se)),
      rg_p     = suppressWarnings(as.numeric(tab$rg_p)),
      rg_p_fdr = suppressWarnings(as.numeric(tab$rg_p_fdr)),
      n_snps   = suppressWarnings(as.integer(tab$n_snps))
    )
  })
  parts <- Filter(function(x) !is.null(x) && nrow(x) > 0, parts)
  if (length(parts) == 0) return(empty)
  long <- data.table::rbindlist(parts, use.names = TRUE, fill = TRUE)

  # Symmetrise: for every (i, j) present, also emit (j, i) with the same
  # values; where both directions exist naturally we average them.
  long[, pair_key := paste(pmin(gwas_i, gwas_j), pmax(gwas_i, gwas_j), sep = "|")]
  avg <- long[, .(
    rg       = mean(rg, na.rm = TRUE),
    rg_se    = mean(rg_se, na.rm = TRUE),
    rg_p     = mean(rg_p, na.rm = TRUE),
    rg_p_fdr = mean(rg_p_fdr, na.rm = TRUE),
    n_snps   = suppressWarnings(as.integer(mean(n_snps, na.rm = TRUE)))
  ), by = pair_key]
  # Expand each averaged pair back into two directional rows.
  parts_ij <- strsplit(avg$pair_key, "|", fixed = TRUE)
  a <- vapply(parts_ij, `[[`, character(1), 1L)
  b <- vapply(parts_ij, `[[`, character(1), 2L)
  out <- rbind(
    data.table::data.table(gwas_i = a, gwas_j = b, rg = avg$rg, rg_se = avg$rg_se,
                            rg_p = avg$rg_p, rg_p_fdr = avg$rg_p_fdr, n_snps = avg$n_snps),
    data.table::data.table(gwas_i = b, gwas_j = a, rg = avg$rg, rg_se = avg$rg_se,
                            rg_p = avg$rg_p, rg_p_fdr = avg$rg_p_fdr, n_snps = avg$n_snps)
  )
  # Add the diagonal (rg = 1 for self-correlation) so heatmaps render a full
  # square without gaps on the identity.
  diag_rows <- data.table::data.table(
    gwas_i = gwas_vec, gwas_j = gwas_vec,
    rg = 1, rg_se = NA_real_, rg_p = NA_real_, rg_p_fdr = NA_real_,
    n_snps = NA_integer_
  )
  data.table::rbindlist(list(out, diag_rows), use.names = TRUE, fill = TRUE)
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
