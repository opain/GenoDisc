#' Safely access nested list elements
#'
#' @param data A list
#' @param ... Character keys to traverse
#' @return The value at the path, or NULL if any level is missing
safe_access <- function(data, ...) {
  keys <- list(...)
  result <- data
  for (key in keys) {
    if (is.null(result) || !is.list(result) || !(key %in% names(result))) {
      return(NULL)
    }
    result <- result[[key]]
  }
  result
}

#' Build molecular association summary data
#'
#' @param gd A gd_result opened with gd_open()
#' @param gwas GWAS name
#' @param cf Named list of config flags from parse_config_flags()
#' @return data.frame with columns: Panel, ID, Z, Sig, Coloc, Method, Type
build_mol_assoc_data <- function(gd, gwas, cf) {
  all_func_res <- NULL

  if (cf$finemap) {
    finemap <- gd_read(gd, gwas, "mol_assoc/finemap")
    finemap_ids <- safe_access(finemap, "L1")
    if (!is.null(finemap_ids) && length(finemap_ids) > 0) {
      all_func_res <- rbind(all_func_res, data.frame(
        Panel = "SuSie (L=1)", ID = finemap_ids, Z = 1, Sig = F, Coloc = F,
        Method = "SNP\nFine-mapping", Type = ''))
    }
  }

  if (cf$twas) {
    fusion_res <- safe_access(gd_read(gd, gwas, "mol_assoc/exp/fusion"), "res")
    if (!is.null(fusion_res)) {
      twas_tmp <- data.table(
        Panel = fusion_res$PANEL,
        ID = fusion_res$`Gene Symbol`,
        Z = fusion_res$TWAS.Z,
        Sig = fusion_res$TWAS.P.FDR < 0.05,
        Coloc = fusion_res$COLOC_logical)
      twas_tmp$Method <- 'FUSION'
      twas_tmp$Type <- 'Expr.'
      twas_tmp$Type[grepl('SPLIC', twas_tmp$Panel, ignore.case = T)] <- 'Splice'
      twas_tmp <- twas_tmp[order(-abs(twas_tmp$Z)), ]
      twas_tmp <- twas_tmp[!duplicated(paste0(twas_tmp$Panel, twas_tmp$ID)), ]
      twas_tmp$Type <- factor(twas_tmp$Type, levels = c('Expr.', 'Splice'))
      twas_tmp <- twas_tmp[order(twas_tmp$Type), ]
      all_func_res <- rbind(all_func_res, twas_tmp)
    }
  }

  if (cf$smr_expression) {
    smr_res <- safe_access(gd_read(gd, gwas, "mol_assoc/exp/smr"), "results")
    if (!is.null(smr_res)) {
      smr_expr_id <- smr_res$`Gene Symbol`
      smr_expr_id[is.na(smr_expr_id)] <- smr_res$`Ensembl ID`[is.na(smr_expr_id)]
      smr_tmp <- data.table(
        Panel = smr_res$PANEL, ID = smr_expr_id,
        Z = smr_res$b_SMR / smr_res$se_SMR,
        Sig = smr_res$p_SMR.FDR < 0.05,
        Coloc = smr_res$p_HEIDI > 0.05)
      smr_tmp$Method <- 'SMR'
      smr_tmp$Type <- 'Expr.'
      all_func_res <- rbind(all_func_res, smr_tmp)
    }
  }

  if (any(cf$pwas_panel_rosmap, cf$pwas_panel_banner)) {
    pwas_res <- safe_access(gd_read(gd, gwas, "mol_assoc/protein/fusion"), "results")
    if (!is.null(pwas_res)) {
      pwas_tmp <- data.table(
        Panel = pwas_res$PANEL, ID = pwas_res$`Gene Symbol`,
        Z = pwas_res$pwas_all.Z,
        Sig = pwas_res$pwas_all.P.FDR < 0.05,
        Coloc = pwas_res$COLOC_logical)
      pwas_tmp <- pwas_tmp[order(-abs(pwas_tmp$Z)), ]
      pwas_tmp <- pwas_tmp[!duplicated(paste0(pwas_tmp$Panel, pwas_tmp$ID)), ]
      pwas_tmp$Method <- 'FUSION'
      pwas_tmp$Type <- 'Protein'
      all_func_res <- rbind(all_func_res, pwas_tmp)
    }
  }

  if (cf$smr_protein_panel_rosmap) {
    smr_prot <- safe_access(gd_read(gd, gwas, "mol_assoc/protein/smr"), "results")
    if (!is.null(smr_prot)) {
      smr_prot_tmp <- data.table(
        Panel = smr_prot$PANEL, ID = smr_prot$`Gene Symbol`,
        Z = smr_prot$b_SMR / smr_prot$se_SMR,
        Sig = smr_prot$p_SMR.FDR < 0.05,
        Coloc = smr_prot$p_HEIDI > 0.05)
      smr_prot_tmp <- smr_prot_tmp[order(-abs(smr_prot_tmp$Z)), ]
      smr_prot_tmp <- smr_prot_tmp[!duplicated(paste0(smr_prot_tmp$Panel, smr_prot_tmp$ID)), ]
      smr_prot_tmp$Method <- 'SMR'
      smr_prot_tmp$Type <- 'Protein'
      all_func_res <- rbind(all_func_res, smr_prot_tmp)
    }
  }

  if (cf$magma_gene) {
    magma <- gd_read(gd, gwas, "mol_assoc/magma")
    if (!is.null(magma)) {
      all_func_res <- rbind(all_func_res, data.frame(
        Panel = 'MAGMA', ID = magma$ID,
        Z = abs(qnorm(as.numeric(magma$P))),
        Sig = as.numeric(magma$P.FDR) < 0.05,
        Coloc = F, Method = 'MAGMA', Type = ''))
    }
  }

  if (cf$clump) {
    nearest <- safe_access(gd_read(gd, gwas, "mol_assoc/nearest"), "clump")
    if (!is.null(nearest) && length(nearest) > 0) {
      all_func_res <- rbind(all_func_res, data.frame(
        Panel = 'NearestGene', ID = nearest, Z = 1, Sig = F, Coloc = F,
        Method = 'Nearest\nGene', Type = ''))
    }
  }

  all_func_res
}

#' Build a gene-symbol -> genomic position lookup
#'
#' build_mol_assoc_data() keeps only the gene symbol (ID), not position. The
#' underlying method blocks that the results tables read (via $res) do carry a
#' chromosome and position, so we harvest those here to give each symbol one
#' canonical position for locus grouping. SuSiE fine-mapping and Nearest-gene
#' slots hold symbols only and contribute no position (those genes stay
#' unplaced downstream).
#'
#' @param gd A gd_result opened with gd_open()
#' @param gwas GWAS name
#' @param cf Named list of config flags from parse_config_flags()
#' @return data.frame with columns ID, CHR, BP (one row per symbol), or an
#'   empty data.frame if no position-bearing method is present.
build_gene_position_map <- function(gd, gwas, cf) {
  pos <- NULL

  add_pos <- function(id, chr, bp) {
    ok <- !is.na(id) & !is.na(chr) & !is.na(bp)
    if (!any(ok)) return(invisible(NULL))
    pos <<- rbind(pos, data.frame(
      ID = as.character(id)[ok],
      CHR = suppressWarnings(as.numeric(chr))[ok],
      BP = suppressWarnings(as.numeric(bp))[ok],
      stringsAsFactors = FALSE))
  }

  # The FUSION/SMR result table lives under $res in the table renderers but
  # $results in build_mol_assoc_data(); read whichever this package uses.
  get_res <- function(blk) {
    o <- gd_read(gd, gwas, blk)
    r <- safe_access(o, "res")
    if (is.null(r)) r <- safe_access(o, "results")
    r
  }

  # MAGMA gene table: ID = symbol, gene bounds START/STOP -> midpoint.
  if (isTRUE(cf$magma_gene)) {
    m <- gd_read(gd, gwas, "mol_assoc/magma")
    if (!is.null(m) && all(c("ID", "CHR", "START", "STOP") %in% names(m))) {
      add_pos(m$ID, m$CHR, (as.numeric(m$START) + as.numeric(m$STOP)) / 2)
    }
  }

  # FUSION expr/protein: Gene Symbol, gene bounds P0/P1 -> midpoint.
  fusion_blocks <- c()
  if (isTRUE(cf$twas)) fusion_blocks <- c(fusion_blocks, "mol_assoc/exp/fusion")
  if (any(cf$pwas_panel_rosmap, cf$pwas_panel_banner)) fusion_blocks <- c(fusion_blocks, "mol_assoc/protein/fusion")
  for (blk in fusion_blocks) {
    r <- get_res(blk)
    if (!is.null(r) && all(c("Gene Symbol", "CHR", "P0", "P1") %in% names(r))) {
      add_pos(r$`Gene Symbol`, r$CHR, (as.numeric(r$P0) + as.numeric(r$P1)) / 2)
    }
  }

  # SMR expr/protein: Gene Symbol, single position BP.
  smr_blocks <- c()
  if (isTRUE(cf$smr_expression)) smr_blocks <- c(smr_blocks, "mol_assoc/exp/smr")
  if (isTRUE(cf$smr_protein_panel_rosmap)) smr_blocks <- c(smr_blocks, "mol_assoc/protein/smr")
  for (blk in smr_blocks) {
    r <- get_res(blk)
    if (!is.null(r) && all(c("Gene Symbol", "CHR", "BP") %in% names(r))) {
      add_pos(r$`Gene Symbol`, r$CHR, r$BP)
    }
  }

  if (is.null(pos) || nrow(pos) == 0) {
    return(data.frame(ID = character(0), CHR = numeric(0), BP = numeric(0)))
  }

  # Collapse to one canonical position per symbol (CHR = first, BP = median).
  pos <- pos[!is.na(pos$CHR) & !is.na(pos$BP), ]
  agg <- aggregate(BP ~ ID, data = pos, FUN = median)
  chr <- aggregate(CHR ~ ID, data = pos, FUN = function(x) x[1])
  merge(chr, agg, by = "ID")
}

#' Load the bundled GRCh37 gene-position reference (cached)
#'
#' Reads data/gene_positions.rds (built by data/make_gene_positions.R) once and
#' caches it. Returns a list with three symbol-keyed data.frames (ID, CHR, BP):
#' by_symbol, by_synonym, by_ensembl. Returns NULL if the file is absent, in
#' which case the app falls back to per-method harvested positions.
load_gene_positions <- local({
  cache <- NULL
  loaded <- FALSE
  function() {
    if (loaded) return(cache)
    loaded <<- TRUE
    here <- tryCatch(dirname(sys.frame(1L)$ofile), error = function(e) getwd())
    candidates <- c(
      file.path("data", "gene_positions.rds"),
      file.path(here, "..", "data", "gene_positions.rds"),
      file.path(here, "data", "gene_positions.rds")
    )
    hit <- candidates[file.exists(candidates)][1]
    cache <<- if (is.na(hit)) NULL else readRDS(hit)
    cache
  }
})

#' Resolve a genomic position for each feature via one canonical reference
#'
#' Deterministically maps each gene id to a single GRCh37 position using the
#' bundled reference (symbol -> synonym -> Ensembl), falling back to positions
#' harvested from the method blocks for anything the reference does not cover.
#' Only ids that resolve are returned (unresolved features stay "Unplaced").
#'
#' @param ids Character vector of feature ids (gene symbols / Ensembl ids)
#' @param gd,gwas,cf As for build_gene_position_map() (fallback source)
#' @return data.frame(ID, CHR, BP) for the resolved subset of `ids`
resolve_feature_positions <- function(ids, gd, gwas, cf) {
  ids <- unique(as.character(ids))
  ids <- ids[!is.na(ids) & ids != "" & ids != "Placeholder"]
  out <- data.frame(ID = ids, CHR = NA_real_, BP = NA_real_, stringsAsFactors = FALSE)
  if (nrow(out) == 0) return(out)

  fill_from <- function(map, eligible) {
    idx <- match(out$ID, map$ID)
    take <- eligible & !is.na(idx)
    out$CHR[take] <<- map$CHR[idx[take]]
    out$BP[take]  <<- map$BP[idx[take]]
  }

  ref <- load_gene_positions()
  if (!is.null(ref)) {
    fill_from(ref$by_symbol,  is.na(out$BP))
    fill_from(ref$by_synonym, is.na(out$BP))
    fill_from(ref$by_ensembl, is.na(out$BP) & grepl("^ENSG", out$ID))
  }

  # Fallback for ids the reference does not cover.
  if (any(is.na(out$BP))) {
    harvest <- build_gene_position_map(gd, gwas, cf)
    if (!is.null(harvest) && nrow(harvest) > 0) fill_from(harvest, is.na(out$BP))
  }

  out[!is.na(out$BP), ]
}

#' Build drug enrichment summary data
#'
#' @param gd A gd_result
#' @param gwas GWAS name
#' @return data.frame with columns: Name, Z, P, P.FDR, ATC Code, Method, Panel
build_drug_summary_data <- function(gd, gwas) {
  drug <- gd_read(gd, gwas, "tx/drug")

  magma_gs <- safe_access(drug, "magma")
  if (!is.null(magma_gs)) {
    magma_gs$Z <- -qnorm(magma_gs$P)
    # P=0 (numerical underflow upstream) makes qnorm return Inf; Inf breaks
    # the ggplot fill scale ("'to' must be a finite number"). Drop to NA so
    # the point still renders (na.value) without dragging the colour limits
    # to infinity.
    magma_gs$Z[!is.finite(magma_gs$Z)] <- NA_real_
    magma_gs <- magma_gs[, c('Name', 'Z', 'P', 'P.FDR', 'ATC Code')]
    magma_gs$Method <- 'MAGMA'
    magma_gs$Panel <- 'MAGMA'
  }

  gcsc_gs <- safe_access(drug, "gcsc")
  if (!is.null(gcsc_gs)) {
    gcsc_gs <- gcsc_gs[, c('Name', 'Z', 'P', 'P.FDR', 'ATC Code')]
    gcsc_gs$Method <- 'GCSC'
    gcsc_gs$Panel <- 'Brain and Blood'
  }

  build_gsea <- function(slot, method_label) {
    g <- safe_access(drug, slot)
    if (is.null(g)) return(NULL)
    g$Method <- method_label
    # Reversal_Z is positive when the drug opposes the disease TWAS signature
    # (candidate therapeutic direction). For directional and non-directional
    # variants alike this is set by the format script / read function, so the
    # Shiny app no longer applies any sign flips here.
    g$Z <- g$Reversal_Z
    # See MAGMA branch above: P=0 upstream produces Reversal_Z=Inf which
    # would crash the ggplot fill scale.
    g$Z[!is.finite(g$Z)] <- NA_real_
    g <- g[, c('Name', 'Z', 'P', 'P.FDR', 'Method', 'Panel', 'ATC Code')]

    g_all <- g
    for (i in unique(g_all$Panel)) {
      g_i <- g[g$Panel == i, ]
      g_other <- g[g$Panel != i, ]
      missing_names <- unique(g_other$Name[!(g_other$Name %in% g_i$Name)])
      if (length(missing_names) > 0) {
        g_rest <- data.frame(
          Name = missing_names,
          Z = NA, P = NA, P.FDR = NA, Method = method_label, Panel = i, ATC_Code = NA)
        names(g_rest) <- gsub('ATC_Code', 'ATC Code', names(g_rest))
        g_all <- rbind(g_all, g_rest)
      }
    }
    g_all
  }

  gsea_gs <- build_gsea("twas_gsea", "TWAS-GSEA")
  gsea_gs_nondir <- build_gsea("twas_gsea_nondir", "TWAS-GSEA (non-dir)")

  do.call(rbind, Filter(Negate(is.null), list(magma_gs, gcsc_gs, gsea_gs, gsea_gs_nondir)))
}

#' Build ATC enrichment summary data
#'
#' @param gd A gd_result
#' @param gwas GWAS name
#' @return data.frame with columns: Name, Z, FDR_Sig, Nom_Sig, Method, Panel
build_atc_summary_data <- function(gd, gwas) {
  atc <- gd_read(gd, gwas, "tx/atc")

  magma_gs_atc <- safe_access(atc, "magma")
  if (!is.null(magma_gs_atc)) {
    magma_gs_atc$Z <- -qnorm(magma_gs_atc$P)
    magma_gs_atc$Z[!is.finite(magma_gs_atc$Z)] <- NA_real_
    magma_gs_atc$FDR_Sig <- magma_gs_atc$P.FDR < 0.05
    magma_gs_atc$Nom_Sig <- magma_gs_atc$P < 0.05
    magma_gs_atc$Name <- paste0(magma_gs_atc$`ATC Code`, ': ', magma_gs_atc$`ATC Description`)
    magma_gs_atc$Method <- 'MAGMA'
    magma_gs_atc$Panel <- 'MAGMA'
    magma_gs_atc <- magma_gs_atc[, c("Name", "Z", "FDR_Sig", "Nom_Sig", "Method", "Panel"), with = F]
  }

  gcsc_gs_atc <- safe_access(atc, "gcsc")
  if (!is.null(gcsc_gs_atc)) {
    gcsc_gs_atc$Z <- -qnorm(gcsc_gs_atc$P)
    gcsc_gs_atc$Z[!is.finite(gcsc_gs_atc$Z)] <- NA_real_
    gcsc_gs_atc$FDR_Sig <- gcsc_gs_atc$P.FDR < 0.05
    gcsc_gs_atc$Nom_Sig <- gcsc_gs_atc$P < 0.05
    gcsc_gs_atc$Name <- paste0(gcsc_gs_atc$`ATC Code`, ': ', gcsc_gs_atc$`ATC Description`)
    gcsc_gs_atc$Method <- 'GCSC'
    gcsc_gs_atc$Panel <- 'GCSC'
    gcsc_gs_atc <- gcsc_gs_atc[, c("Name", "Z", "FDR_Sig", "Nom_Sig", "Method", "Panel"), with = F]
  }

  build_gsea_atc <- function(slot, method_label) {
    g <- safe_access(atc, slot)
    if (is.null(g)) return(NULL)
    # P.FDR is already computed in the RDS (per-panel by the read function);
    # Reversal_Z is positive when the class opposes the disease TWAS signature.
    g$Z <- g$Reversal_Z
    g$Z[!is.finite(g$Z)] <- NA_real_
    g$FDR_Sig <- g$P.FDR < 0.05
    g$Nom_Sig <- g$P < 0.05
    g$Name <- paste0(g$`ATC Code`, ': ', g$`ATC Description`)
    g$Method <- method_label
    g <- g[, c("Name", "Z", "FDR_Sig", "Nom_Sig", "Method", "Panel"), with = F]

    g_all <- g
    for (i in unique(g_all$Panel)) {
      g_i <- g[g$Panel == i, ]
      g_other <- g[g$Panel != i, ]
      missing_atc_names <- unique(g_other$Name[!(g_other$Name %in% g_i$Name)])
      if (length(missing_atc_names) > 0) {
        g_rest <- data.frame(
          Name = missing_atc_names,
          Z = NA, FDR_Sig = NA, Nom_Sig = NA, Method = method_label, Panel = i)
        g_all <- rbind(g_all, g_rest)
      }
    }
    g_all
  }

  gsea_gs_atc <- build_gsea_atc("twas_gsea", "TWAS-GSEA")
  gsea_gs_atc_nondir <- build_gsea_atc("twas_gsea_nondir", "TWAS-GSEA (non-dir)")

  do.call(rbind, Filter(Negate(is.null), list(magma_gs_atc, gcsc_gs_atc, gsea_gs_atc, gsea_gs_atc_nondir)))
}

#' Build CMAP per-signature drug summary data
#'
#' One row per (cmap_name x cell_iname x pert_itime x pert_idose x weight panel).
#' Reversal_Z is positive when the perturbation opposes the disease TWAS
#' signature (candidate therapeutic). Used directly as the colour aesthetic
#' for the heatmap.
build_cmap_drug_summary_data <- function(gd, gwas) {
  d <- safe_access(gd_read(gd, gwas, "tx/cmap"), "drug")
  if (is.null(d) || nrow(d) == 0) return(NULL)
  d <- as.data.frame(d)
  d$Name    <- paste(d$cmap_name, d$pert_itime, d$pert_idose, sep = ' / ')
  d$Z       <- d$Reversal_Z
  d$FDR_Sig <- !is.na(d$P.FDR) & d$P.FDR < 0.05
  d$Nom_Sig <- !is.na(d$P) & d$P < 0.05
  d$Method  <- 'CMAP'
  d
}

#' Build tissue-specific enrichment summary data
#'
#' Reads MAGMA tissue-specific results (already FDR-adjusted and relabelled in
#' the packaging step). Adds a Retained flag indicating whether the tissue
#' survived the upstream conditional analysis.
#'
#' @param gd A gd_result
#' @param gwas GWAS name
#' @return data.frame ordered by P, or NULL if tissue data absent
build_tissue_data <- function(gd, gwas) {
  spec <- safe_access(gd_read(gd, gwas, "tissue"), "specific")
  if (is.null(spec) || is.null(spec$res) || nrow(spec$res) == 0) return(NULL)
  d <- as.data.frame(spec$res)
  keep <- if (is.null(spec$keep)) character(0) else spec$keep
  d$Retained  <- d$Tissue %in% keep
  d$FDR_Sig   <- d$P.FDR < 0.05
  d$Nom_Sig   <- d$P     < 0.05
  d$negLog10P <- -log10(d$P)
  d[order(d$P), ]
}

#' Build CMAP per-MOA enrichment summary data
build_cmap_moa_summary_data <- function(gd, gwas) {
  d <- safe_access(gd_read(gd, gwas, "tx/cmap"), "moa")
  if (is.null(d) || nrow(d) == 0) return(NULL)
  d <- as.data.frame(d)
  d$Name    <- d$MOA
  # Reversal_Z is positive when the MOA opposes the disease TWAS signature.
  # The MOA Wilcoxon's HL = in - out (opposite to the DrugTargetor ATC HL =
  # out - in), but the sign convention is normalised at the format / read
  # layer, so no flipping is needed here.
  d$Z       <- d$Reversal_Z
  d$FDR_Sig <- !is.na(d$P.FDR) & d$P.FDR < 0.05
  d$Nom_Sig <- !is.na(d$P) & d$P < 0.05
  d$Method  <- 'CMAP'
  d
}
