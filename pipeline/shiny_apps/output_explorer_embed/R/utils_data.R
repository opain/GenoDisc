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
#' Reads data/gene_positions.rds (built by pipeline/scripts/make_shiny_gene_positions.R) once and
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
  # data.frame() can't recycle length-1 NA_real_ against a length-0 ids
  # vector — construct the empty frame explicitly.
  if (length(ids) == 0) {
    return(data.frame(ID = character(0), CHR = numeric(0), BP = numeric(0),
                      stringsAsFactors = FALSE))
  }
  out <- data.frame(ID = ids, CHR = NA_real_, BP = NA_real_, stringsAsFactors = FALSE)

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
    # Reversal_Z is positive when the drug opposes the trait's TWAS signature
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
build_atc_summary_data <- function(gd, gwas, level = "L3", atc_source = "gls") {
  atc       <- gd_read(gd, gwas, "tx/atc")            # legacy Wilcoxon (+ magma/gcsc slots)
  atc_gls   <- gd_read(gd, gwas, "tx/atc_gls")        # TWAS-GSEA GLS (directional, per panel)
  atc_gls_m <- gd_read(gd, gwas, "tx/atc_gls_magma")  # MAGMA GLS (non-directional, genome-wide)
  atc_vif   <- gd_read(gd, gwas, "tx/atc_vif")        # TWAS-GSEA VIF-OLS (directional, per panel)
  atc_vif_m <- gd_read(gd, gwas, "tx/atc_vif_magma")  # MAGMA VIF-OLS (non-directional, genome-wide)

  # MAGMA (Wilcoxon) / GCSC ATC results exist only at L3; only overlay them on an L3 view.
  show_drugset <- identical(level, "L3")

  # Standardise any ATC block to the summary columns (Name, Z, FDR_Sig, Nom_Sig, Method, Panel).
  # by_level filters to the selected ATC level for blocks that carry a Level column (the GLS blocks);
  # Z uses Reversal_Z when present (signed/magnitude), else -qnorm(P).
  std_atc <- function(g, method_label, panel_label = NULL, by_level = TRUE) {
    if (is.null(g) || nrow(g) == 0) return(NULL)
    g <- data.table::as.data.table(g)
    if (by_level && "Level" %in% names(g)) g <- g[g$Level == level, ]
    if (nrow(g) == 0) return(NULL)
    g$Z <- if ("Reversal_Z" %in% names(g)) g$Reversal_Z else -qnorm(g$P)
    g$Z[!is.finite(g$Z)] <- NA_real_
    g$FDR_Sig <- g$P.FDR < 0.05
    g$Nom_Sig <- g$P < 0.05
    g$Name <- paste0(g$`ATC Code`, ': ', g$`ATC Description`)
    g$Method <- method_label
    if (!is.null(panel_label)) g$Panel <- panel_label
    g[, c("Name", "Z", "FDR_Sig", "Nom_Sig", "Method", "Panel"), with = F]
  }

  # Pad missing class x panel cells so the heatmap grid is complete.
  pad_panels <- function(g) {
    if (is.null(g) || nrow(g) == 0) return(g)
    g_all <- g
    for (i in unique(g$Panel)) {
      g_i <- g[g$Panel == i, ]; g_other <- g[g$Panel != i, ]
      miss <- unique(g_other$Name[!(g_other$Name %in% g_i$Name)])
      if (length(miss) > 0)
        g_all <- rbind(g_all, data.frame(Name = miss, Z = NA, FDR_Sig = NA, Nom_Sig = NA,
                                         Method = g$Method[1], Panel = i))
    }
    g_all
  }

  if (identical(atc_source, "vif")) {
    # Recommended: drug-level VIF-OLS. TWAS-GSEA (directional, per eQTL panel) + MAGMA (non-directional).
    tw  <- pad_panels(std_atc(atc_vif, "TWAS-GSEA (VIF-OLS)"))
    mag <- std_atc(atc_vif_m, "MAGMA (VIF-OLS)", panel_label = "MAGMA")
    out <- do.call(rbind, Filter(Negate(is.null), list(tw, mag)))
  } else if (identical(atc_source, "gls")) {
    # Drug-level GLS (DRUGSETS-style). TWAS-GSEA (directional, per eQTL panel) + MAGMA (non-directional).
    tw  <- pad_panels(std_atc(atc_gls, "TWAS-GSEA (GLS)"))
    mag <- std_atc(atc_gls_m, "MAGMA (GLS)", panel_label = "MAGMA")
    out <- do.call(rbind, Filter(Negate(is.null), list(tw, mag)))
  } else {
    # Legacy: the per-drug Wilcoxon tests (TWAS-GSEA per panel + MAGMA/GCSC, L3 only).
    magma <- if (show_drugset) std_atc(safe_access(atc, "magma"), "MAGMA (Wilcoxon)", panel_label = "MAGMA", by_level = FALSE) else NULL
    gcsc  <- if (show_drugset) std_atc(safe_access(atc, "gcsc"),  "GCSC",             panel_label = "GCSC",  by_level = FALSE) else NULL
    tw    <- pad_panels(std_atc(safe_access(atc, "twas_gsea"),        "TWAS-GSEA (Wilcoxon)",          by_level = FALSE))
    twn   <- pad_panels(std_atc(safe_access(atc, "twas_gsea_nondir"), "TWAS-GSEA (Wilcoxon, non-dir)", by_level = FALSE))
    out <- do.call(rbind, Filter(Negate(is.null), list(magma, gcsc, tw, twn)))
  }
  out
}

#' Is the drug-level GLS (DRUGSETS-style) TWAS-GSEA ATC block present in this bundle?
has_atc_gls <- function(gd, gwas) {
  !is.null(gd_read(gd, gwas, "tx/atc_gls"))
}

#' Is the MAGMA drug-level GLS ATC block present in this bundle?
has_atc_gls_magma <- function(gd, gwas) {
  !is.null(gd_read(gd, gwas, "tx/atc_gls_magma"))
}

#' Is the drug-level VIF-OLS (CAMERA-style) TWAS-GSEA ATC block present in this bundle?
has_atc_vif <- function(gd, gwas) {
  !is.null(gd_read(gd, gwas, "tx/atc_vif"))
}

#' Is the MAGMA drug-level VIF-OLS ATC block present in this bundle?
has_atc_vif_magma <- function(gd, gwas) {
  !is.null(gd_read(gd, gwas, "tx/atc_vif_magma"))
}

#' Is any legacy (Wilcoxon) ATC block present? (tx/atc with >=1 non-null slot).
#' Older bundles carry only this block (no model-based VIF/GLS blocks).
has_atc_legacy <- function(gd, gwas) {
  a <- gd_read(gd, gwas, "tx/atc")
  !is.null(a) && length(Filter(Negate(is.null), a)) > 0
}

#' The default ATC method source for a bundle: the first method actually present
#' (VIF-OLS recommended -> GLS -> legacy Wilcoxon). Used to seed the filter choices
#' and the summary/table reactives so legacy-only bundles still render.
default_atc_source <- function(gd, gwas) {
  if (has_atc_vif(gd, gwas) || has_atc_vif_magma(gd, gwas)) "vif"
  else if (has_atc_gls(gd, gwas) || has_atc_gls_magma(gd, gwas)) "gls"
  else "legacy"
}

#' ATC levels available across the model-based drug-level blocks (VIF-OLS + GLS, both engines).
atc_gls_levels <- function(gd, gwas) {
  lv <- character(0)
  for (blk in c("tx/atc_vif", "tx/atc_vif_magma", "tx/atc_gls", "tx/atc_gls_magma")) {
    g <- gd_read(gd, gwas, blk)
    if (!is.null(g) && "Level" %in% names(g)) lv <- union(lv, unique(g$Level))
  }
  intersect(c("L2","L3","L4"), lv)
}

#' Per-gene evidence tables (membership x TWAS Z x provenance) for the drill-down.
#' Return NULL when the bundle predates the evidence feature.
get_evidence_block <- function(gd, gwas) gd_read(gd, gwas, "tx/evidence")

#' Is the MAGMA per-gene evidence block present in this bundle?
has_magma_evidence <- function(gd, gwas) !is.null(gd_read(gd, gwas, "tx/evidence_magma"))

#' Member-gene MAGMA evidence for one drug (by drug_name): the drug's target genes
#' with their MAGMA gene-level association (ZSTAT/P) + provenance, ordered by ZSTAT.
get_magma_evidence_rows <- function(gd, gwas, drug_name) {
  d <- safe_access(gd_read(gd, gwas, "tx/evidence_magma"), "drug")
  if (is.null(d) || nrow(d) == 0) return(NULL)
  d <- as.data.frame(d)                       # base subset: avoid data.table arg/column scoping
  if (!is.null(drug_name)) d <- d[toupper(d$drug_name) == toupper(drug_name), , drop = FALSE]
  if (nrow(d) == 0) return(NULL)
  d[order(-d$MAGMA.Z), , drop = FALSE]
}

#' Member-gene evidence for one ATC class (code, level, panel) or one drug (name,
#' panel), ordered by |contribution|. `what` is "class" or "drug".
get_evidence_rows <- function(gd, gwas, what = c("class","drug"),
                              id = NULL, panel = NULL, level = NULL) {
  what <- match.arg(what)
  ev <- get_evidence_block(gd, gwas)
  d <- safe_access(ev, what)
  if (is.null(d) || nrow(d) == 0) return(NULL)
  d <- data.table::as.data.table(d)
  if (!is.null(panel)) d <- d[d$Panel == panel, ]
  if (what == "class") {
    if (!is.null(level)) d <- d[d$Level == level, ]
    if (!is.null(id))    d <- d[d$code == id, ]
  } else {
    if (!is.null(id))    d <- d[d$drug_name == id, ]
  }
  if (nrow(d) == 0) return(NULL)
  d[order(-abs(d$contribution)), ]
}

#' Per-drug results for the drugs in one ATC class — the drug-level observations
#' behind a drug-level GLS ATC coefficient. Reads tx/drug (the per-drug DrugTargetor
#' results for the given engine + panel) and keeps drugs whose ATC code falls under
#' the clicked class, ordered by |per-drug T|. The tx/drug ATC code is ATC level 3
#' (4-char); L4 classes are matched at their level-3 parent.
get_atc_class_drugs <- function(gd, gwas, code, level = "L3", panel = NULL,
                                engine = c("twas_gsea","magma")) {
  engine <- match.arg(engine)
  d <- safe_access(gd_read(gd, gwas, "tx/drug"), engine)
  if (is.null(d) || nrow(d) == 0 || !("ATC Code" %in% names(d))) return(NULL)
  d <- data.table::as.data.table(d)
  if (!is.null(panel) && "Panel" %in% names(d)) d <- d[d$Panel == panel, ]
  # tx/drug 'ATC Code' is a ';'-joined list of a drug's L3 codes (multi-ATC drugs carry
  # several). Match if ANY code falls under the clicked class, so multi-code drugs appear
  # in every class they belong to (consistent with the VIF/GLS multi-ATC class explosion).
  k <- if (identical(level, "L2")) 3L else 4L
  tgt <- substr(code, 1, k)
  has_code <- vapply(strsplit(d$`ATC Code`, ";", fixed = TRUE),
                     function(cs) any(substr(cs, 1, k) == tgt), logical(1))
  d <- d[has_code, ]
  if (nrow(d) == 0) return(NULL)
  # Per-drug T from the enrichment coefficient/SE (TWAS uses Estimate; MAGMA uses BETA).
  if (all(c("Estimate","SE") %in% names(d)))   d$T <- d$Estimate / d$SE
  else if (all(c("BETA","SE") %in% names(d)))  d$T <- d$BETA / d$SE
  ord <- if ("T" %in% names(d)) order(-abs(d$T)) else if ("P" %in% names(d)) order(d$P) else seq_len(nrow(d))
  d[ord, ]
}

#' Build CMAP per-signature drug summary data
#'
#' One row per (cmap_name x cell_iname x pert_itime x pert_idose x weight panel).
#' Reversal_Z is positive when the perturbation opposes the trait's TWAS
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
  # Reversal_Z is positive when the MOA opposes the trait's TWAS signature.
  # The MOA Wilcoxon's HL = in - out (opposite to the DrugTargetor ATC HL =
  # out - in), but the sign convention is normalised at the format / read
  # layer, so no flipping is needed here.
  d$Z       <- d$Reversal_Z
  d$FDR_Sig <- !is.na(d$P.FDR) & d$P.FDR < 0.05
  d$Nom_Sig <- !is.na(d$P) & d$P < 0.05
  d$Method  <- 'CMAP'
  d
}
