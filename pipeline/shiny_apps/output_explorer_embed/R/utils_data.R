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

#' Validate RDS structure for GenoDisc results
#'
#' @param rds The object read from readRDS
#' @return Character error message, or NULL if valid
validate_rds <- function(rds) {
  if (!is.list(rds)) {
    return("Uploaded file does not contain a valid GenoDisc results object.")
  }
  if (!("configuration" %in% names(rds))) {
    return("Results object is missing 'configuration'. Was this file produced by the GenoDisc pipeline?")
  }
  conf <- rds$configuration
  if (is.null(conf$config) || is.null(conf$gwas_list)) {
    return("Results object has an incomplete 'configuration' section.")
  }
  gwas_names <- setdiff(names(rds), "configuration")
  if (length(gwas_names) == 0) {
    return("Results object contains no GWAS results.")
  }
  NULL
}

#' Build molecular association summary data
#'
#' @param gwas_results The per-GWAS results list (e.g., rds[["MILEAGE01"]])
#' @param cf Named list of config flags from parse_config_flags()
#' @return data.frame with columns: Panel, ID, Z, Sig, Coloc, Method, Type
build_mol_assoc_data <- function(gwas_results, cf) {
  all_func_res <- NULL

  if (cf$finemap) {
    finemap_ids <- safe_access(gwas_results, "mol_assoc", "finemap", "L1")
    if (!is.null(finemap_ids) && length(finemap_ids) > 0) {
      all_func_res <- rbind(all_func_res, data.frame(
        Panel = "SuSie (L=1)", ID = finemap_ids, Z = 1, Sig = F, Coloc = F,
        Method = "SNP\nFine-mapping", Type = ''))
    }
  }

  if (cf$twas) {
    fusion_res <- safe_access(gwas_results, "mol_assoc", "exp", "fusion", "res")
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
    smr_res <- safe_access(gwas_results, "mol_assoc", "exp", "smr", "results")
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
    pwas_res <- safe_access(gwas_results, "mol_assoc", "protein", "fusion", "results")
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
    smr_prot <- safe_access(gwas_results, "mol_assoc", "protein", "smr", "results")
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
    magma <- safe_access(gwas_results, "mol_assoc", "magma")
    if (!is.null(magma)) {
      all_func_res <- rbind(all_func_res, data.frame(
        Panel = 'MAGMA', ID = magma$ID,
        Z = abs(qnorm(as.numeric(magma$P))),
        Sig = as.numeric(magma$P.FDR) < 0.05,
        Coloc = F, Method = 'MAGMA', Type = ''))
    }
  }

  if (cf$clump) {
    nearest <- safe_access(gwas_results, "mol_assoc", "nearest", "clump")
    if (!is.null(nearest) && length(nearest) > 0) {
      all_func_res <- rbind(all_func_res, data.frame(
        Panel = 'NearestGene', ID = nearest, Z = 1, Sig = F, Coloc = F,
        Method = 'Nearest\nGene', Type = ''))
    }
  }

  all_func_res
}

#' Build drug enrichment summary data
#'
#' @param gwas_results The per-GWAS results list
#' @return data.frame with columns: Name, Z, P, P.FDR, ATC Code, Method, Panel
build_drug_summary_data <- function(gwas_results) {
  magma_gs <- safe_access(gwas_results, "tx", "drug", "magma")
  if (!is.null(magma_gs)) {
    magma_gs$Z <- -qnorm(magma_gs$P)
    magma_gs <- magma_gs[, c('Name', 'Z', 'P', 'P.FDR', 'ATC Code')]
    magma_gs$Method <- 'MAGMA'
    magma_gs$Panel <- 'MAGMA'
  }

  gcsc_gs <- safe_access(gwas_results, "tx", "drug", "gcsc")
  if (!is.null(gcsc_gs)) {
    gcsc_gs <- gcsc_gs[, c('Name', 'Z', 'P', 'P.FDR', 'ATC Code')]
    gcsc_gs$Method <- 'GCSC'
    gcsc_gs$Panel <- 'Brain and Blood'
  }

  gsea_gs <- safe_access(gwas_results, "tx", "drug", "twas_gsea")
  if (!is.null(gsea_gs)) {
    gsea_gs$Method <- 'TWAS-GSEA'
    gsea_gs$Z <- gsea_gs$Estimate / gsea_gs$SE
    gsea_gs <- gsea_gs[, c('Name', 'Z', 'P', 'P.FDR', 'Method', 'Panel', 'ATC Code')]
    gsea_gs$Z <- -gsea_gs$Z

    gsea_gs_all <- gsea_gs
    for (i in unique(gsea_gs_all$Panel)) {
      gsea_gs_i <- gsea_gs[gsea_gs$Panel == i, ]
      gsea_gs_other <- gsea_gs[gsea_gs$Panel != i, ]
      gsea_gs_rest <- data.frame(
        Name = unique(gsea_gs_other$Name[!(gsea_gs_other$Name %in% gsea_gs_i$Name)]),
        Z = NA, P = NA, P.FDR = NA, Method = 'TWAS-GSEA', Panel = i, ATC_Code = NA)
      names(gsea_gs_rest) <- gsub('ATC_Code', 'ATC Code', names(gsea_gs_rest))
      gsea_gs_all <- rbind(gsea_gs_all, gsea_gs_rest)
    }
    gsea_gs <- gsea_gs_all
  }

  do.call(rbind, Filter(Negate(is.null), list(magma_gs, gcsc_gs, gsea_gs)))
}

#' Build ATC enrichment summary data
#'
#' @param gwas_results The per-GWAS results list
#' @return data.frame with columns: Name, Z, FDR_Sig, Nom_Sig, Method, Panel
build_atc_summary_data <- function(gwas_results) {
  magma_gs_atc <- safe_access(gwas_results, "tx", "atc", "magma")
  if (!is.null(magma_gs_atc)) {
    magma_gs_atc$Z <- -qnorm(magma_gs_atc$P)
    magma_gs_atc$FDR_Sig <- magma_gs_atc$P.FDR < 0.05
    magma_gs_atc$Nom_Sig <- magma_gs_atc$P < 0.05
    magma_gs_atc$Name <- paste0(magma_gs_atc$`ATC Code`, ': ', magma_gs_atc$`ATC Description`)
    magma_gs_atc$Method <- 'MAGMA'
    magma_gs_atc$Panel <- 'MAGMA'
    magma_gs_atc <- magma_gs_atc[, c("Name", "Z", "FDR_Sig", "Nom_Sig", "Method", "Panel"), with = F]
  }

  gcsc_gs_atc <- safe_access(gwas_results, "tx", "atc", "gcsc")
  if (!is.null(gcsc_gs_atc)) {
    gcsc_gs_atc$Z <- -qnorm(gcsc_gs_atc$P)
    gcsc_gs_atc$FDR_Sig <- gcsc_gs_atc$P.FDR < 0.05
    gcsc_gs_atc$Nom_Sig <- gcsc_gs_atc$P < 0.05
    gcsc_gs_atc$Name <- paste0(gcsc_gs_atc$`ATC Code`, ': ', gcsc_gs_atc$`ATC Description`)
    gcsc_gs_atc$Method <- 'GCSC'
    gcsc_gs_atc$Panel <- 'GCSC'
    gcsc_gs_atc <- gcsc_gs_atc[, c("Name", "Z", "FDR_Sig", "Nom_Sig", "Method", "Panel"), with = F]
  }

  gsea_gs_atc <- safe_access(gwas_results, "tx", "atc", "twas_gsea")
  if (!is.null(gsea_gs_atc)) {
    gsea_gs_atc$P.FDR_all <- p.adjust(gsea_gs_atc$P, method = 'fdr')
    gsea_gs_atc$P.FDR.onside_all <- p.adjust(gsea_gs_atc$P.oneside, method = 'fdr')
    gsea_gs_atc$Z <- -qnorm(gsea_gs_atc$P)
    gsea_gs_atc$Z <- gsea_gs_atc$Z * sign(gsea_gs_atc$Estimate)
    gsea_gs_atc$FDR_Sig <- gsea_gs_atc$P.FDR_all < 0.05
    gsea_gs_atc$Nom_Sig <- gsea_gs_atc$P < 0.05
    gsea_gs_atc$Name <- paste0(gsea_gs_atc$`ATC Code`, ': ', gsea_gs_atc$`ATC Description`)
    gsea_gs_atc$Method <- 'TWAS-GSEA'
    gsea_gs_atc <- gsea_gs_atc[, c("Name", "Z", "FDR_Sig", "Nom_Sig", "Method", "Panel"), with = F]

    gsea_gs_atc_all <- gsea_gs_atc
    for (i in unique(gsea_gs_atc_all$Panel)) {
      gsea_gs_atc_i <- gsea_gs_atc[gsea_gs_atc$Panel == i, ]
      gsea_gs_atc_other <- gsea_gs_atc[gsea_gs_atc$Panel != i, ]
      gsea_gs_atc_rest <- data.frame(
        Name = unique(gsea_gs_atc_other$Name[!(gsea_gs_atc_other$Name %in% gsea_gs_atc_i$Name)]),
        Z = NA, FDR_Sig = NA, Nom_Sig = NA, Method = 'TWAS-GSEA', Panel = i)
      gsea_gs_atc_all <- rbind(gsea_gs_atc_all, gsea_gs_atc_rest)
    }
    gsea_gs_atc <- gsea_gs_atc_all
  }

  do.call(rbind, Filter(Negate(is.null), list(magma_gs_atc, gcsc_gs_atc, gsea_gs_atc)))
}
