#' Parse config flags from a config character vector
#'
#' @param config Character vector from readLines of config file
#' @return Named list of logicals for each config flag and composite flags
parse_config_flags <- function(config) {

  flag <- function(config, param) {
    match <- config[grepl(paste0('^', param, ':'), config)]
    if (length(match) == 0) return(FALSE)
    match[1] == paste0(param, ": T")
  }

  clump <- flag(config, "clump")
  cojo <- flag(config, "cojo")
  finemap <- flag(config, "finemap")
  ldsc <- flag(config, "ldsc")
  magma_gene <- flag(config, "magma_gene")

  twas_panel_psychencode <- flag(config, "twas_panel_psychencode")
  twas_panel_fusion <- flag(config, "twas_panel_fusion")
  twas <- any(twas_panel_psychencode, twas_panel_fusion)
  twas_conditional <- flag(config, "twas_conditional") & twas == TRUE

  smr_expression_panel_psychencode <- flag(config, "smr_expression_panel_psychencode")
  smr_expression_panel_metabrain_basalganglia <- flag(config, "smr_expression_panel_metabrain_basalganglia")
  smr_expression_panel_metabrain_cerebellum <- flag(config, "smr_expression_panel_metabrain_cerebellum")
  smr_expression_panel_metabrain_cortex <- flag(config, "smr_expression_panel_metabrain_cortex")
  smr_expression_panel_metabrain_hippocampus <- flag(config, "smr_expression_panel_metabrain_hippocampus")
  smr_expression_panel_metabrain_spinalcord <- flag(config, "smr_expression_panel_metabrain_spinalcord")
  smr_expression_panel_eqtlgen <- flag(config, "smr_expression_panel_eqtlgen")

  pwas_panel_rosmap <- flag(config, "pwas_panel_rosmap")
  pwas_panel_banner <- flag(config, "pwas_panel_banner")

  smr_protein_panel_rosmap <- flag(config, "smr_protein_panel_rosmap")

  gcsc <- flag(config, "gcsc")
  magma_drugtargetor <- flag(config, "magma_drugtargetor")
  twas_gsea_lincs <- flag(config, "twas_gsea_lincs")
  twas_so_lincs <- flag(config, "twas_so_lincs")
  dgi_db_comp <- flag(config, "dgi_db_comp")
  twas_gsea_drugtargetor <- flag(config, "twas_gsea_drugtargetor")
  twas_gsea_drugtargetor_nondirectional <- flag(config, "twas_gsea_drugtargetor_nondirectional")
  twas_gsea_cmap <- flag(config, "twas_gsea_cmap")

  # Composite flags
  mol_assoc <- any(magma_gene,
                   twas,
                   smr_expression_panel_psychencode,
                   smr_expression_panel_metabrain_basalganglia,
                   smr_expression_panel_metabrain_cerebellum,
                   smr_expression_panel_metabrain_cortex,
                   smr_expression_panel_metabrain_hippocampus,
                   smr_expression_panel_metabrain_spinalcord,
                   smr_expression_panel_eqtlgen,
                   pwas_panel_rosmap,
                   pwas_panel_banner,
                   smr_protein_panel_rosmap)

  metabrain <- any(smr_expression_panel_metabrain_basalganglia,
                   smr_expression_panel_metabrain_cerebellum,
                   smr_expression_panel_metabrain_cortex,
                   smr_expression_panel_metabrain_hippocampus,
                   smr_expression_panel_metabrain_spinalcord)

  smr_expression <- any(smr_expression_panel_psychencode,
                        metabrain,
                        smr_expression_panel_eqtlgen)

  expression <- any(twas,
                    smr_expression)

  psychencode <- any(twas_panel_psychencode,
                     smr_expression_panel_psychencode)

  protein <- any(pwas_panel_rosmap,
                 pwas_panel_banner,
                 smr_protein_panel_rosmap)

  drug <- any(dgi_db_comp,
              magma_drugtargetor,
              twas_gsea_lincs,
              twas_so_lincs,
              twas_gsea_drugtargetor,
              twas_gsea_drugtargetor_nondirectional,
              twas_gsea_cmap)

  list(
    clump = clump,
    cojo = cojo,
    finemap = finemap,
    ldsc = ldsc,
    magma_gene = magma_gene,
    twas_panel_psychencode = twas_panel_psychencode,
    twas_panel_fusion = twas_panel_fusion,
    twas = twas,
    twas_conditional = twas_conditional,
    smr_expression_panel_psychencode = smr_expression_panel_psychencode,
    smr_expression_panel_metabrain_basalganglia = smr_expression_panel_metabrain_basalganglia,
    smr_expression_panel_metabrain_cerebellum = smr_expression_panel_metabrain_cerebellum,
    smr_expression_panel_metabrain_cortex = smr_expression_panel_metabrain_cortex,
    smr_expression_panel_metabrain_hippocampus = smr_expression_panel_metabrain_hippocampus,
    smr_expression_panel_metabrain_spinalcord = smr_expression_panel_metabrain_spinalcord,
    smr_expression_panel_eqtlgen = smr_expression_panel_eqtlgen,
    pwas_panel_rosmap = pwas_panel_rosmap,
    pwas_panel_banner = pwas_panel_banner,
    smr_protein_panel_rosmap = smr_protein_panel_rosmap,
    gcsc = gcsc,
    magma_drugtargetor = magma_drugtargetor,
    twas_gsea_lincs = twas_gsea_lincs,
    twas_so_lincs = twas_so_lincs,
    dgi_db_comp = dgi_db_comp,
    twas_gsea_drugtargetor = twas_gsea_drugtargetor,
    twas_gsea_drugtargetor_nondirectional = twas_gsea_drugtargetor_nondirectional,
    twas_gsea_cmap = twas_gsea_cmap,
    mol_assoc = mol_assoc,
    metabrain = metabrain,
    smr_expression = smr_expression,
    expression = expression,
    psychencode = psychencode,
    protein = protein,
    drug = drug
  )
}
