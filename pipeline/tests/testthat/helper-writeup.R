# Shared setup + fixtures for the write-up generator tests.
# Run with:  Rscript -e "testthat::test_dir('tests/testthat')"  from pipeline/.

# Locate the pipeline directory (the one containing writeup/registry.yaml) by
# walking up from the test file, so tests work from any CWD.
wu_pipeline_dir <- function() {
  cand <- c(
    Sys.getenv("GENODISC_PIPELINE_DIR", unset = NA),
    normalizePath(file.path(getwd()), mustWork = FALSE),
    normalizePath(file.path(getwd(), "..", ".."), mustWork = FALSE),
    normalizePath(file.path(dirname(getwd())), mustWork = FALSE)
  )
  for (d in cand) {
    if (!is.na(d) && file.exists(file.path(d, "writeup", "registry.yaml"))) return(d)
  }
  # final attempt: search upward
  d <- normalizePath(getwd(), mustWork = FALSE)
  for (i in 1:6) {
    if (file.exists(file.path(d, "writeup", "registry.yaml"))) return(d)
    d <- dirname(d)
  }
  stop("could not locate pipeline dir (writeup/registry.yaml)")
}

wu_source <- function() {
  pd <- wu_pipeline_dir()
  for (f in c("format.R", "registry.R", "methods.R", "results.R", "writeup.R")) {
    source(file.path(pd, "scripts", "functions", "writeup", f))
  }
  invisible(pd)
}

# ---- fixtures ---------------------------------------------------------------

fx_twas_res <- function(fdr = c(1e-7, 0.03, 0.5)) {
  data.frame(CHR = c(1, 2, 3),
    `Gene Symbol` = c("G1", "G2", "G3"), PANEL = "Blood (GTEx)",
    TWAS.Z = c(6, -3, 2), TWAS.P = c(1e-9, 2e-3, 0.2), TWAS.P.FDR = fdr,
    COLOC.PP3 = c(0.1, 0.1, 0.1), COLOC.PP4 = c(0.85, 0.2, 0.1),
    COLOC_logical = c(TRUE, FALSE, FALSE), check.names = FALSE, stringsAsFactors = FALSE)
}

fx_panels <- function() data.frame(Type = "Expression", Software = "FUSION",
  Panel = c("Blood (GTEx)", "Brain - Cortex (GTEx)"), N_indiv = c(500, 200),
  N_gene = c(5000, 4000), stringsAsFactors = FALSE)

# A full EUR GWAS exercising h2 (LDSC), clump, MAGMA, TWAS (with sig + coloc),
# plus a present-but-empty SMR block and an absent (not-run) cojo block.
fx_gwas_full <- function() list(
  gwas_qc = list(
    cross_ancestry_weights = FALSE,
    cleaner_dat = list(val = list(n_var_orig = 9000000,
                                  build = list(build = "GRCh37", overlap = 0.99),
                                  n_snp_final = 7000000, lambda_gc = 1.05,
                                  max_chi2 = 120, n_sig_snp = 1234)),
    ldsc_dat = list(val = list(obs_h2_est = 0.066, obs_h2_se = 0.009,
                               int_est = 1.01, int_se = 0.008))
  ),
  snp_assoc = list(
    clump = data.frame(CHR = c(1, 2), BP = c(100, 200), P = c(1e-12, 3e-9),
                       SNP = c("rs1", "rs2"), NearestGene = c("GENEA", "GENEB"),
                       stringsAsFactors = FALSE),
    cojo = NULL
  ),
  mol_assoc = list(
    exp = list(fusion = list(panels = fx_panels(), res = fx_twas_res()),
               smr = list(panels = NULL, results = NULL)),
    protein = list(fusion = NULL, smr = NULL),
    magma = data.frame(CHR = c(1, 2), ID = c("GENE1", "GENE2"),
                       P = c(1e-8, 0.5), P.FDR = c(1e-6, 0.9), stringsAsFactors = FALSE)
  ),
  tx = list(drug = list(), atc = list(), cmap = list(), pathway = list()),
  tissue = list(specific = NULL)
)

# A cross-ancestry GWAS: no LDSC (not produced), a null-result TWAS, no MAGMA.
fx_gwas_null <- function() {
  g <- fx_gwas_full()
  g$gwas_qc$cross_ancestry_weights <- TRUE
  g$gwas_qc$ldsc_dat <- NULL
  g$mol_assoc$exp$fusion$res <- fx_twas_res(fdr = c(0.4, 0.3, 0.5))
  g$mol_assoc$magma <- NULL
  g
}

fx_output <- function() list(TRAIT1 = fx_gwas_full(), TRAIT2 = fx_gwas_null())

fx_config <- function() list(
  info = 0.9, maf = 0.01, maf_diff = 0.2,
  ldsc = "T", clump = "T", cojo = "T", manc_cojo = "F", finemap = "F",
  sbayesrc = "F", sbayesrc_r = "F", sumher = "F",
  gencor_gwas_list = "NA", gencor_within_gwas_list = "F",
  magma_gene = "T", twas_panel_fusion = "T",
  gtex_weights = c("Whole_Blood", "Brain_Cortex"),
  twas_panel_psychencode = "F", twas_panel_fetal_brain = "F", twas_panel_herv = "F",
  external_weights = "F", pwas_panel_rosmap = "F", pwas_panel_banner = "F",
  smr_expression_panel_psychencode = "F", smr_protein_panel_rosmap = "F",
  tissue_magma = "F", magma_drugtargetor = "F")

fx_manifest <- function() list(
  schema_version = 2L, pipeline_version = "1.4.0",
  pipeline_tag = NULL, pipeline_commit = NULL,
  gwas = list(
    list(name = "TRAIT1", label = "Trait One", population = "EUR", n = 100000),
    list(name = "TRAIT2", label = "Trait Two", population = "AFR", n = 50000)))
