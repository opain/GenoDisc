#!/usr/bin/Rscript
library("optparse")

option_list <- list(
  make_option("--config", action = "store", default = NA, type = "character",
              help = "Path to config file [required]")
)

option_list <- c(option_list, list(
  make_option("--pipeline_dir", action="store", default=NA, type="character",
              help="Path to the pipeline directory [required]")
))

opt = parse_args(OptionParser(option_list=option_list))
options(pipeline_dir = opt$pipeline_dir)

library(data.table)
library(jsonlite)
source(file.path(opt$pipeline_dir, 'scripts', 'functions', 'utils_functions.R'))
source_all(file.path(opt$pipeline_dir, 'scripts', 'functions'))
source(file.path(opt$pipeline_dir, 'scripts', 'reader.R'))

# Read in config: merge default config with user config (user takes priority)
library(yaml)
default_config <- read_yaml(file.path(opt$pipeline_dir, 'config.yaml'))
user_config <- read_yaml(opt$config)
merged_config <- default_config
merged_config[names(user_config)] <- user_config

# Convert to readLines-style character vector for downstream compatibility
config <- vapply(names(merged_config), function(k) {
  v <- merged_config[[k]]
  if (is.null(v) || (length(v) == 1 && is.na(v))) {
    paste0(k, ": NA")
  } else if (length(v) > 1) {
    paste0(k, ": [", paste0("'", v, "'", collapse = ","), "]")
  } else {
    paste0(k, ": ", v)
  }
}, character(1), USE.NAMES = FALSE)

# Read in config parameters
gwas_list<-read_param(config = opt$config, param = 'gwas_list', return_obj = T)
outdir <- read_param(config = opt$config, param = 'outdir', return_obj = F)

output<-list()
for(gwas_i in gwas_list$name){

  ##################
  # GWAS QC
  ##################

  gwas_qc<-list()

  gwas_qc$cleaner_dat<-process_cleaner_log(config=opt$config, gwas=gwas_i)
  gwas_qc$ldsc_dat<-process_ldsc_log(config=opt$config, gwas=gwas_i)
  gwas_qc$ldsc_gencor_dat<-process_ldsc_gencor(config=opt$config, gwas=gwas_i)

  # Read MAF plot as base64 (if it exists)
  maf_plot_path <- paste0(outdir, '/results/', gwas_i, '/gwas_sumstat/', gwas_i, '.cleaned.MAF_plot.png')
  if (file.exists(maf_plot_path)) {
    gwas_qc$maf_plot_base64 <- base64enc::base64encode(maf_plot_path)
  } else {
    gwas_qc$maf_plot_base64 <- NULL
  }

  qq_plot_path         <- paste0(outdir, '/results/', gwas_i, '/gwas_sumstat/', gwas_i, '.qq_plot.png')
  manhattan_unlab_path <- paste0(outdir, '/results/', gwas_i, '/gwas_sumstat/', gwas_i, '.manhattan_plot.unlabelled.png')
  manhattan_lab_path   <- paste0(outdir, '/results/', gwas_i, '/gwas_sumstat/', gwas_i, '.manhattan_plot.labelled.png')

  gwas_qc$qq_plot_base64                   <- if (file.exists(qq_plot_path))         base64enc::base64encode(qq_plot_path)         else NULL
  gwas_qc$manhattan_plot_unlabelled_base64 <- if (file.exists(manhattan_unlab_path)) base64enc::base64encode(manhattan_unlab_path) else NULL
  gwas_qc$manhattan_plot_labelled_base64   <- if (file.exists(manhattan_lab_path))   base64enc::base64encode(manhattan_lab_path)   else NULL

  # Raw filtered variant data (variants, offsets, indexes, effective P-threshold)
  # for the Shiny app to rebuild the Manhattan plot on the fly with user-facing
  # knobs (theme, thresholds, gene labels, ...).
  manhattan_data_path <- paste0(outdir, '/results/', gwas_i, '/gwas_sumstat/', gwas_i, '.manhattan_data.rds')
  gwas_qc$manhattan_data <- if (file.exists(manhattan_data_path)) readRDS(manhattan_data_path) else NULL

  ##################
  # SNP associations
  ##################

  snp_assoc<-list()

  if(read_param(config = opt$config, param = 'clump', return_obj = F) == "T"){
    snp_assoc$clump<-fread(paste0(outdir,'/results/',gwas_i,'/clump/',gwas_i,'.GW.clump.clean.csv'))

    locus_plots_dir <- paste0(outdir, '/results/', gwas_i, '/locus_plots/')
    if (dir.exists(locus_plots_dir)) {
      locus_pngs <- list.files(locus_plots_dir, pattern = '\\.png$', full.names = TRUE)
      if (length(locus_pngs) > 0) {
        snp_ids <- sub(paste0('^', gwas_i, '\\.locus_plot\\.'), '', basename(locus_pngs))
        snp_ids <- sub('\\.png$', '', snp_ids)
        snp_assoc$locus_plots <- setNames(
          lapply(locus_pngs, base64enc::base64encode),
          snp_ids
        )
      }
    }
  } else {
    snp_assoc$clump<-NULL
  }

  if(read_param(config = opt$config, param = 'cojo', return_obj = F) == "T"){
    cojo_csv<-paste0(outdir,'/results/',gwas_i,'/cojo/',gwas_i,'.GW.cojo.clean.csv')
    if(file.exists(cojo_csv)){
      snp_assoc$cojo<-fread(cojo_csv)
      snp_assoc$cojo<-snp_assoc$cojo[order(snp_assoc$cojo$CHR, snp_assoc$cojo$BP),]
    } else {
      snp_assoc$cojo<-NULL
    }

    # Carry the per-chromosome COJO status so the app can report which chromosomes (if any)
    # could not be analysed because their independent signals exceeded the LD reference size.
    cojo_status_csv<-paste0(outdir,'/results/',gwas_i,'/cojo/',gwas_i,'.GW.cojo.status.csv')
    if(file.exists(cojo_status_csv)){
      cojo_status_tab<-fread(cojo_status_csv)
      failed_chrs<-cojo_status_tab$chr[cojo_status_tab$status == 'reference_too_small']
      snp_assoc$cojo_status<-list(
        failed_chrs = failed_chrs,
        n_failed    = length(failed_chrs),
        any_failed  = length(failed_chrs) > 0,
        all_failed  = length(failed_chrs) > 0 && (is.null(snp_assoc$cojo) || nrow(snp_assoc$cojo) == 0)
      )
    } else {
      snp_assoc$cojo_status<-NULL
    }
  } else {
    snp_assoc$cojo<-NULL
    snp_assoc$cojo_status<-NULL
  }

  if(read_param(config = opt$config, param = 'finemap', return_obj = F) == "T"){
    snp_assoc$susie<-list()
    snp_assoc$susie$L10<-process_susie(outdir=outdir, gwas=gwas_i, L=10)
    snp_assoc$susie$L1<-process_susie(outdir=outdir, gwas=gwas_i, L=1)
  } else {
    snp_assoc$susie<-NULL
  }

  #################
  # Molecular associations
  #################

  mol_assoc<-list()

  ######
  # Expression
  ######

  mol_assoc$exp<-list()

  ###
  # FUSION
  ###

  mol_assoc$exp$fusion<-read_fusion_exp(config=opt$config, gwas=gwas_i)

  ###
  # SMR
  ###

  mol_assoc$exp$smr<-read_smr_exp(config=opt$config, gwas=gwas_i)

  ######
  # Protein
  ######

  mol_assoc$protein<-list()

  ###
  # FUSION
  ###

  mol_assoc$protein$fusion<-read_fusion_protein(config=opt$config, gwas=gwas_i)

  ###
  # SMR
  ###

  mol_assoc$protein$smr<-read_smr_protein(config=opt$config, gwas=gwas_i)

  ######
  # MAGMA
  ######

  mol_assoc$magma<-read_magma_gene(config=opt$config, gwas=gwas_i)

  ######
  # Nearest
  ######

  mol_assoc$nearest<-list()

  if(read_param(config = opt$config, param = 'clump', return_obj = F) == "T"){
    mol_assoc$nearest$clump<-identify_nearest(snp_assoc$clump$NearestGene)
  }

  if(read_param(config = opt$config, param = 'cojo', return_obj = F) == "T" && !is.null(snp_assoc$cojo) && nrow(snp_assoc$cojo) > 0){
    mol_assoc$nearest$cojo<-identify_nearest(snp_assoc$cojo$NearestGene)
  }

  ######
  # Finemapping
  ######

  mol_assoc$finemap<-list()

  if(read_param(config = opt$config, param = 'finemap', return_obj = F) == "T"){
    mol_assoc$finemap$L1<-unlist(strsplit(snp_assoc$susie$L1$Gene, ', '))
    mol_assoc$finemap$L1<-mol_assoc$finemap$L1[mol_assoc$finemap$L1 != 'None']
    mol_assoc$finemap$L10<-unlist(strsplit(snp_assoc$susie$L10$Gene, ', '))
    mol_assoc$finemap$L10<-mol_assoc$finemap$L10[mol_assoc$finemap$L10 != 'None']
  } else {
    mol_assoc$finemap<-NULL
  }

  #################
  # Drug repruposing
  #################

  tx<-list()

  ######
  # Drug-specific
  ######

  tx$drug<-list()

  ###
  # MAGMA
  ###

  tx$drug$magma<-read_magma_drug(config=opt$config, gwas=gwas_i)

  ###
  # GCSC
  ###

  tx$drug$gcsc <- read_gcsc(config = opt$config, gwas = gwas_i)

  ###
  # TWAS-GSEA
  ###

  tx$drug$twas_gsea<-read_twas_gsea_drug(config=opt$config, gwas=gwas_i, mode='directional')
  tx$drug$twas_gsea_nondir<-read_twas_gsea_drug(config=opt$config, gwas=gwas_i, mode='nondirectional')

  ######
  # ATC
  ######

  tx$atc<-list()

  ###
  # MAGMA
  ###

  tx$atc$magma<-read_magma_drug_atc(config=opt$config, gwas=gwas_i)

  ###
  # GCSC
  ###

  tx$atc$gcsc<-read_gcsc_atc(config=opt$config, gwas=gwas_i)

  ###
  # TWAS-GSEA
  ###

  tx$atc$twas_gsea<-read_twas_gsea_atc(config=opt$config, gwas=gwas_i, mode='directional')
  tx$atc$twas_gsea_nondir<-read_twas_gsea_atc(config=opt$config, gwas=gwas_i, mode='nondirectional')

  ######
  # CMAP TWAS-GSEA (per-signature drug + per-MOA enrichment)
  ######

  tx$cmap<-list()
  tx$cmap$drug<-read_twas_gsea_cmap_drug(config=opt$config, gwas=gwas_i)
  tx$cmap$moa<-read_twas_gsea_cmap_moa(config=opt$config, gwas=gwas_i)

  #################
  # Tissue Enrichment
  #################

  tissue<-list()

  ######
  # Tissue-specific
  ######

  tissue$specific<-read_magma_tissue(config=opt$config, gwas=gwas_i, type='specific')

  ################
  # Package results
  ################

  output[[gwas_i]]<-list( gwas_qc=gwas_qc,
                          snp_assoc=snp_assoc,
                          mol_assoc=mol_assoc,
                          tx=tx,
                          tissue=tissue)

}

################
# Write split bundle
################

# `.gd_block_ids`, `.gd_traverse`, `.gd_row_count` come from reader.R.
# The block IDs there are the canonical set; we walk them per-GWAS and
# write only the blocks whose value is non-NULL.

pkg_dir     <- file.path(outdir, 'results', 'package')
bundle_path <- file.path(outdir, 'results', 'bundle.tar.gz')

# Snakemake declares pkg_dir/manifest.json as the output — the directory
# already exists via `directory()`-style implicit creation, but rerun after
# manual deletion should still succeed. If a previous package exists, clear
# it to avoid stale block files (e.g. a block that was previously present
# but is now empty because config was toggled).
if (dir.exists(pkg_dir)) unlink(pkg_dir, recursive = TRUE)
dir.create(pkg_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(pkg_dir, 'configuration'), showWarnings = FALSE)
dir.create(file.path(pkg_dir, 'gwas'),          showWarnings = FALSE)

# Read pipeline version from VERSION at repo root. Previously derived from
# `git describe --tags` at runtime, which failed under non-owner UIDs
# (CVE-2022-24765) and needed a scoped safe.directory workaround.
version_path <- file.path(dirname(opt$pipeline_dir), 'VERSION')
pipeline_version <- if (file.exists(version_path)) trimws(readLines(version_path, n = 1L)) else NA_character_

# --- Provenance beyond the VERSION string ------------------------------------
# The VERSION string alone can't distinguish two different dev commits that both
# read e.g. "1.1.0". Record enough to make each bundle self-describing and to
# DETECT reference-data drift on a rerun (see the pipeline-versioning work).
#
# Exact commit + tag come from the release worktree's RELEASE.json (written by
# bin/cut_release.sh). For a fallback run off repo/current there is no
# RELEASE.json, so these stay NA - we deliberately do NOT shell out to git (same
# non-owner-UID reason the VERSION read replaced `git describe`).
release_path <- file.path(dirname(opt$pipeline_dir), 'RELEASE.json')
release <- if (file.exists(release_path)) jsonlite::read_json(release_path) else list()
pipeline_commit <- if (!is.null(release$commit))  as.character(release$commit)  else NA_character_
pipeline_tag    <- if (!is.null(release$version)) as.character(release$version) else NA_character_

# Conda env specs the pinned code declares (md5 of each envs/*.yaml). Combined
# with the pinned commit these determine the built environments, since Snakemake
# content-addresses conda envs by their YAML.
env_files <- sort(Sys.glob(file.path(opt$pipeline_dir, 'envs', '*.yaml')))
conda_envs <- if (length(env_files))
  as.list(setNames(unname(tools::md5sum(env_files)), basename(env_files))) else NULL

# Reference-data stamp: pipeline_resources is shared + mutable, so record enough
# to flag drift on a rerun - especially the secondary-GWAS catalog, which grew
# from 3 to 63 traits mid-project and directly changes gencor results.
resource_stamp <- NULL
# Read the raw path from merged_config: `config` is an UNNAMED, "readLines-style"
# character vector of "key: value" lines (vapply USE.NAMES=FALSE above), so
# config[['resdir']] never resolves. merged_config is the parsed yaml list.
resdir <- merged_config[['resdir']]
if (!is.null(resdir) && !is.na(resdir) && nzchar(resdir) && dir.exists(resdir)) {
  lv_path <- file.path(resdir, 'last_version.txt')
  sec_dir <- file.path(resdir, 'data', 'secondary_gwas')
  sec_codes <- if (dir.exists(sec_dir))
    sort(sub('\\.sumstats\\.gz$', '', basename(Sys.glob(file.path(sec_dir, '*.sumstats.gz'))))) else character(0)
  resource_stamp <- list(
    resdir               = resdir,
    last_version         = if (file.exists(lv_path)) trimws(suppressWarnings(readLines(lv_path, n = 1L))) else NA_character_,
    secondary_gwas_count = length(sec_codes),
    secondary_gwas_sha1  = if (length(sec_codes))
      digest::digest(paste(sec_codes, collapse = '\n'), algo = 'sha1', serialize = FALSE) else NA_character_
  )
}

blocks_meta <- setNames(vector('list', length(gwas_list$name)), gwas_list$name)
for (gwas_i in gwas_list$name) {
  gdata <- output[[gwas_i]]
  per <- list()
  for (bid in .gd_block_ids) {
    keys  <- strsplit(bid, '/', fixed = TRUE)[[1L]]
    value <- .gd_traverse(gdata, keys)
    if (is.null(value)) next

    block_file <- file.path(pkg_dir, 'gwas', gwas_i, paste0(bid, '.rds'))
    dir.create(dirname(block_file), recursive = TRUE, showWarnings = FALSE)
    saveRDS(value, file = block_file)

    per[[bid]] <- list(
      bytes = as.integer(file.info(block_file)$size),
      rows  = .gd_row_count(value)
    )
  }
  blocks_meta[[gwas_i]] <- per
}

# Configuration files (character vector + full gwas_list data.frame). These
# stay as .rds because the character vector is what the Shiny app / report
# already consume, and the gwas_list is a data.frame that doesn't round-trip
# through JSON cleanly.
saveRDS(config,    file = file.path(pkg_dir, 'configuration', 'config.rds'))
saveRDS(gwas_list, file = file.path(pkg_dir, 'configuration', 'gwas_list.rds'))

# Copy VERSION and reader.R into the bundle so a downloaded bundle is
# self-contained.
if (file.exists(version_path)) file.copy(version_path, file.path(pkg_dir, 'VERSION'), overwrite = TRUE)
file.copy(file.path(opt$pipeline_dir, 'scripts', 'reader.R'),
          file.path(pkg_dir, 'reader.R'), overwrite = TRUE)

# Manifest — small JSON summary loaded first by the viewer.
gwas_meta <- lapply(seq_len(nrow(gwas_list)), function(i) {
  row <- gwas_list[i, ]
  as.list(row)
})

# schema_version 2 adds pipeline_tag/pipeline_commit/conda_envs/resources. The
# additions are purely additive (older readers ignore unknown keys) and nothing
# validates schema_version today, so this doesn't break existing bundles.
manifest <- list(
  schema_version   = 2L,
  pipeline_version = pipeline_version,
  pipeline_tag     = pipeline_tag,
  pipeline_commit  = pipeline_commit,
  conda_envs       = conda_envs,
  resources        = resource_stamp,
  created_at       = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z", tz = "UTC"),
  config_file      = opt$config,
  gwas             = gwas_meta,
  block_ids        = .gd_block_ids,
  blocks           = blocks_meta
)

# auto_unbox = TRUE turns length-1 vectors into JSON scalars (so
# "pipeline_version": "1.0.0" not ["1.0.0"]). null=null keeps explicit NULLs.
jsonlite::write_json(
  manifest,
  path        = file.path(pkg_dir, 'manifest.json'),
  auto_unbox  = TRUE,
  pretty      = TRUE,
  null        = 'null',
  na          = 'null'
)

################
# Bundle as tarball
################

# tar from parent so the archive contains a `package/` top-level directory.
if (file.exists(bundle_path)) unlink(bundle_path)
tar_status <- system2(
  'tar',
  args = c('-czf', shQuote(bundle_path), '-C', shQuote(dirname(pkg_dir)), basename(pkg_dir))
)
if (!identical(as.integer(tar_status), 0L)) {
  stop("failed to create bundle tarball (tar exit ", tar_status, ")")
}
