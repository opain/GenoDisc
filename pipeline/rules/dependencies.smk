########
# Import required packages
########

import pandas as pd
from pathlib import Path
import multiprocessing
import hashlib
import sys
import tempfile
import os
import subprocess
import re
import glob

######
# Check genodisc conda env is activated
######

conda_env_name = os.getenv('CONDA_DEFAULT_ENV')
if not conda_env_name == 'genodisc':
  print("Error: The genodisc conda environment must be active when running the pipeline.\nFor more information: https://opain.github.io/GenoDisc/pipeline_readme.html")
  sys.exit(1)

######
# Check config file
######

# Check for missing required configuration parameters
required_config_params = ['outdir', 'config_file', 'gwas_list']
missing_or_invalid_config_params = [param for param in required_config_params if param not in config or config[param] == 'NA']

if missing_or_invalid_config_params:
    # Print an informative message
    print(f"Missing or invalid (set to 'NA') required configuration parameters: {', '.join(missing_or_invalid_config_params)}. Please specify these in the configuration file.")

    # Exit Snakemake gracefully
    sys.exit(1)

# Bivariate LDSC precondition (runs after Snakemake config merge)
gencor_gwas_list_val = config.get('gencor_gwas_list', 'NA')
if gencor_gwas_list_val not in (None, 'NA'):
    if config.get('ldsc', 'F') != 'T':
        print("Error: gencor_gwas_list is set but ldsc is not 'T'; genetic correlation requires LDSC heritability - set ldsc: 'T'.")
        sys.exit(1)

# Set outdir parameter
outdir=config['outdir']

# Set resource directory
resdir = config.get('resdir', None)
if resdir is None or resdir == 'NA':
    resdir = 'resources'

# Set chromosomes to analyse
chromosomes = config.get("chromosomes", list(range(1, 23)))

########
# Create required functions
########

def get_current_version():
    out = subprocess.run(
        ["git", "-c", "safe.directory=*", "describe", "--tags"],
        cwd=workflow.basedir,
        capture_output=True, text=True,
    )
    if out.returncode != 0:
        raise ValueError(f"git describe failed: {out.stderr.strip()}")
    tag = out.stdout.strip()
    m = re.match(r"v?(\d+)\.(\d+)", tag)
    if not m:
        raise ValueError(f"Git tag {tag!r} has no valid version format.")
    return int(m.group(1)), int(m.group(2))

def read_last_version():
    if os.path.exists(last_version_file):
        with open(last_version_file, "r") as file:
            major, minor = file.read().strip().split('.')
            return int(major), int(minor)
    return 0, 0  # Default to 0.0 if file does not exist

def write_last_version(major, minor):
    with open(last_version_file, "w") as file:
        file.write(f"{major}.{minor}")

########
# Check for repo version updates
########

# If there has been a change to the major or minor version numbers, we will rerun the entire pipeline

# Define the path for storing the last known version. Stored PER JOB (outdir),
# not in the shared resdir: with per-job pipeline-version pinning (the web app
# pins each submission to a release and reruns reuse it), a job dir always
# re-runs the version it was created with, so a per-job stamp keeps this guard
# correct AND stops it firing spuriously when different pinned versions run
# concurrently against the same shared resdir. (resdir/last_version.txt is left
# alone as a separate resource-set stamp, recorded in the results manifest.)
os.makedirs(resdir, exist_ok=True)
os.makedirs(outdir, exist_ok=True)
last_version_file = f"{outdir}/last_version.txt"

# Access overwrite flag from config
overwrite = config.get("overwrite", "false").lower() == "true"

# Main logic to check version and decide on execution
current_major, current_minor = get_current_version()
last_major, last_minor = read_last_version()

# Check if the last version is 0.0, proceed without requiring overwrite
if last_major == 0 and last_minor == 0:
    print(f"Initial version setup detected. Updating to v{current_major}.{current_minor}.")
    write_last_version(current_major, current_minor)
else:
    # Check for both major and minor version changes
    if current_major != last_major or current_minor != last_minor:
        if not overwrite:
            print(f"Change in version of GenoDisc detected from v{last_major}.{last_minor} to v{current_major}.{current_minor}. Use --config overwrite=true to proceed.")
            sys.exit(1)
        else:
            print("Proceeding with version update due to overwrite=true config.")
            write_last_version(current_major, current_minor)  # Update the stored version

####
# Download BioMart gene annotations
####

rule download_biomart:
  output:
    f"{resdir}/data/biomart/biomart_genes_grch37.tsv",
    f"{resdir}/data/biomart/gene_locations.tsv"
  benchmark:
    f"{resdir}/benchmarks/download_biomart.tsv"
  params:
    resdir=resdir
  conda:
    "../envs/main.yaml"
  log:
    f"{resdir}/logs/download_biomart.log"
  shell:
    "Rscript --vanilla {workflow.basedir}/scripts/download_biomart.R --pipeline_dir {workflow.basedir} --resdir {params.resdir} > {log} 2>&1"

####
# Download recombination rate maps (Pickrell HapMap-interpolated, GRCh37)
####

rule download_recomb_map:
  output:
    f"{resdir}/data/recomb_maps/chr{{chr}}.interpolated_genetic_map.gz"
  benchmark:
    f"{resdir}/benchmarks/download_recomb_map_chr{{chr}}.tsv"
  params:
    resdir=resdir
  log:
    f"{resdir}/logs/download_recomb_map_chr{{chr}}.log"
  shell:
    "(mkdir -p {params.resdir}/data/recomb_maps/; "
    "wget --no-check-certificate -O {output} "
    "https://raw.githubusercontent.com/joepickrell/1000-genomes-genetic-maps/master/interpolated_from_hapmap/chr{wildcards.chr}.interpolated_genetic_map.gz) > {log} 2>&1"

####
# Download liftover and GRCh 37 to 38 and 36 track
####

# Download liftover
rule install_liftover:
  output:
    touch(f"{resdir}/software/install_liftover.done")
  benchmark:
    f"{resdir}/benchmarks/install_liftover.tsv"
  conda:
    "../envs/main.yaml"
  params:
    resdir=resdir
  log:
    f"{resdir}/logs/install_liftover.log"
  shell:
    "(mkdir -p {params.resdir}/software/liftover/; wget --no-check-certificate -O {params.resdir}/software/liftover/liftover https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver) > {log} 2>&1"

# Download liftover track
rule download_liftover_track:
  output:
    touch(f"{resdir}/software/download_liftover_track.done")
  benchmark:
    f"{resdir}/benchmarks/download_liftover_track.tsv"
  conda:
    "../envs/main.yaml"
  params:
    resdir=resdir
  log:
    f"{resdir}/logs/download_liftover_track.log"
  shell:
    "(mkdir -p {params.resdir}/data/liftover/; wget --no-check-certificate -O {params.resdir}/data/liftover/hg19ToHg38.over.chain.gz ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz; wget --no-check-certificate -O {params.resdir}/data/liftover/hg19ToHg18.over.chain.gz ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg18.over.chain.gz) > {log} 2>&1"

####
# Download and format 1000 Genomes reference data
####

rule prep_1kg:
  input:
    rules.install_liftover.output,
    rules.download_liftover_track.output
  resources:
    mem_mb=20000
  output:
    touch(f"{resdir}/data/prep_1kg.done")
  benchmark:
    f"{resdir}/benchmarks/prep_1kg.tsv"
  conda:
    "../envs/main.yaml"
  log:
    f"{resdir}/logs/prep_1kg.log"
  shell:
    "Rscript --vanilla {workflow.basedir}/scripts/prep_1kg.R --pipeline_dir {workflow.basedir} > {log} 2>&1"

####
# Download LDSC
####

# Install LDSC
rule install_ldsc:
  output:
    directory(f"{resdir}/software/ldsc/")
  benchmark:
    f"{resdir}/benchmarks/install_ldsc.tsv"
  conda:
    "../envs/main.yaml"
  log:
    f"{resdir}/logs/install_ldsc.log"
  shell:
    "(git clone https://github.com/bulik/ldsc.git {output}; \
    cd {output}; \
    git reset --hard aa33296abac9569a6422ee6ba7eb4b902422cc74) > {log} 2>&1"

# Download LDSC reference data
rule download_ldsc_scores:
  output:
    f"{resdir}/data/ldsc/eur_w_ld_chr/10.l2.ldscore.gz"
  benchmark:
    f"{resdir}/benchmarks/download_ldsc_scores.tsv"
  conda:
    "../envs/main.yaml"
  params:
    resdir=resdir
  log:
    f"{resdir}/logs/download_ldsc_scores.log"
  shell:
    "(mkdir -p {params.resdir}/data/ldsc; wget --no-check-certificate -O {params.resdir}/data/ldsc/eur_w_ld_chr.tar.gz https://zenodo.org/record/8182036/files/eur_w_ld_chr.tar.gz?download=1; tar -xf {params.resdir}/data/ldsc/eur_w_ld_chr.tar.gz -C {params.resdir}/data/ldsc; rm {params.resdir}/data/ldsc/eur_w_ld_chr.tar.gz) > {log} 2>&1"

rule download_ldsc_hm3:
  output:
    f"{resdir}/data/ldsc/w_hm3.snplist"
  benchmark:
    f"{resdir}/benchmarks/download_ldsc_hm3.tsv"
  conda:
    "../envs/main.yaml"
  params:
    resdir=resdir
  log:
    f"{resdir}/logs/download_ldsc_hm3.log"
  shell:
    "(mkdir -p {params.resdir}/data/ldsc; wget --no-check-certificate -O {params.resdir}/data/ldsc/w_hm3.snplist.gz https://zenodo.org/record/7773502/files/w_hm3.snplist.gz?download=1; gzip -d {params.resdir}/data/ldsc/w_hm3.snplist.gz) > {log} 2>&1"

# Download GCTA
rule download_gcta:
  output:
    directory(f'{resdir}/software/gcta/gcta_v1.94.0Beta_linux_kernel_3_x86_64/')
  benchmark:
    f"{resdir}/benchmarks/download_gcta.tsv"
  conda:
    "../envs/main.yaml"
  params:
    resdir=resdir
  log:
    f"{resdir}/logs/download_gcta.log"
  shell:
    "(mkdir -p {params.resdir}/software/gcta; wget -O {params.resdir}/software/gcta/gcta_v1.94.0Beta_linux_kernel_3_x86_64.zip https://yanglab.westlake.edu.cn/software/gcta/bin/gcta_v1.94.0Beta_linux_kernel_3_x86_64.zip; unzip {params.resdir}/software/gcta/gcta_v1.94.0Beta_linux_kernel_3_x86_64.zip -d {params.resdir}/software/gcta; rm {params.resdir}/software/gcta/gcta_v1.94.0Beta_linux_kernel_3_x86_64.zip) > {log} 2>&1"

# Install GenoUtils
rule install_genoutils:
  input:
    f"{workflow.basedir}/envs/main.yaml"
  output:
    touch(f"{resdir}/software/install_genoutils.done")
  benchmark:
    f"{resdir}/benchmarks/install_genoutils.tsv"
  conda:
    "../envs/main.yaml"
  log:
    f"{resdir}/logs/install_genoutils.log"
  shell:
    "Rscript --vanilla -e 'devtools::install_github(\"opain/GenoUtils@fb6c1f4711b37451a078ac0438714f14db693dbd\", upgrade = \"never\")' > {log} 2>&1"

####
# Download MAGMA
####

rule download_magma:
  output:
    f"{resdir}/software/magma/magma"
  benchmark:
    f"{resdir}/benchmarks/download_magma.tsv"
  conda:
    "../envs/main.yaml"
  params:
    resdir=resdir
  log:
    f"{resdir}/logs/download_magma.log"
  shell:
    "(rm -rf {params.resdir}/software/magma; \
    wget -O {params.resdir}/software/magma.zip https://vu.data.surfsara.nl/index.php/s/zkKbNeNOZAhFXZB/download; \
    unzip {params.resdir}/software/magma.zip -d {params.resdir}/software/magma; \
    rm {params.resdir}/software/magma.zip) > {log} 2>&1"

####
# Download MAGMA gene locations
####

rule download_magma_gene_loc:
  output:
    f"{resdir}/data/magma/NCBI37.3.gene.loc"
  benchmark:
    f"{resdir}/benchmarks/download_magma_gene_loc.tsv"
  conda:
    "../envs/main.yaml"
  params:
    resdir=resdir
  log:
    f"{resdir}/logs/download_magma_gene_loc.log"
  shell:
    "(rm -rf {params.resdir}/data/magma; \
    wget -O {params.resdir}/data/magma.zip https://vu.data.surfsara.nl/index.php/s/Pj2orwuF2JYyKxq/download; \
    unzip {params.resdir}/data/magma.zip -d {params.resdir}/data/magma; \
    rm {params.resdir}/data/magma.zip) > {log} 2>&1"

####
# Download MAGMA reference
####

rule download_magma_ref:
  output:
    f"{resdir}/data/magma_ref/g1000_eur.bed"
  benchmark:
    f"{resdir}/benchmarks/download_magma_ref.tsv"
  conda:
    "../envs/main.yaml"
  params:
    resdir=resdir
  log:
    f"{resdir}/logs/download_magma_ref.log"
  shell:
    "(rm -rf {params.resdir}/data/magma_ref; \
    wget -O {params.resdir}/data/magma.zip https://vu.data.surfsara.nl/index.php/s/VZNByNwpD8qqINe/download; \
    unzip {params.resdir}/data/magma.zip -d {params.resdir}/data/magma_ref; \
    rm {params.resdir}/data/magma.zip) > {log} 2>&1"

####
# Create MAGMA annotation file
####

rule magma_annot:
  input:
    rules.download_magma.output,
    rules.download_magma_gene_loc.output,
    rules.download_magma_ref.output
  output:
    f"{resdir}/data/magma/NCBI37.3.genes.annot"
  benchmark:
    f"{resdir}/benchmarks/magma_annot.tsv"
  conda:
    "../envs/main.yaml"
  params:
    resdir=resdir
  log:
    f"{resdir}/logs/magma_annot.log"
  shell:
    "{params.resdir}/software/magma/magma \
      --annotate window=35,10 \
    	--snp-loc {params.resdir}/data/magma_ref/g1000_eur.bim \
    	--gene-loc {params.resdir}/data/magma/NCBI37.3.gene.loc \
    	--out {params.resdir}/data/magma/NCBI37.3 > {log} 2>&1"

####
# Download ATC codes
####

rule download_atc:
  output:
    f"{resdir}/data/atc/atc_20220201.txt"
  benchmark:
    f"{resdir}/benchmarks/download_atc.tsv"
  conda:
    "../envs/main.yaml"
  params:
    resdir=resdir
  log:
    f"{resdir}/logs/download_atc.log"
  shell:
    "(rm -rf {params.resdir}/data/atc; \
    wget -O {params.resdir}/data/2022-02-01-v3extracts.zip https://www.pbs.gov.au/downloads/2022/02/2022-02-01-v3extracts.zip; \
    mkdir -p {params.resdir}/data/atc; \
    unzip {params.resdir}/data/2022-02-01-v3extracts.zip -d {params.resdir}/data/atc; \
    rm {params.resdir}/data/2022-02-01-v3extracts.zip) > {log} 2>&1"

####
# Download and format DrugTargetor database
####

rule download_drug_targetor:
  output:
    f"{resdir}/data/drug_targetor/wholedatabase_for_targetor"
  benchmark:
    f"{resdir}/benchmarks/download_drug_targetor.tsv"
  conda:
    "../envs/main.yaml"
  params:
    resdir=resdir
  log:
    f"{resdir}/logs/download_drug_targetor.log"
  shell:
    "(mkdir -p {params.resdir}/data/drug_targetor/; \
    wget -O {params.resdir}/data/drug_targetor/wholedatabase_for_targetor https://github.com/hagax8/drugtargetor/raw/master/wholedatabase_for_targetor) > {log} 2>&1"

rule format_drug_targetor:
  input:
    rules.download_drug_targetor.output,
    rules.download_magma_gene_loc.output
  output:
    f"{resdir}/data/drug_targetor/wholedatabase_for_targetor.gmt",
    f"{resdir}/data/drug_targetor/wholedatabase_for_targetor_symbols.gmt"
  benchmark:
    f"{resdir}/benchmarks/format_drug_targetor.tsv"
  conda:
    "../envs/main.yaml"
  log:
    f"{resdir}/logs/format_drug_targetor.log"
  shell:
    "Rscript --vanilla {workflow.basedir}/scripts/format_drug_targetor.R --pipeline_dir {workflow.basedir} > {log} 2>&1"

####
# Download and format GTEx TPM data
####

rule prep_tissue_exp:
  input:
    rules.download_magma_gene_loc.output,
    f"{workflow.basedir}/scripts/prep_tissue_exp.R"
  output:
    f"{resdir}/data/gtex/GTEx_v8_group.tsv"
  benchmark:
    f"{resdir}/benchmarks/prep_tissue_exp.tsv"
  conda:
    "../envs/main.yaml"
  log:
    f"{resdir}/logs/prep_tissue_exp.log"
  shell:
    "Rscript --vanilla {workflow.basedir}/scripts/prep_tissue_exp.R --pipeline_dir {workflow.basedir} > {log} 2>&1"

####
# Download software required for TWAS-related analysis
####

# Install fusion
rule install_fusion:
  output:
    directory(f"{resdir}/software/fusion/")
  benchmark:
    f"{resdir}/benchmarks/install_fusion.tsv"
  conda:
    "../envs/main.yaml"
  log:
    f"{resdir}/logs/install_fusion.log"
  shell:
    "(git clone https://github.com/opain/fusion_twas.git {output}; \
    cd {output}; \
    git reset --hard 4635dd1aeafabafd5c062d2e9002e37e129a043d) > {log} 2>&1"

# Download plink2R
rule download_plink2R:
  input:
    rules.install_fusion.output
  output:
    f"{resdir}/software/plink2R/plink2R-master/data.bed"
  benchmark:
    f"{resdir}/benchmarks/download_plink2R.tsv"
  conda:
    "../envs/main.yaml"
  params:
    resdir=resdir
  log:
    f"{resdir}/logs/download_plink2R.log"
  shell:
    "(rm -rf {params.resdir}/software/plink2R; \
    mkdir -p {params.resdir}/software/plink2R; \
    wget -O {params.resdir}/software/plink2R/master.zip https://github.com/gabraham/plink2R/archive/master.zip; \
    unzip {params.resdir}/software/plink2R/master.zip -d {params.resdir}/software/plink2R) > {log} 2>&1"

# Install plink2R
rule install_plink2R:
  input:
    rules.download_plink2R.output,
    f"{workflow.basedir}/envs/main.yaml"
  output:
    touch(f"{resdir}/software/install_plink2R")
  benchmark:
    f"{resdir}/benchmarks/install_plink2R.tsv"
  conda:
    "../envs/main.yaml"
  params:
    resdir=resdir
  log:
    f"{resdir}/logs/install_plink2R.log"
  shell:
    "Rscript --vanilla -e 'install.packages(\"{params.resdir}/software/plink2R/plink2R-master/plink2R/\",repos=NULL)' > {log} 2>&1"

# Install SNP-weights pipeline repo
rule install_snp_weight_pipe:
  output:
    directory(f"{resdir}/software/Calculating-FUSION-TWAS-weights-pipeline/")
  benchmark:
    f"{resdir}/benchmarks/install_snp_weight_pipe.tsv"
  conda:
    "../envs/main.yaml"
  log:
    f"{resdir}/logs/install_snp_weight_pipe.log"
  shell:
    "(git clone https://github.com/opain/Calculating-FUSION-TWAS-weights-pipeline.git {output}; \
    cd {output}; \
    git reset --hard ab15a41e4568107f29bc5a538ea016a554d58589) > {log} 2>&1"

####
# Download data for TWAS related analysis
####

# Download PsychENCODE SNP-weights
rule download_psychENCODE_weights:
  output:
    touch(f"{resdir}/data/download_psychENCODE_weights.done")
  benchmark:
    f"{resdir}/benchmarks/download_psychENCODE_weights.tsv"
  conda:
    "../envs/main.yaml"
  params:
    resdir=resdir
  log:
    f"{resdir}/logs/download_psychENCODE_weights.log"
  shell:
    "(mkdir -p {params.resdir}/data/fusion_snp_weights/psychencode; \
    wget -O {params.resdir}/data/fusion_snp_weights/psychencode/PEC_TWAS_weights.tar.gz http://resource.psychencode.org/Datasets/Derived/PEC_TWAS_weights.tar.gz; \
    mkdir -p {params.resdir}/data/fusion_snp_weights/psychencode/psychencode; \
    tar xvzf {params.resdir}/data/fusion_snp_weights/psychencode/PEC_TWAS_weights.tar.gz -C {params.resdir}/data/fusion_snp_weights/psychencode/psychencode; \
    rm {params.resdir}/data/fusion_snp_weights/psychencode/PEC_TWAS_weights.tar.gz) > {log} 2>&1"

# Format PsychENCODE SNP-weights
rule format_psychencode:
  input:
    psychencode_data=rules.download_psychENCODE_weights.output,
    weights_pipe=rules.install_snp_weight_pipe.output,
    biomart=rules.download_biomart.output
  output:
    f"{resdir}/data/format_psychencode.done"
  benchmark:
    f"{resdir}/benchmarks/format_psychencode.tsv"
  conda:
    "../envs/main.yaml"
  log:
    f"{resdir}/logs/format_psychencode.log"
  shell:
    "Rscript --vanilla {workflow.basedir}/scripts/format_psychENCODE.R --pipeline_dir {workflow.basedir} > {log} 2>&1"

# Download FUSION GTEx v8 EUR SNP-weights
# I am using EUR instead of full sample to avoid LD mismatch
gtex_weights=config["gtex_weights"]

rule download_gtex_weights:
  output:
    touch(f"{resdir}/data/download_fusion_gtex_{{weight}}_weights.done")
  benchmark:
    f"{resdir}/benchmarks/download_gtex_weights_{{weight}}.tsv"
  conda:
    "../envs/main.yaml"
  params:
    resdir=resdir
  log:
    f"{resdir}/logs/download_gtex_weights-{{weight}}.log"
  shell:
    "(mkdir -p {params.resdir}/data/fusion_snp_weights/{wildcards.weight}; wget -O {params.resdir}/data/fusion_snp_weights/GTExv8.EUR.{wildcards.weight}.tar.gz https://s3.us-west-1.amazonaws.com/gtex.v8.fusion/EUR/GTExv8.EUR.{wildcards.weight}.tar.gz; tar xf {params.resdir}/data/fusion_snp_weights/GTExv8.EUR.{wildcards.weight}.tar.gz -C {params.resdir}/data/fusion_snp_weights/{wildcards.weight}; rm {params.resdir}/data/fusion_snp_weights/GTExv8.EUR.{wildcards.weight}.tar.gz; mv {params.resdir}/data/fusion_snp_weights/{wildcards.weight}/GTExv8.EUR.{wildcards.weight}.pos {params.resdir}/data/fusion_snp_weights/{wildcards.weight}/{wildcards.weight}.pos; mv {params.resdir}/data/fusion_snp_weights/{wildcards.weight}/GTExv8.EUR.{wildcards.weight} {params.resdir}/data/fusion_snp_weights/{wildcards.weight}/{wildcards.weight}) > {log} 2>&1"

# Update GTEx v8 P0 and P1 to build GRCh 37
rule update_gtex_coord:
  input:
    f"{resdir}/data/download_fusion_gtex_{{weight}}_weights.done",
    rules.download_biomart.output
  output:
    touch(f"{resdir}/data/update_gtex_coord_{{weight}}.done")
  benchmark:
    f"{resdir}/benchmarks/update_gtex_coord_{{weight}}.tsv"
  conda:
    "../envs/main.yaml"
  log:
    f"{resdir}/logs/update_gtex_coord-{{weight}}.log"
  shell:
    "Rscript --vanilla {workflow.basedir}/scripts/update_gtex_coord.R --pipeline_dir {workflow.basedir} \
      --panel {wildcards.weight} > {log} 2>&1"

rule update_gtex_coord_all_panel:
    input: expand(f"{resdir}/data/update_gtex_coord_{{weight}}.done", weight=gtex_weights)

# Download FUSION non-GTEx SNP-weights
non_gtex_weights=config["non_gtex_weights"]

rule download_non_gtex_weights:
  output:
    touch(f"{resdir}/data/download_non_gtex_{{weight}}_weights.done")
  benchmark:
    f"{resdir}/benchmarks/download_non_gtex_weights_{{weight}}.tsv"
  conda:
    "../envs/main.yaml"
  params:
    resdir=resdir
  log:
    f"{resdir}/logs/download_non_gtex_weights-{{weight}}.log"
  shell:
    "(mkdir -p {params.resdir}/data/fusion_snp_weights/{wildcards.weight}; wget --no-check-certificate -O {params.resdir}/data/fusion_snp_weights/{wildcards.weight}.tar.bz2 https://data.broadinstitute.org/alkesgroup/FUSION/WGT/{wildcards.weight}.tar.bz2; tar xvjf {params.resdir}/data/fusion_snp_weights/{wildcards.weight}.tar.bz2 -C {params.resdir}/data/fusion_snp_weights/{wildcards.weight}; rm {params.resdir}/data/fusion_snp_weights/{wildcards.weight}.tar.bz2) > {log} 2>&1"

# Insert N into non-GTEX SNP-weights
rule insert_n_nongtex:
  input:
    f"{resdir}/data/download_non_gtex_{{weight}}_weights.done",
    rules.download_biomart.output
  output:
    touch(f"{resdir}/data/insert_n_nongtex_{{weight}}.done")
  benchmark:
    f"{resdir}/benchmarks/insert_n_nongtex_{{weight}}.tsv"
  conda:
    "../envs/main.yaml"
  log:
    f"{resdir}/logs/insert_n_nongtex-{{weight}}.log"
  shell:
    "Rscript --vanilla {workflow.basedir}/scripts/insert_n_nongtex.R --pipeline_dir {workflow.basedir} \
      --panel {wildcards.weight} > {log} 2>&1"

rule insert_n_nongtex_all_panel:
    input: expand(f"{resdir}/data/insert_n_nongtex_{{weight}}.done", weight=non_gtex_weights)

# Download glist file
rule download_glist:
  output:
    f"{resdir}/data/glist-hg19"
  benchmark:
    f"{resdir}/benchmarks/download_glist.tsv"
  conda:
    "../envs/main.yaml"
  params:
    resdir=resdir
  log:
    f"{resdir}/logs/download_glist.log"
  shell:
    "wget -P {params.resdir}/data/ https://www.cog-genomics.org/static/bin/plink/glist-hg19 > {log} 2>&1"

####
# Download TWAS-GSEA
####

rule install_twas_gsea:
  output:
    directory(f"{resdir}/software/TWAS-GSEA/")
  benchmark:
    f"{resdir}/benchmarks/install_twas_gsea.tsv"
  conda:
    "../envs/main.yaml"
  log:
    f"{resdir}/logs/install_twas_gsea.log"
  shell:
    "(git clone https://github.com/opain/TWAS-GSEA.git {output}; \
    cd {output}; \
    git checkout optimisation; \
    git reset --hard b47a3a36375b4420645de810c0c5ef782c8ad44d) > {log} 2>&1"

####
# Download FeaturePred
####

rule install_feature_pred:
  output:
    directory(f"{resdir}/software/Predicting-TWAS-features/")
  benchmark:
    f"{resdir}/benchmarks/install_feature_pred.tsv"
  conda:
    "../envs/main.yaml"
  log:
    f"{resdir}/logs/install_feature_pred.log"
  shell:
    "(git clone https://github.com/opain/Predicting-TWAS-features.git {output}; \
    cd {output}; \
    git reset --hard b9defcf3c96145ab86f605629c48e0d29daebe0c) > {log} 2>&1"

####
# Download pigz
####

rule install_pigz:
  output:
    f"{resdir}/software/pigz/pigz/pigz"
  benchmark:
    f"{resdir}/benchmarks/install_pigz.tsv"
  conda:
    "../envs/main.yaml"
  params:
    resdir=resdir
  log:
    f"{resdir}/logs/install_pigz.log"
  shell:
    "(wget -O {params.resdir}/software/pigz.tar.gz https://zlib.net/pigz/pigz.tar.gz; mkdir -p {params.resdir}/software/pigz; tar xvzf {params.resdir}/software/pigz.tar.gz -C {params.resdir}/software/pigz; rm {params.resdir}/software/pigz.tar.gz; cd {params.resdir}/software/pigz/pigz; make) > {log} 2>&1"

####
# Format the external SNP-weights for TWAS
####

if config["external_weights"] == "T":
  external_weights_list=config["external_weights_pos_path"]
  import os.path
  external_weights_path_list=[os.path.dirname(x) for x in external_weights_list]
  external_weights_id_list=[os.path.basename(x) for x in external_weights_list]
  external_weights_id_list=[re.sub(".pos", "", x) for x in external_weights_id_list]

  import os
  for x in list(range(0, len(external_weights_path_list))):
    if not os.path.isdir("".join([resdir, "/data/fusion_snp_weights/",external_weights_id_list[x]])):
      os.system("".join(["mkdir ", resdir, "/data/fusion_snp_weights/",external_weights_id_list[x]]))
      os.system("".join(["cp -r ",external_weights_path_list[x],"/* ", resdir, "/data/fusion_snp_weights/", external_weights_id_list[x],"/"]))

####
# Predict features into 1kg sample
####

# Make complete list of panels without Splicing
weights=gtex_weights + non_gtex_weights
if config["twas_panel_psychencode"] == "T":
  weights.append("psychencode")

if config["external_weights"] == "T":
  weights=weights + external_weights_id_list

import copy
weights_nosplice=copy.copy(weights)
if "CMC.BRAIN.RNASEQ_SPLICING" in weights_nosplice:
    weights_nosplice.remove("CMC.BRAIN.RNASEQ_SPLICING")

def feature_pred_input(wildcards):
    inputs = [
        f"{resdir}/software/Predicting-TWAS-features/",
        f"{resdir}/software/pigz/pigz/pigz"
    ]
    w = wildcards.weight
    if w == "psychencode":
        inputs.append(f"{resdir}/data/format_psychencode.done")
    elif w in gtex_weights:
        inputs.append(f"{resdir}/data/update_gtex_coord_{w}.done")
    elif w in non_gtex_weights:
        inputs.append(f"{resdir}/data/insert_n_nongtex_{w}.done")
    return inputs

# Modify panel column in .pos file
rule feature_pred:
  resources:
    mem_mb=50000,
    cpus=5
  input:
    feature_pred_input
  output:
    f"{resdir}/data/predicted_expression/{{weight}}/Reference_Expression/Reference_Expression_{{weight}}.txt.gz"
  benchmark:
    f"{resdir}/benchmarks/feature_pred_{{weight}}.tsv"
  conda:
    "../envs/main.yaml"
  params:
    resdir=resdir
  log:
    f"{resdir}/logs/feature_pred-{{weight}}.log"
  shell:
    "Rscript --vanilla {params.resdir}/software/Predicting-TWAS-features/FeaturePred.V2.0.R \
    	--PLINK_prefix_chr {params.resdir}/data/1kg/1KG.Phase3.EUR.MAF_001.chr \
    	--weights {params.resdir}/data/fusion_snp_weights/{wildcards.weight}/{wildcards.weight}.pos \
    	--weights_dir {params.resdir}/data/fusion_snp_weights/{wildcards.weight} \
    	--ref_ld_chr {params.resdir}/data/1kg/1KG.Phase3.EUR.MAF_001.chr \
    	--targ_pred F \
    	--save_ref_expr T \
    	--save_score F \
    	--plink plink \
    	--ref_maf {params.resdir}/data/1kg/1KG.Phase3.EUR.MAF_001.chr \
    	--pigz {params.resdir}/software/pigz/pigz/pigz \
    	--memory 40000 \
      --n_cores 5 \
    	--output {params.resdir}/data/predicted_expression/{wildcards.weight} > {log} 2>&1"

# Format expression data for TWAS-GSEA (i.e. remove PANEL from column names)
rule format_pred:
  input:
    f"{resdir}/data/predicted_expression/{{weight}}/Reference_Expression/Reference_Expression_{{weight}}.txt.gz"
  output:
    touch(f"{resdir}/data/predicted_expression/format_pred_{{weight}}.done")
  benchmark:
    f"{resdir}/benchmarks/format_pred_{{weight}}.tsv"
  conda:
    "../envs/main.yaml"
  params:
    resdir=resdir
  log:
    f"{resdir}/logs/format_pred-{{weight}}.log"
  shell:
    "(zcat {params.resdir}/data/predicted_expression/{wildcards.weight}/Reference_Expression/Reference_Expression_{wildcards.weight}.txt.gz | sed -e s/{wildcards.weight}.//g | gzip > {params.resdir}/data/predicted_expression/{wildcards.weight}/Reference_Expression/Reference_Expression_{wildcards.weight}_mod.txt.gz; mv {params.resdir}/data/predicted_expression/{wildcards.weight}/Reference_Expression/Reference_Expression_{wildcards.weight}_mod.txt.gz {params.resdir}/data/predicted_expression/{wildcards.weight}/Reference_Expression/Reference_Expression_{wildcards.weight}.txt.gz) > {log} 2>&1"

####
# Install lme4qtl
####
# Note the version of conda was not working in R 4.0.2

rule install_lme4qtl:
  input:
    f"{workflow.basedir}/envs/main.yaml"
  output:
    touch(f"{resdir}/software/install_lme4qtl.done")
  benchmark:
    f"{resdir}/benchmarks/install_lme4qtl.tsv"
  conda:
    "../envs/main.yaml"
  log:
    f"{resdir}/logs/install_lme4qtl.log"
  shell:
    "Rscript --vanilla -e 'devtools::install_github(\"variani/lme4qtl@0c173ea8d8386b205f62ad642698519a861650b4\", upgrade = \"never\")' > {log} 2>&1"

####
# Format ROSMAP and Banner PWAS data
####

rule format_pwas_data:
  output:
    f"{resdir}/data/banner_twas/Banner.n152.fusion.WEIGHTS/train_weights_withN.pos"
  benchmark:
    f"{resdir}/benchmarks/format_pwas_data.tsv"
  conda:
    "../envs/main.yaml"
  params:
    rosmap_fusion= config["rosmap_fusion"],
    banner_fusion= config["banner_fusion"]
  log:
    f"{resdir}/logs/format_pwas_data.log"
  shell:
    "Rscript --vanilla {workflow.basedir}/scripts/format_pwas_data.R --pipeline_dir {workflow.basedir} \
      --rosmap {params.rosmap_fusion} \
      --banner {params.banner_fusion} > {log} 2>&1"

# Format drugtargetor database for TWAS-GSEA
rule format_drug_targetor_for_twas_gsea:
  input:
    rules.download_drug_targetor.output,
    rules.download_magma_gene_loc.output
  output:
    f"{resdir}/data/drug_targetor/wholedatabase_for_targetor_directional.prop"
  benchmark:
    f"{resdir}/benchmarks/format_drug_targetor_for_twas_gsea.tsv"
  conda:
    "../envs/main.yaml"
  log:
    f"{resdir}/logs/format_drug_targetor_for_twas_gsea.log"
  shell:
    "Rscript --vanilla {workflow.basedir}/scripts/format_drug_targetor_for_twas_gsea.R --pipeline_dir {workflow.basedir} > {log} 2>&1"

####
# Download SMR
####

rule download_smr:
  output:
    f"{resdir}/software/smr/smr_linux_x86_64"
  benchmark:
    f"{resdir}/benchmarks/download_smr.tsv"
  conda:
    "../envs/main.yaml"
  params:
    resdir=resdir
  log:
    f"{resdir}/logs/download_smr.log"
  shell:
    "(rm -rf {params.resdir}/software/smr; \
    mkdir -p {params.resdir}/software/smr; \
    wget -O {params.resdir}/software/smr/smr_Linux.zip https://yanglab.westlake.edu.cn/software/smr/download/smr_Linux.zip; \
    unzip {params.resdir}/software/smr/smr_Linux.zip -d {params.resdir}/software/smr; \
    rm {params.resdir}/software/smr/smr_Linux.zip) > {log} 2>&1"

####
# Format ROSMAP SMR data
####

rule format_rosmap_smr_data:
  input:
    rules.download_smr.output,
    rules.prep_1kg.output
  output:
    f"{resdir}/data/rosmap_smr/ROSMAP.n376.pQTL.MatrixQTL.txt.besd.epi"
  benchmark:
    f"{resdir}/benchmarks/format_rosmap_smr_data.tsv"
  conda:
    "../envs/main.yaml"
  params:
    rosmap_smr= config["rosmap_smr"],
  log:
    f"{resdir}/logs/format_rosmap_smr_data.log"
  shell:
    "Rscript --vanilla {workflow.basedir}/scripts/format_rosmap_smr_data.R --pipeline_dir {workflow.basedir} \
      --rosmap {params.rosmap_smr} > {log} 2>&1"

####
# Download PsychENCODE data for SMR
####

rule download_psychencode_smr:
  output:
    directory(f"{resdir}/data/psychencode_smr/PsychENCODE_cis_eqtl_HCP100_summary/")
  benchmark:
    f"{resdir}/benchmarks/download_psychencode_smr.tsv"
  conda:
    "../envs/main.yaml"
  params:
    resdir=resdir
  log:
    f"{resdir}/logs/download_psychencode_smr.log"
  shell:
    "(rm -rf {params.resdir}/data/psychencode_smr; \
    mkdir -p {params.resdir}/data/psychencode_smr; \
    wget --no-check-certificate -O {params.resdir}/data/psychencode_smr/PsychENCODE_cis_eqtl_HCP100_summary.tar.gz http://cnsgenomics.com/data/SMR/PsychENCODE_cis_eqtl_HCP100_summary.tar.gz; \
    tar -xvzf {params.resdir}/data/psychencode_smr/PsychENCODE_cis_eqtl_HCP100_summary.tar.gz -C {params.resdir}/data/psychencode_smr; \
    rm {params.resdir}/data/psychencode_smr/PsychENCODE_cis_eqtl_HCP100_summary.tar.gz) > {log} 2>&1"

##
# Download MetaBrain data in SMR format
##

# Basalganglia
rule download_MetaBrain_Basalganglia:
  output:
    directory(f"{resdir}/data/MetaBrain/Basalganglia")
  benchmark:
    f"{resdir}/benchmarks/download_MetaBrain_Basalganglia.tsv"
  conda:
    "../envs/main.yaml"
  params:
    resdir=resdir
  log:
    f"{resdir}/logs/download_MetaBrain_Basalganglia.log"
  shell:
    "(mkdir -p {params.resdir}/data/MetaBrain/Basalganglia; \
    wget -O {params.resdir}/data/MetaBrain/Basalganglia/2020-05-26-Basalganglia-EUR-smr.zip https://download.metabrain.nl/2020-05-26-release/2020-05-26-CisEQTLSummaryStats/2020-05-26-Basalganglia-EUR/2020-05-26-Basalganglia-EUR-smr.zip; \
    unzip -d {params.resdir}/data/MetaBrain/Basalganglia/ {params.resdir}/data/MetaBrain/Basalganglia/2020-05-26-Basalganglia-EUR-smr.zip; \
    rm {params.resdir}/data/MetaBrain/Basalganglia/2020-05-26-Basalganglia-EUR-smr.zip) > {log} 2>&1"

# Cerebellum
rule download_MetaBrain_Cerebellum:
  output:
    directory(f"{resdir}/data/MetaBrain/Cerebellum")
  benchmark:
    f"{resdir}/benchmarks/download_MetaBrain_Cerebellum.tsv"
  conda:
    "../envs/main.yaml"
  params:
    resdir=resdir
  log:
    f"{resdir}/logs/download_MetaBrain_Cerebellum.log"
  shell:
    "(mkdir -p {params.resdir}/data/MetaBrain/Cerebellum; \
    wget -O {params.resdir}/data/MetaBrain/Cerebellum/2020-05-26-Cerebellum-EUR-smr.zip https://download.metabrain.nl/2020-05-26-release/2020-05-26-CisEQTLSummaryStats/2020-05-26-Cerebellum-EUR/2020-05-26-Cerebellum-EUR-smr.zip; \
    unzip -d {params.resdir}/data/MetaBrain/Cerebellum/ {params.resdir}/data/MetaBrain/Cerebellum/2020-05-26-Cerebellum-EUR-smr.zip; \
    rm {params.resdir}/data/MetaBrain/Cerebellum/2020-05-26-Cerebellum-EUR-smr.zip) > {log} 2>&1"

# Cortex
rule download_MetaBrain_Cortex:
  output:
    directory(f"{resdir}/data/MetaBrain/Cortex")
  benchmark:
    f"{resdir}/benchmarks/download_MetaBrain_Cortex.tsv"
  conda:
    "../envs/main.yaml"
  params:
    resdir=resdir
  log:
    f"{resdir}/logs/download_MetaBrain_Cortex.log"
  shell:
    "(mkdir -p {params.resdir}/data/MetaBrain/Cortex; \
    wget -O {params.resdir}/data/MetaBrain/Cortex/2020-05-26-Cortex-EUR-smr.zip https://download.metabrain.nl/2020-05-26-release/2020-05-26-CisEQTLSummaryStats/2020-05-26-Cortex-EUR/2020-05-26-Cortex-EUR-smr.zip; \
    unzip -d {params.resdir}/data/MetaBrain/Cortex/ {params.resdir}/data/MetaBrain/Cortex/2020-05-26-Cortex-EUR-smr.zip; \
    rm {params.resdir}/data/MetaBrain/Cortex/2020-05-26-Cortex-EUR-smr.zip) > {log} 2>&1"

# Hippocampus
rule download_MetaBrain_Hippocampus:
  output:
    directory(f"{resdir}/data/MetaBrain/Hippocampus")
  benchmark:
    f"{resdir}/benchmarks/download_MetaBrain_Hippocampus.tsv"
  conda:
    "../envs/main.yaml"
  params:
    resdir=resdir
  log:
    f"{resdir}/logs/download_MetaBrain_Hippocampus.log"
  shell:
    "(mkdir -p {params.resdir}/data/MetaBrain/Hippocampus; \
    wget -O {params.resdir}/data/MetaBrain/Hippocampus/2020-05-26-Hippocampus-EUR-smr.zip https://download.metabrain.nl/2020-05-26-release/2020-05-26-CisEQTLSummaryStats/2020-05-26-Hippocampus-EUR/2020-05-26-Hippocampus-EUR-smr.zip; \
    unzip -d {params.resdir}/data/MetaBrain/Hippocampus/ {params.resdir}/data/MetaBrain/Hippocampus/2020-05-26-Hippocampus-EUR-smr.zip; \
    rm {params.resdir}/data/MetaBrain/Hippocampus/2020-05-26-Hippocampus-EUR-smr.zip) > {log} 2>&1"

# Spinalcord
rule download_MetaBrain_Spinalcord:
  output:
    directory(f"{resdir}/data/MetaBrain/Spinalcord")
  benchmark:
    f"{resdir}/benchmarks/download_MetaBrain_Spinalcord.tsv"
  conda:
    "../envs/main.yaml"
  params:
    resdir=resdir
  log:
    f"{resdir}/logs/download_MetaBrain_Spinalcord.log"
  shell:
    "(mkdir -p {params.resdir}/data/MetaBrain/Spinalcord; \
    wget -O {params.resdir}/data/MetaBrain/Spinalcord/2020-05-26-Spinalcord-EUR-smr.zip https://download.metabrain.nl/2020-05-26-release/2020-05-26-CisEQTLSummaryStats/2020-05-26-Spinalcord-EUR/2020-05-26-Spinalcord-EUR-smr.zip; \
    unzip -d {params.resdir}/data/MetaBrain/Spinalcord/ {params.resdir}/data/MetaBrain/Spinalcord/2020-05-26-Spinalcord-EUR-smr.zip; \
    rm {params.resdir}/data/MetaBrain/Spinalcord/2020-05-26-Spinalcord-EUR-smr.zip) > {log} 2>&1"

rule download_MetaBrain_all:
  input:
    rules.download_MetaBrain_Basalganglia.output,
    rules.download_MetaBrain_Cerebellum.output,
    rules.download_MetaBrain_Cortex.output,
    rules.download_MetaBrain_Hippocampus.output,
    rules.download_MetaBrain_Spinalcord.output
  output:
    touch(f'{resdir}/data/MetaBrain_download.out')

# Update variant IDs in MetaBrain SMR files
rule format_metabrain_esi:
  input:
    f"{resdir}/data/MetaBrain_download.out"
  output:
    touch(f"{resdir}/data/MetaBrain/format_MetaBrain_esi.out")
  benchmark:
    f"{resdir}/benchmarks/format_metabrain_esi.tsv"
  conda:
    "../envs/main.yaml"
  log:
    f"{resdir}/logs/format_metabrain_esi.log"
  shell:
    "Rscript --vanilla {workflow.basedir}/scripts/format_metabrain_esi.R --pipeline_dir {workflow.basedir} > {log} 2>&1"

# Download eQTLGen data in SMR format
rule download_eqtlgen:
  output:
    touch(f"{resdir}/data/eqtlgen.done")
  benchmark:
    f"{resdir}/benchmarks/download_eqtlgen.tsv"
  conda:
    "../envs/main.yaml"
  params:
    resdir=resdir
  log:
    f"{resdir}/logs/download_eqtlgen.log"
  shell:
    "(rm -rf {params.resdir}/data/eqtlgen; \
    mkdir {params.resdir}/data/eqtlgen; \
    wget -O {params.resdir}/data/eqtlgen/cis-eQTL-SMR_20191212.tar.gz https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/SMR_formatted/cis-eQTL-SMR_20191212.tar.gz; \
    tar -xvzf {params.resdir}/data/eqtlgen/cis-eQTL-SMR_20191212.tar.gz -C {params.resdir}/data/eqtlgen/; \
    rm {params.resdir}/data/eqtlgen/cis-eQTL-SMR_20191212.tar.gz; \
    gunzip {params.resdir}/data/eqtlgen/*) > {log} 2>&1"

####
# Download GCSC
####

rule install_gcsc:
  output:
    directory(f"{resdir}/software/GCSC/")
  benchmark:
    f"{resdir}/benchmarks/install_gcsc.tsv"
  conda:
    "../envs/main.yaml"
  log:
    f"{resdir}/logs/install_gcsc.log"
  shell:
    "(git clone https://github.com/ksiewert/GCSC.git {output}; \
     cd {output}; \
     git reset --hard b10ea77b9a43399801b46ef70c80516599264123) > {log} 2>&1"

####
# Download GCSC gene co-regulation scores
####

gcsc_tissues=config["gcsc_tissues"]

rule download_gcsc_coreg:
  output:
    f"{resdir}/data/GCSC/coreg/{{gcsc_tissue}}_geneNames.txt"
  benchmark:
    f"{resdir}/benchmarks/download_gcsc_coreg_{{gcsc_tissue}}.tsv"
  conda:
    "../envs/main.yaml"
  params:
    resdir=resdir
  log:
    f"{resdir}/logs/download_gcsc_coreg-{{gcsc_tissue}}.log"
  shell:
    "(wget -O {params.resdir}/data/GCSC/coreg/{wildcards.gcsc_tissue}_coregscores.npz https://storage.googleapis.com/broad-alkesgroup-public/GCSC/Coreg_scores/{wildcards.gcsc_tissue}_coregscores.npz; \
    wget -O {params.resdir}/data/GCSC/coreg/{wildcards.gcsc_tissue}_geneNames.txt https://storage.googleapis.com/broad-alkesgroup-public/GCSC/Coreg_scores/{wildcards.gcsc_tissue}_geneNames.txt) > {log} 2>&1"

rule download_gcsc_coreg_all_tissue:
    input: expand(f"{resdir}/data/GCSC/coreg/{{gcsc_tissue}}_geneNames.txt", gcsc_tissue=gcsc_tissues)

####
# Download corresponding GTEx v7 TWAS weights
####

rule download_gcsc_twas_weights:
  output:
    directory(f"{resdir}/data/GCSC/twas_weights/GTEx.{{gcsc_tissue}}.P01")
  benchmark:
    f"{resdir}/benchmarks/download_gcsc_twas_weights_{{gcsc_tissue}}.tsv"
  conda:
    "../envs/main.yaml"
  params:
    resdir=resdir
  log:
    f"{resdir}/logs/download_gcsc_twas_weights-{{gcsc_tissue}}.log"
  shell:
    "(mkdir {params.resdir}/data/GCSC/twas_weights/GTEx.{wildcards.gcsc_tissue}.P01; \
    wget -O {params.resdir}/data/GCSC/twas_weights/GTEx.{wildcards.gcsc_tissue}.P01.tar.bz2 http://gusevlab.org/projects/fusion/weights/GTEx.{wildcards.gcsc_tissue}.P01.tar.bz2; \
    tar xjvf {params.resdir}/data/GCSC/twas_weights/GTEx.{wildcards.gcsc_tissue}.P01.tar.bz2 -C {params.resdir}/data/GCSC/twas_weights/GTEx.{wildcards.gcsc_tissue}.P01; \
    rm {params.resdir}/data/GCSC/twas_weights/GTEx.{wildcards.gcsc_tissue}.P01.tar.bz2) > {log} 2>&1"

rule download_gcsc_twas_weights_all_tissue:
    input: expand(f"{resdir}/data/GCSC/twas_weights/GTEx.{{gcsc_tissue}}.P01", gcsc_tissue=gcsc_tissues)

####
# Pre-compute the gene-gene predicted-expression correlation matrix once per
# weight panel. The result is reused across every (gwas, gene-set) call to
# TWAS-GSEA-fast.R for that panel, so this rule has no {gwas} wildcard.
####

rule build_twas_gsea_cormat:
  resources:
    mem_mb=50000,
    cpus=5
  input:
    rules.install_twas_gsea.output,
    f"{resdir}/data/predicted_expression/format_pred_{{weight}}.done"
  output:
    f"{resdir}/data/predicted_expression/{{weight}}/Reference_Expression/{{weight}}.CorMat.RDS"
  benchmark:
    f"{resdir}/benchmarks/build_twas_gsea_cormat_{{weight}}.tsv"
  conda:
    "../envs/main.yaml"
  params:
    resdir=resdir
  log:
    f"{resdir}/logs/build_twas_gsea_cormat-{{weight}}.log"
  shell:
    "Rscript --vanilla {params.resdir}/software/TWAS-GSEA/build_cor_matrix.R \
      --expression_ref {params.resdir}/data/predicted_expression/{wildcards.weight}/Reference_Expression/Reference_Expression_{wildcards.weight}.txt.gz \
      --pos {params.resdir}/data/fusion_snp_weights/{wildcards.weight}/{wildcards.weight}.pos \
      --min_r2 0.01 \
      --n_cores 5 \
      --output {params.resdir}/data/predicted_expression/{wildcards.weight}/Reference_Expression/{wildcards.weight} > {log} 2>&1"


####
# Prepare all resources
####

resource_inputs = list()

# Core infrastructure (always needed)
resource_inputs.extend([
    rules.install_liftover.output,
    rules.download_liftover_track.output,
    rules.prep_1kg.output,
    rules.install_genoutils.output,
    rules.download_biomart.output
])

# LDSC
if config["ldsc"] == "T":
    resource_inputs.extend([
        rules.install_ldsc.output,
        rules.download_ldsc_scores.output,
        rules.download_ldsc_hm3.output
    ])

# COJO
if config["cojo"] == "T":
    resource_inputs.append(rules.download_gcta.output)

# MAGMA
if config["magma_gene"] == "T" or config["magma_drugtargetor"] == "T" or config["tissue_magma"] == "T":
    resource_inputs.append(rules.magma_annot.output)

if config["magma_drugtargetor"] == "T":
    resource_inputs.extend([
        rules.download_atc.output,
        rules.format_drug_targetor.output
    ])

if config["tissue_magma"] == "T":
    resource_inputs.append(rules.prep_tissue_exp.output)

# TWAS
if config["twas_panel_fusion"] == "T" or config["twas_panel_psychencode"] == "T":
    resource_inputs.extend([
        rules.install_fusion.output,
        rules.install_plink2R.output
    ])

if config["twas_panel_psychencode"] == "T":
    resource_inputs.append(rules.format_psychencode.output)

if config["twas_panel_fusion"] == "T":
    resource_inputs.extend([
        rules.update_gtex_coord_all_panel.input,
        rules.insert_n_nongtex_all_panel.input
    ])

if config["twas_conditional"] == "T":
    resource_inputs.append(rules.download_glist.output)

# TWAS-GSEA
if config["twas_gsea_drugtargetor"] == "T":
    resource_inputs.extend([
        rules.install_twas_gsea.output,
        rules.install_feature_pred.output,
        rules.install_pigz.output,
        rules.install_lme4qtl.output,
        rules.format_drug_targetor_for_twas_gsea.output,
        expand(f"{resdir}/data/predicted_expression/format_pred_{{weight}}.done", weight=weights_nosplice),
        expand(f"{resdir}/data/predicted_expression/{{weight}}/Reference_Expression/{{weight}}.CorMat.RDS", weight=weights_nosplice),
    ])

# PWAS
if config["pwas_panel_rosmap"] == "T" or config["pwas_panel_banner"] == "T":
    resource_inputs.append(rules.format_pwas_data.output)

# SMR
any_smr = any(config.get(k) == "T" for k in [
    "smr_expression_panel_psychencode",
    "smr_expression_panel_metabrain_basalganglia",
    "smr_expression_panel_metabrain_cerebellum",
    "smr_expression_panel_metabrain_cortex",
    "smr_expression_panel_metabrain_hippocampus",
    "smr_expression_panel_metabrain_spinalcord",
    "smr_expression_panel_eqtlgen",
    "smr_protein_panel_rosmap"
])

if any_smr:
    resource_inputs.append(rules.download_smr.output)

if config["smr_expression_panel_psychencode"] == "T":
    resource_inputs.append(rules.download_psychencode_smr.output)

any_metabrain = any(config.get(k) == "T" for k in [
    "smr_expression_panel_metabrain_basalganglia",
    "smr_expression_panel_metabrain_cerebellum",
    "smr_expression_panel_metabrain_cortex",
    "smr_expression_panel_metabrain_hippocampus",
    "smr_expression_panel_metabrain_spinalcord"
])

if any_metabrain:
    resource_inputs.append(rules.format_metabrain_esi.output)

if config["smr_expression_panel_eqtlgen"] == "T":
    resource_inputs.append(rules.download_eqtlgen.output)

if config["smr_protein_panel_rosmap"] == "T":
    resource_inputs.append(rules.format_rosmap_smr_data.output)

# GCSC
if config["gcsc"] == "T":
    resource_inputs.extend([
        rules.install_gcsc.output,
        rules.download_gcsc_coreg_all_tissue.input,
        rules.download_gcsc_twas_weights_all_tissue.input
    ])

rule prepare_resources:
  input:
    resource_inputs

