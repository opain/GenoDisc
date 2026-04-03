####
# Download MAGMA
####

rule download_magma:
  output:
    f"{resdir}/software/magma/magma"
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
    f"{resdir}/data/drug_targetor/wholedatabase_for_targetor.gmt"
  conda:
    "../envs/main.yaml"
  log:
    f"{resdir}/logs/format_drug_targetor.log"
  shell:
    "Rscript scripts/format_drug_targetor.R > {log} 2>&1"

####
# Download and format GTEx TPM data
####

rule prep_tissue_exp:
  input:
    rules.download_magma_gene_loc.output,
    "scripts/prep_tissue_exp.R"
  output:
    f"{resdir}/data/gtex/GTEx_v8_group.tsv"
  conda:
    "../envs/main.yaml"
  log:
    f"{resdir}/logs/prep_tissue_exp.log"
  shell:
    "Rscript scripts/prep_tissue_exp.R > {log} 2>&1"

##########
# Analyse GWAS summary statistics
##########

####
# Run MAGMA
####

# Run gene level association analysis
rule magma_gene_level:
  input:
    "{outdir}/data/gwas_sumstat/{gwas}/{gwas}.cleaned.gz",
    rules.magma_annot.output
  output:
    "{outdir}/results/{gwas}/magma/magma_gene_level.genes.raw"
  conda:
    "../envs/main.yaml"
  params:
    resdir=resdir
  log:
    "{outdir}/logs/magma_gene_level-{gwas}.log"
  shell:
    "(gzip -f -d -c {outdir}/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}.cleaned.gz > {outdir}/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}.cleaned; \
    {params.resdir}/software/magma/magma \
      --bfile {params.resdir}/data/magma_ref/g1000_eur \
      --pval {outdir}/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}.cleaned use=SNP,P ncol=N \
      --gene-annot {params.resdir}/data/magma/NCBI37.3.genes.annot \
      --out {outdir}/results/{wildcards.gwas}/magma/magma_gene_level; rm {outdir}/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}.cleaned) > {log} 2>&1"

# Format the MAGMA gene results 
rule format_magma_gene_results:
  input:
    "{outdir}/results/{gwas}/magma/magma_gene_level.genes.raw"
  output:
    "{outdir}/results/{gwas}/magma/magma_gene_level.clean.csv"
  conda:
    "../envs/main.yaml"
  params:
    config_file=config['config_file']
  log:
    "{outdir}/logs/format_magma_gene_results-{gwas}.log"
  shell:
    "Rscript scripts/format_magma_gene_results.R \
      --gwas {wildcards.gwas} \
      --config_file {params.config_file} > {log} 2>&1"

# Run Drug Targetor enrichment analysis
rule magma_drug_targetor:
  input:
    "{outdir}/results/{gwas}/magma/magma_gene_level.genes.raw",
    f"{resdir}/data/drug_targetor/wholedatabase_for_targetor.gmt"
  output:
    "{outdir}/results/{gwas}/magma/magma_drug_targetor.gsa.out"
  conda:
    "../envs/main.yaml"
  params:
    resdir=resdir
  log:
    "{outdir}/logs/magma_drug_targetor-{gwas}.log"
  shell:
    "{params.resdir}/software/magma/magma \
      --gene-results {outdir}/results/{wildcards.gwas}/magma/magma_gene_level.genes.raw \
      --set-annot {params.resdir}/data/drug_targetor/wholedatabase_for_targetor.gmt \
      --out {outdir}/results/{wildcards.gwas}/magma/magma_drug_targetor > {log} 2>&1"

# Format the MAGMA GSEA results 
rule format_magma_results:
  input:
    "{outdir}/results/{gwas}/magma/magma_drug_targetor.gsa.out",
    rules.download_atc.output
  output:
    "{outdir}/results/{gwas}/magma/magma_drug_targetor_atc_res.csv"
  conda:
    "../envs/main.yaml"
  params:
    config_file=config['config_file']
  log:
    "{outdir}/logs/format_magma_results-{gwas}.log"
  shell:
    "Rscript scripts/format_magma_gsea_results.R \
      --gwas {wildcards.gwas} \
      --config_file {params.config_file} > {log} 2>&1"

# Compare TWAS signiture compared to enriched drugs in MAGMA
rule comp_magma_gsea_twas_results:
  input:
    "{outdir}/results/{gwas}/magma/magma_drug_targetor_atc_res.csv",
    "{outdir}/results/{gwas}/twas/{gwas}_twas_GW_clean.txt.gz"
  output:
    "{outdir}/results/{gwas}/magma/magma_drug_targetor_twas_comp.csv"
  conda:
    "../envs/main.yaml"
  params:
    config_file=config['config_file']
  log:
    "{outdir}/logs/comp_magma_gsea_twas_results-{gwas}.log"
  shell:
    "Rscript scripts/comp_magma_gsea_twas_results.R \
      --gwas {wildcards.gwas} \
      --config_file {params.config_file} > {log} 2>&1"

# Perform tissue specific enrichment analysis
rule magma_tissue_spec:
  input:
    "{outdir}/results/{gwas}/magma/magma_gene_level.genes.raw",
    f"{resdir}/data/gtex/GTEx_v8_group.tsv"
  output:
    "{outdir}/results/{gwas}/magma/magma_tissue_spec.gsa.out"
  conda:
    "../envs/main.yaml"
  params:
    resdir=resdir
  log:
    "{outdir}/logs/magma_tissue_spec-{gwas}.log"
  shell:
    "{params.resdir}/software/magma/magma \
      --gene-results {outdir}/results/{wildcards.gwas}/magma/magma_gene_level.genes.raw \
      --gene-covar {params.resdir}/data/gtex/GTEx_v8_tissue.tsv \
      --model direction-covar=greater condition-hide=Average \
      --out {outdir}/results/{wildcards.gwas}/magma/magma_tissue_spec > {log} 2>&1"

# Perform tissue group enrichment analysis
rule magma_tissue_group:
  input:
    "{outdir}/results/{gwas}/magma/magma_gene_level.genes.raw",
    f"{resdir}/data/gtex/GTEx_v8_group.tsv"
  output:
    "{outdir}/results/{gwas}/magma/magma_tissue_group.gsa.out"
  conda:
    "../envs/main.yaml"
  params:
    resdir=resdir
  log:
    "{outdir}/logs/magma_tissue_group-{gwas}.log"
  shell:
    "{params.resdir}/software/magma/magma \
      --gene-results {outdir}/results/{wildcards.gwas}/magma/magma_gene_level.genes.raw \
      --gene-covar {params.resdir}/data/gtex/GTEx_v8_group.tsv \
      --model direction-covar=greater condition-hide=Average \
      --out {outdir}/results/{wildcards.gwas}/magma/magma_tissue_group > {log} 2>&1"

# Perform conditional analysis of tissues
rule magma_tissue_conditional:
  input:
    "{outdir}/results/{gwas}/magma/magma_tissue_spec.gsa.out",
    "scripts/magma_tissue_conditional.R"
  output:
    touch("{outdir}/results/{gwas}/magma/magma_property_conditional.done")
  conda:
    "../envs/main.yaml"
  params:
    config_file= config['config_file']
  log:
    "{outdir}/logs/magma_tissue_conditional-{gwas}.log"
  shell:
    "Rscript scripts/magma_tissue_conditional.R \
      --config_file {params.config_file} \
      --gwas {wildcards.gwas} > {log} 2>&1"



