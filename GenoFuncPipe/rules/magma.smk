####
# Download MAGMA
####

rule download_magma:
  output:
    "resources/software/magma/magma"
  conda:
    "../envs/GenoFunc.yaml"
  shell:
    "wget -O resources/software/magma.zip https://ctg.cncr.nl/software/MAGMA/prog/magma_v1.10.zip; unzip resources/software/magma.zip -d resources/software/magma; rm resources/software/magma.zip"

####
# Download MAGMA gene locations
####

rule download_magma_gene_loc:
  output:
    "resources/data/magma/NCBI37.3.gene.loc"
  conda:
    "../envs/GenoFunc.yaml"
  shell:
    "wget -O resources/data/magma.zip https://ctg.cncr.nl/software/MAGMA/aux_files/NCBI37.3.zip; unzip resources/data/magma.zip -d resources/data/magma; rm resources/data/magma.zip"

####
# Download MAGMA reference
####

rule download_magma_ref:
  output:
    "resources/data/magma_ref/g1000_eur.bed"
  conda:
    "../envs/GenoFunc.yaml"
  shell:
    "wget -O resources/data/magma.zip https://ctg.cncr.nl/software/MAGMA/ref_data/g1000_eur.zip; unzip resources/data/magma.zip -d resources/data/magma_ref; rm resources/data/magma.zip"

####
# Create MAGMA annotation file
####

rule magma_annot:
  input:
    rules.download_magma.output,
    rules.download_magma_gene_loc.output,
    rules.download_magma_ref.output
  output:
    "resources/data/magma/NCBI37.3.genes.annot"
  conda: 
    "../envs/GenoFunc.yaml"
  shell:
    "resources/software/magma/magma \
      --annotate window=35,10 \
    	--snp-loc resources/data/magma_ref/g1000_eur.bim \
    	--gene-loc resources/data/magma/NCBI37.3.gene.loc \
    	--out resources/data/magma/NCBI37.3"

####
# Download ATC codes
####

rule download_atc:
  output:
    "resources/data/atc/atc_20220201.txt"
  conda:
    "../envs/GenoFunc.yaml"
  shell:
    "wget -O resources/data/2022-02-01-v3extracts.zip https://www.pbs.gov.au/downloads/2022/02/2022-02-01-v3extracts.zip; mkdir -p resources/data/atc; unzip resources/data/2022-02-01-v3extracts.zip -d resources/data/atc; rm resources/data/2022-02-01-v3extracts.zip"

####
# Download and format DrugTargetor database
####

rule download_drug_targetor:
  output:
    "resources/data/drug_targetor/wholedatabase_for_targetor"
  conda:
    "../envs/GenoFunc.yaml"
  shell:
    "mkdir -p resources/data/drug_targetor/ ; wget -O resources/data/drug_targetor/wholedatabase_for_targetor https://github.com/hagax8/drugtargetor/raw/master/wholedatabase_for_targetor"

rule format_drug_targetor:
  input:
    rules.download_drug_targetor.output,
    rules.download_magma_gene_loc.output
  output:
    "resources/data/drug_targetor/wholedatabase_for_targetor.gmt"
  conda:
    "../envs/GenoFunc.yaml"
  shell:
    "Rscript scripts/format_drug_targetor.R"

##########
# Analyse GWAS summary statistics
##########

####
# Run MAGMA
####

# Run gene level association analysis
rule magma_gene_level:
  input:
    "resources/data/gwas_sumstat/{gwas}/{gwas}.cleaned.gz",
    rules.magma_annot.output
  output:
    "results/{gwas}/magma/magma_gene_level.genes.raw"
  conda: 
    "../envs/GenoFunc.yaml"
  shell:
    "gzip -f -d -c resources/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}.cleaned.gz > resources/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}.cleaned; resources/software/magma/magma \
      --bfile resources/data/magma_ref/g1000_eur \
      --pval resources/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}.cleaned use=SNP,P ncol=N \
      --gene-annot resources/data/magma/NCBI37.3.genes.annot \
      --out results/{wildcards.gwas}/magma/magma_gene_level; rm resources/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}.cleaned"

# Format the MAGMA gene results 
rule format_magma_gene_results:
  input:
    "results/{gwas}/magma/magma_gene_level.genes.raw"
  output:
    "results/{gwas}/magma/magma_gene_level.clean.csv"
  conda: 
    "../envs/GenoFunc.yaml"
  params:
    lincs_siginfo= config["lincs_siginfo"]
  shell:
    "Rscript scripts/format_magma_gene_results.R \
    --gwas {wildcards.gwas}"

# Run Drug Targetor enrichment analysis
rule magma_drug_targetor:
  input:
    "results/{gwas}/magma/magma_gene_level.genes.raw",
    "resources/data/drug_targetor/wholedatabase_for_targetor.gmt"
  output:
    "results/{gwas}/magma/magma_drug_targetor.gsa.out"
  conda: 
    "../envs/GenoFunc.yaml"
  shell:
    "resources/software/magma/magma \
      --gene-results results/{wildcards.gwas}/magma/magma_gene_level.genes.raw \
      --set-annot resources/data/drug_targetor/wholedatabase_for_targetor.gmt \
      --out results/{wildcards.gwas}/magma/magma_drug_targetor"

# Format the MAGMA GSEA results 
rule format_magma_results:
  input:
    "results/{gwas}/magma/magma_drug_targetor.gsa.out",
    rules.download_atc.output
  output:
    "results/{gwas}/magma/magma_drug_targetor_atc_res.csv"
  conda: 
    "../envs/GenoFunc.yaml"
  params:
    lincs_siginfo= config["lincs_siginfo"]
  shell:
    "Rscript scripts/format_magma_gsea_results.R \
    --gwas {wildcards.gwas}"

