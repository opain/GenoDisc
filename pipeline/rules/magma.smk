##########
# Analyse GWAS summary statistics
##########

####
# Run MAGMA
####

# Run gene level association analysis
rule magma_gene_level:
  input:
    "{outdir}/results/{gwas}/gwas_sumstat/{gwas}.cleaned.gz",
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
    "(gzip -f -d -c {outdir}/results/{wildcards.gwas}/gwas_sumstat/{wildcards.gwas}.cleaned.gz > {outdir}/results/{wildcards.gwas}/gwas_sumstat/{wildcards.gwas}.cleaned; \
    {params.resdir}/software/magma/magma \
      --bfile {params.resdir}/data/magma_ref/g1000_eur \
      --pval {outdir}/results/{wildcards.gwas}/gwas_sumstat/{wildcards.gwas}.cleaned use=SNP,P ncol=N \
      --gene-annot {params.resdir}/data/magma/NCBI37.3.genes.annot \
      --out {outdir}/results/{wildcards.gwas}/magma/magma_gene_level; rm {outdir}/results/{wildcards.gwas}/gwas_sumstat/{wildcards.gwas}.cleaned) > {log} 2>&1"

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
