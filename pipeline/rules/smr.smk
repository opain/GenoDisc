##########
# Analyse GWAS summary statistics
##########

####
# Format sumstats for SMR
####

rule format_sumstats_smr:
  input:
    "{outdir}/results/{gwas}/gwas_sumstat/{gwas}.cleaned.gz"
  output:
    "{outdir}/results/{gwas}/gwas_sumstat/{gwas}.cleaned.cojo"
  conda:
    "../envs/main.yaml"
  params:
    config_file=config['config_file']
  log:
    "{outdir}/logs/format_sumstats_smr-{gwas}.log"
  shell:
    "Rscript scripts/format_sumstats_smr.R \
      --gwas {wildcards.gwas} \
      --config_file {params.config_file} > {log} 2>&1"

####
# Run SMR using eQTL
####

rule run_psychencode_smr:
  input:
    "{outdir}/results/{gwas}/gwas_sumstat/{gwas}.cleaned.cojo",
    rules.download_psychencode_smr.output
  output:
    "{outdir}/results/{gwas}/smr/psychencode/{gwas}_smr_psychencode_chr{chr}.smr"
  conda:
    "../envs/main.yaml"
  params:
    resdir=resdir
  log:
    "{outdir}/logs/run_psychencode_smr-{gwas}-chr{chr}.log"
  shell:
    "{params.resdir}/software/smr/smr_linux_x86_64 \
    --bfile {params.resdir}/data/1kg/1KG.Phase3.EUR.MAF_001.chr{wildcards.chr} \
    --gwas-summary {outdir}/results/{wildcards.gwas}/gwas_sumstat/{wildcards.gwas}.cleaned.cojo \
    --beqtl-summary {params.resdir}/data/psychencode_smr/PsychENCODE_cis_eqtl_HCP100_summary/Gandal_PsychENCODE_eQTL_HCP100+gPCs20_QTLtools \
    --out {outdir}/results/{wildcards.gwas}/smr/psychencode/{wildcards.gwas}_smr_psychencode_chr{wildcards.chr} > {log} 2>&1"

rule run_psychencode_smr_chr:
    input:
      lambda w: expand("{outdir}/results/{gwas}/smr/psychencode/{gwas}_smr_psychencode_chr{chr}.smr", gwas=w.gwas, chr=chromosomes, outdir={outdir})
    output:
      touch("{outdir}/results/{gwas}/checks/psychencode_smr_all_chr.done")

########
# Run SMR with eQTLGen
########

rule run_eqtlgen_smr:
  input:
    "{outdir}/results/{gwas}/gwas_sumstat/{gwas}.cleaned.cojo",
    rules.download_eqtlgen.output
  output:
    "{outdir}/results/{gwas}/smr/eqtlgen/{gwas}_smr_eqtlgen_chr{chr}.smr"
  conda:
    "../envs/main.yaml"
  params:
    resdir=resdir
  log:
    "{outdir}/logs/run_eqtlgen_smr-{gwas}-chr{chr}.log"
  shell:
    "{params.resdir}/software/smr/smr_linux_x86_64 \
    --bfile {params.resdir}/data/1kg/1KG.Phase3.EUR.MAF_001.chr{wildcards.chr} \
    --gwas-summary {outdir}/results/{wildcards.gwas}/gwas_sumstat/{wildcards.gwas}.cleaned.cojo \
    --beqtl-summary {params.resdir}/data/eqtlgen/cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense \
    --out {outdir}/results/{wildcards.gwas}/smr/eqtlgen/{wildcards.gwas}_smr_eqtlgen_chr{wildcards.chr} > {log} 2>&1"

rule run_eqtlgen_smr_chr:
    input:
      lambda w: expand("{outdir}/results/{gwas}/smr/eqtlgen/{gwas}_smr_eqtlgen_chr{chr}.smr", gwas=w.gwas, chr=chromosomes, outdir={outdir})
    output:
      touch("{outdir}/results/{gwas}/checks/eqtlgen_smr_all_chr.done")

# Format SMR eQTLGen results
rule format_eqtlgen_smr:
  input:
    "{outdir}/results/{gwas}/checks/eqtlgen_smr_all_chr.done",
    rules.download_biomart.output
  output:
    "{outdir}/results/{gwas}/smr/eqtlgen/{gwas}_smr_eqtlgen_GW.txt.gz"
  conda:
    "../envs/main.yaml"
  params:
    config_file=config['config_file']
  log:
    "{outdir}/logs/format_eqtlgen_smr-{gwas}.log"
  shell:
    "Rscript scripts/format_eqtlgen_smr.R \
      --gwas {wildcards.gwas} \
      --config_file {params.config_file} > {log} 2>&1"

####
# Run SMR using pQTL
####

rule run_rosmap_smr:
  input:
    "{outdir}/results/{gwas}/gwas_sumstat/{gwas}.cleaned.cojo",
    rules.format_rosmap_smr_data.output
  output:
    "{outdir}/results/{gwas}/smr/rosmap/{gwas}_smr_rosmap_chr{chr}.smr"
  conda:
    "../envs/main.yaml"
  params:
    resdir=resdir
  log:
    "{outdir}/logs/run_rosmap_smr-{gwas}-chr{chr}.log"
  shell:
    "{params.resdir}/software/smr/smr_linux_x86_64 \
    --bfile {params.resdir}/data/1kg/1KG.Phase3.EUR.MAF_001.chr{wildcards.chr} \
    --gwas-summary {outdir}/results/{wildcards.gwas}/gwas_sumstat/{wildcards.gwas}.cleaned.cojo \
    --beqtl-summary {params.resdir}/data/rosmap_smr/ROSMAP.n376.pQTL.MatrixQTL.txt.besd \
    --out {outdir}/results/{wildcards.gwas}/smr/rosmap/{wildcards.gwas}_smr_rosmap_chr{wildcards.chr} > {log} 2>&1"

rule run_rosmap_smr_chr:
    input:
      lambda w: expand("{outdir}/results/{gwas}/smr/rosmap/{gwas}_smr_rosmap_chr{chr}.smr", gwas=w.gwas, chr=chromosomes, outdir={outdir})
    output:
      touch("{outdir}/results/{gwas}/checks/rosmap_smr_all_chr.done")

# Format rosmap smr results
rule process_rosmap_smr:
  input:
    "{outdir}/results/{gwas}/checks/rosmap_smr_all_chr.done",
    rules.download_biomart.output
  output:
    "{outdir}/results/{gwas}/smr/rosmap/{gwas}_smr_rosmap_GW.txt.gz"
  conda:
    "../envs/main.yaml"
  params:
    config_file=config['config_file']
  log:
    "{outdir}/logs/process_rosmap_smr-{gwas}.log"
  shell:
    "Rscript scripts/process_rosmap_smr.R \
      --gwas {wildcards.gwas} \
      --config_file {params.config_file} > {log} 2>&1"


###
# Run SMR with MetaBrain
###

rule smr_analysis_MetaBrain:
  resources:
    mem_mb=15000
  input:
    rules.prep_1kg.output,
    "{outdir}/results/{gwas}/gwas_sumstat/{gwas}.cleaned.cojo",
    rules.format_metabrain_esi.output
  output:
    "{outdir}/results/{gwas}/smr/metabrain/{tissue}/{gwas}_smr_metabrain_{tissue}_chr{chr}.smr"
  conda:
    "../envs/main.yaml"
  params:
    resdir=resdir
  log:
    "{outdir}/logs/smr_analysis_MetaBrain-{gwas}-{tissue}-chr{chr}.log"
  shell:
    "{params.resdir}/software/smr/smr_linux_x86_64 \
      --bfile {params.resdir}/data/1kg/1KG.Phase3.EUR.MAF_001.chr{wildcards.chr} \
      --gwas-summary {outdir}/results/{wildcards.gwas}/gwas_sumstat/{wildcards.gwas}.cleaned.cojo \
      --beqtl-summary {params.resdir}/data/MetaBrain/{wildcards.tissue}/2020-05-26-{wildcards.tissue}-EUR-{wildcards.chr}-SMR-besd \
      --out {outdir}/results/{wildcards.gwas}/smr/metabrain/{wildcards.tissue}/{wildcards.gwas}_smr_metabrain_{wildcards.tissue}_chr{wildcards.chr} \
      --thread-num 1 > {log} 2>&1"

rule run_smr_analysis_MetaBrain_Basalganglia:
    input:
      lambda w: expand("{outdir}/results/{gwas}/smr/metabrain/Basalganglia/{gwas}_smr_metabrain_Basalganglia_chr{chr}.smr", gwas=w.gwas, chr=chromosomes, outdir={outdir})
    output:
      touch('{outdir}/results/{gwas}/checks/metabrain_smr_basalganglia_all_chr.done')

rule run_smr_analysis_MetaBrain_Cerebellum:
    input:
      lambda w: expand("{outdir}/results/{gwas}/smr/metabrain/Cerebellum/{gwas}_smr_metabrain_Cerebellum_chr{chr}.smr", gwas=w.gwas, chr=chromosomes, outdir={outdir})
    output:
      touch('{outdir}/results/{gwas}/checks/metabrain_smr_cerebellum_all_chr.done')

rule run_smr_analysis_MetaBrain_Cortex:
    input:
      lambda w: expand("{outdir}/results/{gwas}/smr/metabrain/Cortex/{gwas}_smr_metabrain_Cortex_chr{chr}.smr", gwas=w.gwas, chr=chromosomes, outdir={outdir})
    output:
      touch('{outdir}/results/{gwas}/checks/metabrain_smr_cortex_all_chr.done')

rule run_smr_analysis_MetaBrain_Hippocampus:
    input:
      lambda w: expand("{outdir}/results/{gwas}/smr/metabrain/Hippocampus/{gwas}_smr_metabrain_Hippocampus_chr{chr}.smr", gwas=w.gwas, chr=chromosomes, outdir={outdir})
    output:
      touch('{outdir}/results/{gwas}/checks/metabrain_smr_hippocampus_all_chr.done')

rule run_smr_analysis_MetaBrain_Spinalcord:
    input:
      lambda w: expand("{outdir}/results/{gwas}/smr/metabrain/Spinalcord/{gwas}_smr_metabrain_Spinalcord_chr{chr}.smr", gwas=w.gwas, chr=chromosomes, outdir={outdir})
    output:
      touch('{outdir}/results/{gwas}/checks/metabrain_smr_spinalcord_all_chr.done')

rule run_smr_analysis_MetaBrain_all_tissue:
    input:
      lambda w: expand("{outdir}/results/{gwas}/smr/metabrain/{tissue}/{gwas}_smr_metabrain_{tissue}_chr{chr}.smr", gwas=w.gwas, chr=chromosomes, tissue=['Basalganglia','Cerebellum','Cortex','Hippocampus','Spinalcord'], outdir={outdir})
    output:
      touch('{outdir}/results/{gwas}/checks/metabrain_smr_all_tissue_all_chr.done')

# Format MetaBrain SMR results
metabrain_output = list()

if config["smr_expression_panel_metabrain_basalganglia"] == "T":
    metabrain_output.append("{outdir}/results/{gwas}/checks/metabrain_smr_basalganglia_all_chr.done")

if config["smr_expression_panel_metabrain_cerebellum"] == "T":
    metabrain_output.append("{outdir}/results/{gwas}/checks/metabrain_smr_cerebellum_all_chr.done")

if config["smr_expression_panel_metabrain_cortex"] == "T":
    metabrain_output.append("{outdir}/results/{gwas}/checks/metabrain_smr_cortex_all_chr.done")

if config["smr_expression_panel_metabrain_hippocampus"] == "T":
    metabrain_output.append("{outdir}/results/{gwas}/checks/metabrain_smr_hippocampus_all_chr.done")

if config["smr_expression_panel_metabrain_spinalcord"] == "T":
    metabrain_output.append("{outdir}/results/{gwas}/checks/metabrain_smr_spinalcord_all_chr.done")

rule format_metabrain_smr:
  input:
    metabrain_output,
    "scripts/format_metabrain_smr.R",
    rules.download_biomart.output
  output:
    "{outdir}/results/{gwas}/smr/metabrain/{gwas}_smr_metabrain_GW.txt.gz"
  conda:
    "../envs/main.yaml"
  params:
    config_file=config['config_file']
  log:
    "{outdir}/logs/format_metabrain_smr-{gwas}.log"
  shell:
    "Rscript scripts/format_metabrain_smr.R \
      --gwas {wildcards.gwas} \
      --config_file {params.config_file} > {log} 2>&1"
