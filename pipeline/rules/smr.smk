####
# Download SMR
####

rule download_smr:
  output:
    f"{resdir}/software/smr/smr_linux_x86_64"
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
  conda:
    "../envs/main.yaml"
  params:
    rosmap_smr= config["rosmap_smr"],
  log:
    f"{resdir}/logs/format_rosmap_smr_data.log"
  shell:
    "Rscript scripts/format_rosmap_smr_data.R \
      --rosmap {params.rosmap_smr} > {log} 2>&1"

####
# Download PsychENCODE data for SMR
####

rule download_psychencode_smr:
  output:
    directory(f"{resdir}/data/psychencode_smr/PsychENCODE_cis_eqtl_HCP100_summary/")
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
  conda:
    "../envs/main.yaml"
  log:
    f"{resdir}/logs/format_metabrain_esi.log"
  shell:
    "Rscript scripts/format_metabrain_esi.R > {log} 2>&1"

# Download eQTLGen data in SMR format
rule download_eqtlgen:
  output:
    touch(f"{resdir}/data/eqtlgen.done")
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

##########
# Analyse GWAS summary statistics
##########

####
# Format sumstats for SMR
####

rule format_sumstats_smr:
  input:
    "{outdir}/data/gwas_sumstat/{gwas}/{gwas}.cleaned.gz"
  output:
    "{outdir}/data/gwas_sumstat/{gwas}/{gwas}.cleaned.cojo"
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
    "{outdir}/data/gwas_sumstat/{gwas}/{gwas}.cleaned.cojo",
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
    --gwas-summary {outdir}/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}.cleaned.cojo \
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
    "{outdir}/data/gwas_sumstat/{gwas}/{gwas}.cleaned.cojo",
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
    --gwas-summary {outdir}/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}.cleaned.cojo \
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
    "{outdir}/data/gwas_sumstat/{gwas}/{gwas}.cleaned.cojo",
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
    --gwas-summary {outdir}/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}.cleaned.cojo \
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
    "{outdir}/data/gwas_sumstat/{gwas}/{gwas}.cleaned.cojo",
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
      --gwas-summary {outdir}/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}.cleaned.cojo \
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

