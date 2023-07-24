####
# Download SMR
####

rule download_smr:
  output:
    "resources/software/smr/smr_Linux"
  conda:
    "../envs/main.yaml"
  shell:
    "mkdir -p resources/software/smr; wget -O resources/software/smr/smr_Linux.zip https://yanglab.westlake.edu.cn/software/smr/download/smr_Linux.zip; unzip resources/software/smr/smr_Linux.zip -d resources/software/smr; rm resources/software/smr/smr_Linux.zip"

####
# Format ROSMAP SMR data
####

rule format_rosmap_smr_data:
  input:
    rules.download_smr.output,
    rules.prep_1kg.output
  output:
    "resources/data/rosmap_smr/ROSMAP.n376.pQTL.MatrixQTL.txt.besd.epi"
  conda:
    "../envs/main.yaml"
  params:
    rosmap_smr= config["rosmap_smr"],
  shell:
    "Rscript scripts/format_rosmap_smr_data.R \
      --rosmap {params.rosmap_smr}"

####
# Download PsychENCODE data for SMR
####

rule download_psychencode_smr:
  output:
    directory("resources/data/psychencode_smr/PsychENCODE_cis_eqtl_HCP100_summary/")
  conda:
    "../envs/main.yaml"
  shell:
    "mkdir -p resources/data/psychencode_smr; wget --no-check-certificate -O resources/data/psychencode_smr/PsychENCODE_cis_eqtl_HCP100_summary.tar.gz http://cnsgenomics.com/data/SMR/PsychENCODE_cis_eqtl_HCP100_summary.tar.gz; tar -xvzf resources/data/psychencode_smr/PsychENCODE_cis_eqtl_HCP100_summary.tar.gz -C resources/data/psychencode_smr; rm resources/data/psychencode_smr/PsychENCODE_cis_eqtl_HCP100_summary.tar.gz"

##
# Download MetaBrain data in SMR format
##

# Basalganglia
rule download_MetaBrain_Basalganglia:
  output: 
    directory("resources/data/MetaBrain/Basalganglia")
  conda: 
    "../envs/main.yaml"
  shell: 
    "mkdir -p resources/data/MetaBrain/Basalganglia; wget -O resources/data/MetaBrain/Basalganglia/2020-05-26-Basalganglia-EUR-smr.zip https://download.metabrain.nl/2020-05-26-CisEQTLSummaryStats/2020-05-26-Basalganglia-EUR/2020-05-26-Basalganglia-EUR-smr.zip; unzip -d resources/data/MetaBrain/Basalganglia/ resources/data/MetaBrain/Basalganglia/2020-05-26-Basalganglia-EUR-smr.zip; rm resources/data/MetaBrain/Basalganglia/2020-05-26-Basalganglia-EUR-smr.zip"

# Cerebellum
rule download_MetaBrain_Cerebellum:
  output: 
    directory("resources/data/MetaBrain/Cerebellum")
  conda: 
    "../envs/main.yaml"
  shell: 
    "mkdir -p resources/data/MetaBrain/Cerebellum; wget -O resources/data/MetaBrain/Cerebellum/2020-05-26-Cerebellum-EUR-smr.zip https://download.metabrain.nl/2020-05-26-CisEQTLSummaryStats/2020-05-26-Cerebellum-EUR/2020-05-26-Cerebellum-EUR-smr.zip; unzip -d resources/data/MetaBrain/Cerebellum/ resources/data/MetaBrain/Cerebellum/2020-05-26-Cerebellum-EUR-smr.zip; rm resources/data/MetaBrain/Cerebellum/2020-05-26-Cerebellum-EUR-smr.zip"

# Cortex
rule download_MetaBrain_Cortex:
  output: 
    directory("resources/data/MetaBrain/Cortex")
  conda: 
    "../envs/main.yaml"
  shell: 
    "mkdir -p resources/data/MetaBrain/Cortex; wget -O resources/data/MetaBrain/Cortex/2020-05-26-Cortex-EUR-smr.zip https://download.metabrain.nl/2020-05-26-CisEQTLSummaryStats/2020-05-26-Cortex-EUR/2020-05-26-Cortex-EUR-smr.zip; unzip -d resources/data/MetaBrain/Cortex/ resources/data/MetaBrain/Cortex/2020-05-26-Cortex-EUR-smr.zip; rm resources/data/MetaBrain/Cortex/2020-05-26-Cortex-EUR-smr.zip"

# Hippocampus
rule download_MetaBrain_Hippocampus:
  output: 
    directory("resources/data/MetaBrain/Hippocampus")
  conda: 
    "../envs/main.yaml"
  shell: 
    "mkdir -p resources/data/MetaBrain/Hippocampus; wget -O resources/data/MetaBrain/Hippocampus/2020-05-26-Hippocampus-EUR-smr.zip https://download.metabrain.nl/2020-05-26-CisEQTLSummaryStats/2020-05-26-Hippocampus-EUR/2020-05-26-Hippocampus-EUR-smr.zip; unzip -d resources/data/MetaBrain/Hippocampus/ resources/data/MetaBrain/Hippocampus/2020-05-26-Hippocampus-EUR-smr.zip; rm resources/data/MetaBrain/Hippocampus/2020-05-26-Hippocampus-EUR-smr.zip"

# Spinalcord
rule download_MetaBrain_Spinalcord:
  output: 
    directory("resources/data/MetaBrain/Spinalcord")
  conda: 
    "../envs/main.yaml"
  shell: 
    "mkdir -p resources/data/MetaBrain/Spinalcord; wget -O resources/data/MetaBrain/Spinalcord/2020-05-26-Spinalcord-EUR-smr.zip https://download.metabrain.nl/2020-05-26-CisEQTLSummaryStats/2020-05-26-Spinalcord-EUR/2020-05-26-Spinalcord-EUR-smr.zip; unzip -d resources/data/MetaBrain/Spinalcord/ resources/data/MetaBrain/Spinalcord/2020-05-26-Spinalcord-EUR-smr.zip; rm resources/data/MetaBrain/Spinalcord/2020-05-26-Spinalcord-EUR-smr.zip"

rule download_MetaBrain_all:
  input: 
    rules.download_MetaBrain_Basalganglia.output,
    rules.download_MetaBrain_Cerebellum.output,
    rules.download_MetaBrain_Cortex.output,
    rules.download_MetaBrain_Hippocampus.output,
    rules.download_MetaBrain_Spinalcord.output
  output:
    touch('resources/data/MetaBrain_download.out')

# Update variant IDs in MetaBrain SMR files
rule format_metabrain_esi:
  output: 
    touch("resources/data/MetaBrain/format_MetaBrain_esi.out")
  conda: 
    "../envs/main.yaml"
  shell: 
    "Rscript scripts/format_metabrain_esi.R"

# Download eQTLGen data in SMR format
rule download_eqtlgen:
  output: 
    touch("resources/data/eqtlgen.done")
  conda: 
    "../envs/main.yaml"
  shell: 
    "mkdir resources/data/eqtlgen; wget -O resources/data/eqtlgen/cis-eQTL-SMR_20191212.tar.gz https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/SMR_formatted/cis-eQTL-SMR_20191212.tar.gz; tar -xvzf resources/data/eqtlgen/cis-eQTL-SMR_20191212.tar.gz -C resources/data/eqtlgen/; rm resources/data/eqtlgen/cis-eQTL-SMR_20191212.tar.gz; gunzip resources/data/eqtlgen/*"

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
  shell:
    "Rscript scripts/format_sumstats_smr.R \
      --gwas {wildcards.gwas}"

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
  shell:
    "resources/software/smr/smr_Linux \
    --bfile resources/data/1kg/1KG.Phase3.EUR.MAF_001.chr{wildcards.chr} \
    --gwas-summary {outdir}/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}.cleaned.cojo \
    --beqtl-summary resources/data/psychencode_smr/PsychENCODE_cis_eqtl_HCP100_summary/Gandal_PsychENCODE_eQTL_HCP100+gPCs20_QTLtools \
    --out {outdir}/results/{wildcards.gwas}/smr/psychencode/{wildcards.gwas}_smr_psychencode_chr{wildcards.chr}"

rule run_psychencode_smr_chr:
    input: 
      lambda w: expand("{outdir}/results/{gwas}/smr/psychencode/{gwas}_smr_psychencode_chr{chr}.smr", gwas=w.gwas, chr=range(1, 23))
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
  shell:
    "resources/software/smr/smr_Linux \
    --bfile resources/data/1kg/1KG.Phase3.EUR.MAF_001.chr{wildcards.chr} \
    --gwas-summary {outdir}/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}.cleaned.cojo \
    --beqtl-summary resources/data/eqtlgen/cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense \
    --out {outdir}/results/{wildcards.gwas}/smr/eqtlgen/{wildcards.gwas}_smr_eqtlgen_chr{wildcards.chr}"

rule run_eqtlgen_smr_chr:
    input: 
      lambda w: expand("{outdir}/results/{gwas}/smr/eqtlgen/{gwas}_smr_eqtlgen_chr{chr}.smr", gwas=w.gwas, chr=range(1, 23))
    output: 
      touch("{outdir}/results/{gwas}/checks/eqtlgen_smr_all_chr.done")
      
# Format SMR eQTLGen results
rule format_eqtlgen_smr:
  input:
    "{outdir}/results/{gwas}/checks/eqtlgen_smr_all_chr.done"
  output:
    "{outdir}/results/{gwas}/smr/eqtlgen/{gwas}_smr_eqtlgen_GW.txt.gz"
  conda: 
    "../envs/main.yaml"
  shell:
    "Rscript scripts/format_eqtlgen_smr.R \
    --gwas {wildcards.gwas}"

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
  shell:
    "resources/software/smr/smr_Linux \
    --bfile resources/data/1kg/1KG.Phase3.EUR.MAF_001.chr{wildcards.chr} \
    --gwas-summary {outdir}/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}.cleaned.cojo \
    --beqtl-summary resources/data/rosmap_smr/ROSMAP.n376.pQTL.MatrixQTL.txt.besd \
    --out {outdir}/results/{wildcards.gwas}/smr/rosmap/{wildcards.gwas}_smr_rosmap_chr{wildcards.chr}"

rule run_rosmap_smr_chr:
    input: 
      lambda w: expand("{outdir}/results/{gwas}/smr/rosmap/{gwas}_smr_rosmap_chr{chr}.smr", gwas=w.gwas, chr=range(1, 23))
    output: 
      touch("{outdir}/results/{gwas}/checks/rosmap_smr_all_chr.done")

# Format rosmap smr results
rule process_rosmap_smr:
  input:
    "{outdir}/results/{gwas}/checks/rosmap_smr_all_chr.done"
  output:
    "{outdir}/results/{gwas}/smr/rosmap/{gwas}_smr_rosmap_GW.txt.gz"
  conda: 
    "../envs/main.yaml"
  shell:
    "Rscript scripts/process_rosmap_smr.R \
    --gwas {wildcards.gwas}"


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
  shell:
    "resources/software/smr/smr_Linux \
      --bfile resources/data/1kg/1KG.Phase3.EUR.MAF_001.chr{wildcards.chr} \
      --gwas-summary {outdir}/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}.cleaned.cojo \
      --beqtl-summary resources/data/MetaBrain/{wildcards.tissue}/2020-05-26-{wildcards.tissue}-EUR-{wildcards.chr}-SMR-besd \
      --out {outdir}/results/{wildcards.gwas}/smr/metabrain/{wildcards.tissue}/{wildcards.gwas}_smr_metabrain_{wildcards.tissue}_chr{wildcards.chr} \
      --thread-num 1"

rule run_smr_analysis_MetaBrain_Basalganglia:
    input: 
      lambda w: expand("{outdir}/results/{gwas}/smr/metabrain/Basalganglia/{gwas}_smr_metabrain_Basalganglia_chr{chr}.smr", gwas=w.gwas, chr=range(1, 23))
    output: 
      touch('{outdir}/results/{gwas}/checks/metabrain_smr_basalganglia_all_chr.done')

rule run_smr_analysis_MetaBrain_Cerebellum:
    input: 
      lambda w: expand("{outdir}/results/{gwas}/smr/metabrain/Cerebellum/{gwas}_smr_metabrain_Cerebellum_chr{chr}.smr", gwas=w.gwas, chr=range(1, 23))
    output: 
      touch('{outdir}/results/{gwas}/checks/metabrain_smr_cerebellum_all_chr.done')

rule run_smr_analysis_MetaBrain_Cortex:
    input: 
      lambda w: expand("{outdir}/results/{gwas}/smr/metabrain/Cortex/{gwas}_smr_metabrain_Cortex_chr{chr}.smr", gwas=w.gwas, chr=range(1, 23))
    output: 
      touch('{outdir}/results/{gwas}/checks/metabrain_smr_cortex_all_chr.done')

rule run_smr_analysis_MetaBrain_Hippocampus:
    input: 
      lambda w: expand("{outdir}/results/{gwas}/smr/metabrain/Hippocampus/{gwas}_smr_metabrain_Hippocampus_chr{chr}.smr", gwas=w.gwas, chr=range(1, 23))
    output: 
      touch('{outdir}/results/{gwas}/checks/metabrain_smr_hippocampus_all_chr.done')

rule run_smr_analysis_MetaBrain_Spinalcord:
    input: 
      lambda w: expand("{outdir}/results/{gwas}/smr/metabrain/Spinalcord/{gwas}_smr_metabrain_Spinalcord_chr{chr}.smr", gwas=w.gwas, chr=range(1, 23))
    output: 
      touch('{outdir}/results/{gwas}/checks/metabrain_smr_spinalcord_all_chr.done')
      
rule run_smr_analysis_MetaBrain_all_tissue:
    input: 
      lambda w: expand("{outdir}/results/{gwas}/smr/metabrain/{tissue}/{gwas}_smr_metabrain_{tissue}_chr{chr}.smr", gwas=w.gwas, chr=range(1, 23), tissue=['Basalganglia','Cerebellum','Cortex','Hippocampus','Spinalcord'])
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
    metabrain_output
  output:
    "{outdir}/results/{gwas}/smr/metabrain/{gwas}_smr_metabrain_GW.txt.gz"
  conda: 
    "../envs/main.yaml"
  shell:
    "Rscript scripts/format_metabrain_smr.R \
    --gwas {wildcards.gwas}"

