################
# GWAS to gene
################

#########
# Prepare required software and reference data
#########

####
# Download liftover and GRCh 37 to 38 and 36 track
####

# Download liftover
rule install_liftover:
  output:
    touch("resources/software/install_liftover.done")
  conda:
    "../envs/GenoFunc.yaml"
  shell:
    "mkdir -p resources/software/liftover/; wget --no-check-certificate -O resources/software/liftover/liftover https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver"

# Download liftover track
rule download_liftover_track:
  output:
    touch("resources/software/download_liftover_track.done")
  conda:
    "../envs/GenoFunc.yaml"
  shell:
    "mkdir -p resources/data/liftover/; wget --no-check-certificate -O resources/data/liftover/hg19ToHg38.over.chain.gz ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz; wget --no-check-certificate -O resources/data/liftover/hg19ToHg18.over.chain.gz ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg18.over.chain.gz"

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
    touch("resources/data/prep_1kg.done")
  conda:
    "../envs/GenoFunc.yaml"
  shell:
    "Rscript scripts/prep_1kg.R"

# Create GW merged version of 1kg refrence data
# Not finished yet. Note. There are variants with duplicate IDs and missing IDs
# We might want to swap over to using CHR:BP:A1:A2 IDs?
rule merge_1kg_GW:
  input:
    rules.prep_1kg.output
  output:
    "resources/data/1kg/1KGPhase3.w_hm3.GW.bed"
  conda:
    "../envs/GenoFunc.yaml"
  shell:
    "for chr in $(seq 1 22);do plink \
       --bfile resources/data/1kg/1KG.Phase3.EUR.MAF_001.chr${chr} \
       --make-bed \
       --set-missing-var-ids @:#\$1:\$2 \
       --out resources/data/1kg/1KG.Phase3.EUR.MAF_001.ID.chr${chr}; done; \
       ls resources/data/1kg/1KG.Phase3.EUR.MAF_001.ID.chr*.bed | sed -e 's/\.bed//g' > resources/data/1kg/merge_list.txt; \
     plink \
       --merge-list resources/data/1kg/merge_list.txt \
       --make-bed \
       --out resources/data/1kg/1KG.Phase3.EUR.MAF_001.GW; \
     rm resources/data/1kg/merge_list.txt"
     
####
# Download software required for TWAS-related analysis
####

# Install fusion
rule install_fusion:
  output:
    directory("resources/software/fusion/")
  conda:
    "../envs/GenoFunc.yaml"
  shell:
    "git clone git@github.com:gusevlab/fusion_twas.git {output}"

# Download plink2R
rule download_plink2R:
  input:
    rules.install_fusion.output
  output:
    "resources/software/plink2R/plink2R-master/data.bed"
  conda:
    "../envs/GenoFunc.yaml"
  shell:
    "wget -O resources/software/plink2R/master.zip https://github.com/gabraham/plink2R/archive/master.zip; unzip resources/software/plink2R/master.zip -d resources/software/plink2R"

# Install plink2R
rule install_plink2R:
  input: 
    rules.download_plink2R.output
  output:
    touch("resources/software/install_plink2R")
  conda:
    "../envs/GenoFunc.yaml"
  shell:
    "Rscript -e 'install.packages(\"resources/software/plink2R/plink2R-master/plink2R/\",repos=NULL)'"

# Install focus
rule install_focus:
  conda:
    "../envs/GenoFunc.yaml"
  output: 
    touch("resources/software/pyfocus")
  shell:
    "pip install pyfocus==0.6.10 --user"

# Install SNP-weights pipeline repo
rule install_snp_weight_pipe:
  output:
    directory("resources/software/Calculating-FUSION-TWAS-weights-pipeline/")
  conda:
    "../envs/GenoFunc.yaml"
  shell:
    "git clone git@github.com:opain/Calculating-FUSION-TWAS-weights-pipeline.git {output}"

####
# Install cmapR R package
####

rule install_cmapR:
  output:
    "resources/software/install_cmapR.done"
  conda: 
    "../envs/GenoFunc.yaml"
  shell:
    "Rscript scripts/install_cmapR.R"


####
# Dowload data for TWAS related analysis
####

# Download PsychENCODE SNP-weights
rule download_psychENCODE_weights:
  output:
    directory("resources/data/psychencode_data/SNP-weights/PEC_TWAS_weights")
  conda:
    "../envs/GenoFunc.yaml"
  shell:
    "wget -O resources/data/psychencode_data/PEC_TWAS_weights.tar.gz http://resource.psychencode.org/Datasets/Derived/PEC_TWAS_weights.tar.gz; mkdir -p resources/data/psychencode_data/SNP-weights/PEC_TWAS_weights; tar xvzf resources/data/psychencode_data/PEC_TWAS_weights.tar.gz -C resources/data/psychencode_data/SNP-weights/PEC_TWAS_weights; rm resources/data/psychencode_data/PEC_TWAS_weights.tar.gz"

# Format PsychENCODE SNP-weights
rule format_psychencode:
  input:
    psychencode_data=rules.download_psychENCODE_weights.output, 
    weights_pipe=rules.install_snp_weight_pipe.output
  output:
    touch("resources/data/psychencode_data/format_psychencode.done")
  conda: 
    "../envs/GenoFunc.yaml"
  shell:
    "Rscript scripts/format_psychENCODE.R"

####
# Download LINCS 1000 data
####

# This should be done in advance as the file is very large.
# The location of the file should be specific in the config.yaml file.
# I have downloaded the the following files: 
# - https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/level5/level5_beta_all_n1201944x12328.gctx
# - https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/siginfo_beta.txt

##########
# Analyse GWAS summary statistics
##########

# For the time being, assume the GWAS sumstats are in Rosalind format
##
# QC and format GWAS summary statistics
##

import pandas as pd
gwas_list_df = pd.read_table(config["gwas_list"], sep=' ')
gwas_list_df_eur = gwas_list_df.loc[gwas_list_df['population'] == 'EUR']

rule sumstat_prep:
  input:
    rules.prep_1kg.output
  output:
    "resources/data/gwas_sumstat/{gwas}/{gwas}.cleaned.gz"
  conda:
    "../envs/GenoFunc.yaml"
  params:
    population= lambda w: gwas_list_df.loc[gwas_list_df['name'] == "{}".format(w.gwas), 'population'].iloc[0],
    path= lambda w: gwas_list_df.loc[gwas_list_df['name'] == "{}".format(w.gwas), 'path'].iloc[0]
  shell:
    "Rscript scripts/sumstat_cleaner.R \
      --sumstats {params.path} \
      --ref_chr resources/data/1kg/1KG.Phase3.EUR.MAF_001.chr \
      --output resources/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}.cleaned"
    
rule run_sumstat_prep:
  input: expand("resources/data/gwas_sumstat/{gwas}/{gwas}.cleaned.gz", gwas=gwas_list_df['name'])

###
# Munge sumstats
###

# munge sumstats using FOCUS munge function
rule focus_munge:
  input:
    premunged="resources/data/gwas_sumstat/{gwas}/{gwas}.cleaned.gz",
    focus=rules.install_focus.output
  output:
    "resources/data/gwas_sumstat/{gwas}/{gwas}.cleaned.munged.sumstats.gz"
  conda:
    "../envs/GenoFunc.yaml"
  shell:
    "focus munge {input.premunged} --output resources/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}.cleaned.munged"

# Calculate median effective sample size
# FUSION requires this parameter to be specified despite having the N column in the sumstats
rule retrieve_N:
  input:
    "resources/data/gwas_sumstat/{gwas}/{gwas}.cleaned.munged.sumstats.gz"
  output:
    "resources/data/gwas_sumstat/{gwas}/{gwas}.cleaned.munged.median_N.txt"
  conda: 
    "../envs/GenoFunc.yaml"
  shell:
    "Rscript scripts/median_n.R --munged {input} --out {output}"

###
# Run TWAS
###

# Run twas using psychENCODE SNP-weights
rule run_psychencode_twas:
  resources: mem_mb=20000 
  input:
    sumstats="resources/data/gwas_sumstat/{gwas}/{gwas}.cleaned.munged.sumstats.gz",
    neff_txt="resources/data/gwas_sumstat/{gwas}/{gwas}.cleaned.munged.median_N.txt",
    fusion=rules.install_fusion.output,
    plink2R=rules.install_plink2R.output,
    format_psychencode=rules.format_psychencode.output,
    prep_1kg=rules.prep_1kg.output
  output:
    "results/{gwas}/twas/psychencode/{gwas}_twas_psychencode_chr{chr}"
  conda: 
    "../envs/GenoFunc.yaml"
  shell:
    "N=$(cat resources/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}.cleaned.munged.median_N.txt); Rscript resources/software/fusion/FUSION.assoc_test.R \
    --sumstats {input.sumstats} \
    --weights resources/data/psychencode_data/PEC_TWAS_weights/PEC_TWAS_weights.pos \
    --weights_dir resources/data/psychencode_data/PEC_TWAS_weights \
    --ref_ld_chr resources/data/1kg/1KG.Phase3.EUR.MAF_001.chr \
    --out {output} \
    --chr {wildcards.chr} \
    --coloc_P 1e-3 \
    --GWASN ${{N}}"

rule psychencode_twas_all_chr:
    input: 
      lambda w: expand("results/{gwas}/twas/psychencode/{gwas}_twas_psychencode_chr{chr}", gwas=w.gwas, chr=range(1, 23))
    output: 
      touch("results/{gwas}/checks/psychencode_twas_all_chr.done")

####
# Run So et al. analysis
####

rule run_psychencode_twas:
  resources: mem_mb=20000 
  input:
    results/{gwas}/checks/psychencode_twas_all_chr.done,
    rules.install_cmapR.output
  output:
    ""
  conda: 
    "../envs/GenoFunc.yaml"
  shell:
    ""

###
# To do:
# Implement So et al script
# Implement CLUE.io analysis
# Implement analysis using TWAS-GSEA
# Input latest GWAS summary statistics for brain-related disorders
###



