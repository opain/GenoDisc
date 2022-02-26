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

####
# Format LINCS data
####
# Subset LINCS data to compounds with dose=10uM and time=24 hour 
# Convert Entrez IDs to ENSEMBL IDs

rule subset_lincs:
  input:
    rules.install_cmapR.output
  output:
    "resources/data/lincs/subset_lincs.done"
  conda: 
    "../envs/GenoFunc.yaml"
  params:
    lincs_level5= config["lincs_level5"],
    lincs_siginfo= config["lincs_siginfo"]
  shell:
    "Rscript scripts/subset_lincs.R \
      --lincs_level5_path {params.lincs_level5} \
      --lincs_siginfo_path {params.lincs_siginfo}"

####
# Download TWAS-GSEA
####

rule install_twas_gsea:
  output:
    directory("resources/software/TWAS-GSEA/")
  conda:
    "../envs/GenoFunc.yaml"
  shell:
    "git clone git@github.com:opain/TWAS-GSEA.git {output}"

####
# Download FeaturePred
####

rule install_feature_pred:
  output:
    directory("resources/software/Predicting-TWAS-features/")
  conda:
    "../envs/GenoFunc.yaml"
  shell:
    "git clone git@github.com:opain/Predicting-TWAS-features.git {output}"

####
# Download pigz
####

rule install_pigz:
  output:
    "resources/software/pigz/pigz/pigz"
  conda:
    "../envs/GenoFunc.yaml"
  shell:
    "wget -O resources/software/pigz.tar.gz https://zlib.net/pigz/pigz.tar.gz; mkdir -p resources/software/pigz; tar xvzf resources/software/pigz.tar.gz -C resources/software/pigz; rm resources/software/pigz.tar.gz; cd resources/software/pigz/pigz; make"

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
# Predict PsychENCODE features into 1kg sample
####

# Modify panel column in .pos file
rule mod_psychencode_pos:
  output:
    "resources/data/psychencode_data/PEC_TWAS_weights/PEC_TWAS_weights_mod.pos"
  conda:
    "../envs/GenoFunc.yaml"
  shell:
    "cut -d' ' -f 1 --complement resources/data/psychencode_data/PEC_TWAS_weights/PEC_TWAS_weights.pos | sed 's/$/ PEC_TWAS_weights/' | sed -e '1s/PEC_TWAS_weights/PANEL/' > resources/data/psychencode_data/PEC_TWAS_weights/PEC_TWAS_weights_mod.pos"

rule feature_pred_psychencode:
  resources: 
    mem_mb=50000,
    cpus=5
  input:
    rules.install_feature_pred.output,
    rules.format_psychencode.output,
    rules.install_pigz.output,
    rules.mod_psychencode_pos.output
  output:
    "resources/data/psychencode_data/PredictedExpression/Reference_Expression/Reference_Expression_PEC_TWAS_weights.txt.gz"
  conda:
    "../envs/GenoFunc.yaml"
  shell:
    "Rscript resources/software/Predicting-TWAS-features/FeaturePred.V2.0.R \
    	--PLINK_prefix_chr resources/data/1kg/1KG.Phase3.EUR.MAF_001.chr \
    	--weights resources/data/psychencode_data/PEC_TWAS_weights/PEC_TWAS_weights_mod.pos \
    	--weights_dir resources/data/psychencode_data/PEC_TWAS_weights \
    	--ref_ld_chr resources/data/1kg/1KG.Phase3.EUR.MAF_001.chr \
    	--targ_pred F \
    	--save_ref_expr T \
    	--plink plink \
    	--ref_maf resources/data/1kg/1KG.Phase3.EUR.MAF_001.chr \
    	--pigz resources/software/pigz/pigz/pigz \
    	--memory 40000 \
      --n_cores 5 \
    	--output resources/data/psychencode_data/PredictedExpression"

# Format expression data for TWAs-GSEA (i.e. remove PANEL from column names)
rule format_pred_psychencode:
  input:
    rules.feature_pred_psychencode.output
  output:
    touch("resources/data/psychencode_data/PredictedExpression/Reference_Expression/format_pred_psychencode.done")
  conda:
    "../envs/GenoFunc.yaml"
  shell:
    "zcat resources/data/psychencode_data/PredictedExpression/Reference_Expression/Reference_Expression_PEC_TWAS_weights.txt.gz | sed -e s/PEC_TWAS_weights.//g | gzip > resources/data/psychencode_data/PredictedExpression/Reference_Expression/Reference_Expression_PEC_TWAS_weights_mod.txt.gz"
    
####
# Install lme4qtl
####
# Note the version of conda was not working in R 4.0.2

rule install_lme4qtl:
  output:
    touch("resources/software/install_lme4qtl.done")
  conda:
    "../envs/GenoFunc.yaml"
  shell:
    "Rscript -e 'devtools::install_github(\"variani/lme4qtl\")'"

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
      --annotate window=10 \
    	--snp-loc resources/data/magma_ref/g1000_eur.bim \
    	--gene-loc resources/data/magma/NCBI37.3.gene.loc \
    	--out resources/data/magma/NCBI37.3"

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

####
# Download LDSC
####

# Install LDSC
rule install_ldsc:
  output:
    directory("resources/software/ldsc/")
  conda:
    "../envs/GenoFunc.yaml"
  shell:
    "git clone git@github.com:bulik/ldsc.git {output}"

# Download LDSC reference data
rule download_ldsc_scores:
  output:
    "resources/software/ldsc/eur_w_ld_chr/10.l2.ldscore.gz"
  conda:
    "../envs/GenoFunc.yaml"
  shell:
    "mkdir -p resources/data/ldsc; wget --no-check-certificate -O resources/data/ldsc/eur_w_ld_chr.tar.bz2 https://data.broadinstitute.org/alkesgroup/LDSCORE/eur_w_ld_chr.tar.bz2; tar -jxvf resources/data/ldsc/eur_w_ld_chr.tar.bz2 -C resources/data/ldsc; rm resources/data/ldsc/eur_w_ld_chr.tar.bz2"

rule download_ldsc_hm3:
  output:
    "resources/software/ldsc/w_hm3.snplist"
  conda:
    "../envs/GenoFunc.yaml"
  shell:
    "mkdir -p resources/data/ldsc; wget --no-check-certificate -O resources/data/ldsc/w_hm3.snplist.bz2 https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2; bunzip2 -d resources/data/ldsc/w_hm3.snplist.bz2"

####
# Format ROSMAP and Banner PWAS data
####

rule format_pwas_data:
  output:
    "resources/data/banner_twas/Banner.n152.fusion.WEIGHTS/train_weights_withN.pos"
  conda:
    "../envs/GenoFunc.yaml"
  params:
    rosmap_fusion= config["rosmap_fusion"],
    banner_fusion= config["banner_fusion"]
  shell:
    "Rscript scripts/format_pwas_data.R \
      --rosmap {params.rosmap_fusion} \
      --banner {params.banner_fusion}"

####
# Download SMR
####

rule download_smr:
  output:
    "resources/software/smr/smr_Linux"
  conda:
    "../envs/GenoFunc.yaml"
  shell:
    "mkdir -p resources/software/smr; wget -O resources/software/smr/smr_Linux.zip https://yanglab.westlake.edu.cn/software/smr/download/smr_Linux.zip; unzip resources/software/smr/smr_Linux.zip -d resources/software/smr; rm resources/software/smr/smr_Linux.zip"

####
# Format ROSMAP SMR data
####

rule format_rosmap_smr_data:
  input:
    rules.download_smr.output
  output:
    "resources/data/rosmap_smr/Banner.n152.fusion.WEIGHTS/train_weights_withN.pos"
  conda:
    "../envs/GenoFunc.yaml"
  params:
    rosmap_smr= config["rosmap_smr"],
  shell:
    "Rscript scripts/format_pwas_data.R \
      --rosmap {params.rosmap_smr}"

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
    population= lambda w: gwas_list_df_eur.loc[gwas_list_df_eur['name'] == "{}".format(w.gwas), 'population'].iloc[0],
    path= lambda w: gwas_list_df_eur.loc[gwas_list_df_eur['name'] == "{}".format(w.gwas), 'path'].iloc[0]
  shell:
    "Rscript scripts/sumstat_cleaner.R \
      --sumstats {params.path} \
      --ref_chr resources/data/1kg/1KG.Phase3.{params.population}.MAF_001.chr \
      --output resources/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}.cleaned"
    
rule run_sumstat_prep:
  input: expand("resources/data/gwas_sumstat/{gwas}/{gwas}.cleaned.gz", gwas=gwas_list_df_eur['name'])

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
# Run LDSC
###

rule ldsc:
  input:
    "resources/data/gwas_sumstat/{gwas}/{gwas}.cleaned.munged.sumstats.gz"
  output:
    "results/{gwas}/twas/ldsc"
  conda:
    "../envs/ldsc.yaml"
  shell:
    "python resources/software/ldsc/ldsc.py \
      --h2 resources/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}.cleaned.munged.sumstats.gz \
      --ref-ld-chr resources/software/ldsc/eur_w_ld_chr/ \
      --w-ld-chr resources/software/ldsc/eur_w_ld_chr/ \
      --out results/{wildcards.gwas}/ldsc/{wildcards.gwas}_ldsc_res"

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

###
# Run PWAS
###

# Run twas using ROSMAP SNP-weights
rule run_rosmap_pwas:
  resources: mem_mb=20000 
  input:
    sumstats="resources/data/gwas_sumstat/{gwas}/{gwas}.cleaned.munged.sumstats.gz",
    neff_txt="resources/data/gwas_sumstat/{gwas}/{gwas}.cleaned.munged.median_N.txt",
    fusion=rules.install_fusion.output,
    plink2R=rules.install_plink2R.output,
    format_psychencode=rules.format_pwas_data.output,
    prep_1kg=rules.prep_1kg.output
  output:
    "results/{gwas}/pwas/rosmap/{gwas}_pwas_rosmap_chr{chr}"
  conda: 
    "../envs/GenoFunc.yaml"
  shell:
    "N=$(cat resources/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}.cleaned.munged.median_N.txt); Rscript resources/software/fusion/FUSION.assoc_test.R \
    --sumstats {input.sumstats} \
    --weights resources/data/rosmap_twas/ROSMAP.n376.fusion.WEIGHTS/train_weights_withN.pos \
    --weights_dir resources/data/rosmap_twas/ROSMAP.n376.fusion.WEIGHTS \
    --ref_ld_chr resources/data/1kg/1KG.Phase3.EUR.MAF_001.chr \
    --out {output} \
    --chr {wildcards.chr} \
    --coloc_P 1e-3 \
    --GWASN ${{N}}"

rule rosmap_pwas_all_chr:
    input: 
      lambda w: expand("results/{gwas}/pwas/rosmap/{gwas}_pwas_rosmap_chr{chr}", gwas=w.gwas, chr=range(1, 23))
    output: 
      touch("results/{gwas}/checks/rosmap_pwas_all_chr.done")

# Run twas using Banner SNP-weights
rule run_banner_pwas:
  resources: mem_mb=20000 
  input:
    sumstats="resources/data/gwas_sumstat/{gwas}/{gwas}.cleaned.munged.sumstats.gz",
    neff_txt="resources/data/gwas_sumstat/{gwas}/{gwas}.cleaned.munged.median_N.txt",
    fusion=rules.install_fusion.output,
    plink2R=rules.install_plink2R.output,
    format_psychencode=rules.format_pwas_data.output,
    prep_1kg=rules.prep_1kg.output
  output:
    "results/{gwas}/pwas/banner/{gwas}_pwas_banner_chr{chr}"
  conda: 
    "../envs/GenoFunc.yaml"
  shell:
    "N=$(cat resources/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}.cleaned.munged.median_N.txt); Rscript resources/software/fusion/FUSION.assoc_test.R \
    --sumstats {input.sumstats} \
    --weights resources/data/banner_twas/Banner.n152.fusion.WEIGHTS/train_weights_withN.pos \
    --weights_dir resources/data/banner_twas/Banner.n152.fusion.WEIGHTS \
    --ref_ld_chr resources/data/1kg/1KG.Phase3.EUR.MAF_001.chr \
    --out {output} \
    --chr {wildcards.chr} \
    --coloc_P 1e-3 \
    --GWASN ${{N}}"

rule banner_pwas_all_chr:
    input: 
      lambda w: expand("results/{gwas}/pwas/banner/{gwas}_pwas_banner_chr{chr}", gwas=w.gwas, chr=range(1, 23))
    output: 
      touch("results/{gwas}/checks/banner_pwas_all_chr.done")

####
# Format sumstats for SMR
####

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

####
# Run SMR using eQTL
####

####
# Run SMR using pQTL
####

####
# Run So et al. analysis
####

# To enable the analysis to be run in parallel for each GWAS, I have split the So et al code into two parts.
# The LINCS data has been slit into 10 batches.
# Step 1 runs the pearson, spearman and KS tests on each batch of pertubagens seperately.
# Step 2 combines the results across all batches, estimates the average ranks, and calculates the p-value.

rule so_cmap_step1:
  input:
    "results/{gwas}/checks/psychencode_twas_all_chr.done",
    rules.install_cmapR.output,
    rules.subset_lincs.output
  output:
    "results/{gwas}/twas/cmap/res_batch_{batch}.RDS"
  conda: 
    "../envs/GenoFunc.yaml"
  shell:
    "Rscript scripts/so_cmap_step1.R \
    --twas {wildcards.gwas} \
    --batch {wildcards.batch}"

rule so_cmap_step1_all_batch:
    input: 
      lambda w: expand("results/{gwas}/twas/cmap/res_batch_{batch}.RDS", gwas=w.gwas, batch=range(1, 11))
    output: 
      touch("results/{gwas}/checks/so_cmap_step1_all_batch.done")


rule so_cmap_step2:
  input:
    "results/{gwas}/checks/so_cmap_step1_all_batch.done"
  output:
    "results/{gwas}/twas/cmap/So_res.csv"
  conda: 
    "../envs/GenoFunc.yaml"
  shell:
    "Rscript scripts/so_cmap_step2.R \
    --twas {wildcards.gwas}"

# Format the so et al results 
rule format_so_results:
  input:
    "results/{gwas}/twas/cmap/So_res.csv",
    rules.download_atc.output
  output:
    "results/{gwas}/twas/cmap/So_res_atc_res.csv"
  conda: 
    "../envs/GenoFunc.yaml"
  params:
    lincs_siginfo= config["lincs_siginfo"]
  shell:
    "Rscript scripts/format_so_results.R \
    --twas {wildcards.gwas} \
    --lincs_siginfo_path {params.lincs_siginfo}"

###
# TWAS-GSEA
###

# Format the CMAP data for TWAS-GSEA
rule format_lincs_for_twas_gsea:
  input:
    rules.subset_lincs.output
  output:
    "resources/data/lincs/lincs_core_subset.txt.gz"
  conda: 
    "../envs/GenoFunc.yaml"
  shell:
    "Rscript scripts/format_lincs_for_twas_gsea.R"

# Concatenate twas results
rule format_twas_for_twas_gsea:
  input:
    rules.subset_lincs.output
  output:
    "results/{gwas}/twas/psychencode/{gwas}_twas_psychencode_GW"
  conda: 
    "../envs/GenoFunc.yaml"
  shell:
    "head -1 results/{wildcards.gwas}/twas/psychencode/{wildcards.gwas}_twas_psychencode_chr1 > results/{wildcards.gwas}/twas/psychencode/{wildcards.gwas}_twas_psychencode_GW; tail -n +2 -q results/{wildcards.gwas}/twas/psychencode/{wildcards.gwas}_twas_psychencode_chr* >> results/{wildcards.gwas}/twas/psychencode/{wildcards.gwas}_twas_psychencode_GW"

# Run TWAS-GSEA
rule run_twas_gsea:
  resources: 
    mem_mb=50000,
    cpus=5
  input:
    rules.install_twas_gsea.output,
    "results/{gwas}/checks/psychencode_twas_all_chr.done",
    rules.feature_pred_psychencode.output,
    rules.format_lincs_for_twas_gsea.output,
    rules.format_twas_for_twas_gsea.output,
    rules.format_pred_psychencode.output,
    rules.install_lme4qtl.output
  output:
    touch("results/{gwas}/twas/cmap/twas_gsea_cmap.done")
  conda:
    "../envs/GenoFunc.yaml"
  shell:
    "Rscript resources/software/TWAS-GSEA/TWAS-GSEA.V1.2.R \
      --twas_results results/{wildcards.gwas}/twas/psychencode/{wildcards.gwas}_twas_psychencode_GW \
    	--pos resources/data/psychencode_data/PEC_TWAS_weights/PEC_TWAS_weights.pos \
    	--prop_file resources/data/lincs/lincs_core_subset.txt.gz \
    	--expression_ref resources/data/psychencode_data/PredictedExpression/Reference_Expression/Reference_Expression_PEC_TWAS_weights_mod.txt.gz \
    	--n_cores 5 \
    	--covar GeneLength,NSNP \
    	--use_alt_id ID \
    	--linear_p_thresh 1 \
    	--save_CorMat F \
    	--min_r2 0.01 \
    	--qqplot F \
    	--output results/{wildcards.gwas}/twas/cmap/twas_gsea_cmap"

# Format TWAS-GSEA results
rule format_twas_gsea_results:
  input:
    "results/{gwas}/twas/cmap/twas_gsea_cmap.done",
    rules.download_atc.output
  output:
    "results/{gwas}/twas/cmap/twas_gsea_res_atc_res.csv"
  conda: 
    "../envs/GenoFunc.yaml"
  params:
    lincs_siginfo= config["lincs_siginfo"]
  shell:
    "Rscript scripts/format_twas_gsea_results.R \
    --twas {wildcards.gwas} \
    --lincs_siginfo_path {params.lincs_siginfo}"

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
    "gzip -f -d resources/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}.cleaned.gz > resources/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}.cleaned; resources/software/magma/magma \
      --bfile resources/data/magma_ref/g1000_eur \
      --pval resources/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}.cleaned use=SNP,P ncol=N \
      --gene-annot resources/data/magma/NCBI37.3.genes.annot \
      --out results/{wildcards.gwas}/magma/magma_gene_level; rm resources/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}.cleaned"

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

# Format the MAGMA results 
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
    "Rscript scripts/format_magma_results.R \
    --gwas {wildcards.gwas}"

####
# Create results report
####
# Create a report for each GWAS

rule create_report:
  input:
    "results/{gwas}/twas/cmap/twas_gsea_res_atc_res.csv",
    "results/{gwas}/twas/cmap/So_res_atc_res.csv",
    "results/{gwas}/magma/magma_drug_targetor_atc_res.csv"
  output:
    "results/{gwas}/reports/{gwas}_report.html"
  conda:
    "../envs/GenoFunc.yaml"
  shell:
    "mkdir -p results/{wildcards.gwas}/reports; Rscript -e \"rmarkdown::render(\'scripts/create_report.Rmd\', \
     output_file = \'../results/{wildcards.gwas}/reports/{wildcards.gwas}_report.html\', \
     params = list(gwas = \'{wildcards.gwas}\'))\""

rule run_create_report:
  input: expand("results/{gwas}/reports/{gwas}_report.html", gwas=gwas_list_df_eur['name'])

###
# To do:
# LDSC SNP-h2
# Manhattan plot
# COJO
# Create lead SNP/locus table
# Include PsychENCODE SMR
# Include MetaBrain SMR
# Include brain PWAS and protein SMR
# Finemapping
# Implement CLUE.io analysis
# Implement GCSC
# Input latest GWAS summary statistics for brain-related disorders
###


