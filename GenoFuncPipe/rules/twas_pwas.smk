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

# We want to download weights and harmonise into a single repository
# First download the weights, then reformat accordingly

# Download PsychENCODE SNP-weights
rule download_psychENCODE_weights:
  output:
    touch("resources/data/download_psychENCODE_weights.done")
  conda:
    "../envs/GenoFunc.yaml"
  shell:
    "mkdir -p resources/data/fusion_snp_weights/psychencode; wget -O resources/data/fusion_snp_weights/psychencode/PEC_TWAS_weights.tar.gz http://resource.psychencode.org/Datasets/Derived/PEC_TWAS_weights.tar.gz; mkdir -p resources/data/fusion_snp_weights/psychencode/psychencode; tar xvzf resources/data/fusion_snp_weights/psychencode/PEC_TWAS_weights.tar.gz -C resources/data/fusion_snp_weights/psychencode/psychencode; rm resources/data/fusion_snp_weights/psychencode/PEC_TWAS_weights.tar.gz"
    
# Format PsychENCODE SNP-weights
rule format_psychencode:
  input:
    psychencode_data=rules.download_psychENCODE_weights.output, 
    weights_pipe=rules.install_snp_weight_pipe.output
  output:
    "resources/data/format_psychencode.done"
  conda: 
    "../envs/GenoFunc.yaml"
  shell:
    "Rscript scripts/format_psychENCODE.R"

# Download FUSION GTEx v8 EUR SNP-weights
# I am using EUR instead of full sample to avoid LD mismatch
gtex_weights=config["gtex_weights"]

rule download_gtex_weights:
  output:
    touch("resources/data/download_fusion_gtex_{weight}_weights.done")
  conda:
    "../envs/GenoFunc.yaml"
  shell:
    "mkdir -p resources/data/fusion_snp_weights/{wildcards.weight}; wget -O resources/data/fusion_snp_weights/GTExv8.EUR.{wildcards.weight}.tar.gz https://s3.us-west-1.amazonaws.com/gtex.v8.fusion/EUR/GTExv8.EUR.{wildcards.weight}.tar.gz; tar xf resources/data/fusion_snp_weights/GTExv8.EUR.{wildcards.weight}.tar.gz -C resources/data/fusion_snp_weights/{wildcards.weight}; rm resources/data/fusion_snp_weights/GTExv8.EUR.{wildcards.weight}.tar.gz; mv resources/data/fusion_snp_weights/{wildcards.weight}/GTExv8.EUR.{wildcards.weight}.pos resources/data/fusion_snp_weights/{wildcards.weight}/{wildcards.weight}.pos; mv resources/data/fusion_snp_weights/{wildcards.weight}/GTExv8.EUR.{wildcards.weight} resources/data/fusion_snp_weights/{wildcards.weight}/{wildcards.weight}"

# Update GTEx v8 P0 and P1 to build GRCh 37
rule update_gtex_coord:
  input:
    "resources/data/download_fusion_gtex_{weight}_weights.done"
  output:
    touch("resources/data/update_gtex_coord_{weight}.done")
  conda:
    "../envs/GenoFunc.yaml"
  shell:
    "Rscript scripts/update_gtex_coord.R \
      --panel {wildcards.weight}"

rule update_gtex_coord_all_panel:
    input: expand("resources/data/update_gtex_coord_{weight}.done", weight=gtex_weights)

# Download FUSION non-GTEx SNP-weights
non_gtex_weights=config["non_gtex_weights"]

rule download_non_gtex_weights:
  output:
    touch("resources/data/download_non_gtex_{weight}_weights.done")
  conda:
    "../envs/GenoFunc.yaml"
  shell:
    "mkdir -p resources/data/fusion_snp_weights/{wildcards.weight}; wget --no-check-certificate -O resources/data/fusion_snp_weights/{wildcards.weight}.tar.bz2 https://data.broadinstitute.org/alkesgroup/FUSION/WGT/{wildcards.weight}.tar.bz2; tar xvjf resources/data/fusion_snp_weights/{wildcards.weight}.tar.bz2 -C resources/data/fusion_snp_weights/{wildcards.weight}; rm resources/data/fusion_snp_weights/{wildcards.weight}.tar.bz2"

# Insert N into non-GTEX SNP-weights
rule insert_n_nongtex:
  input:
    "resources/data/download_non_gtex_{weight}_weights.done"
  output:
    touch("resources/data/insert_n_nongtex_{weight}.done")
  conda:
    "../envs/GenoFunc.yaml"
  shell:
    "Rscript scripts/insert_n_nongtex.R \
      --panel {wildcards.weight}"

rule insert_n_nongtex_all_panel:
    input: expand("resources/data/insert_n_nongtex_{weight}.done", weight=non_gtex_weights)

# Download glist file
rule download_glist:
  output:
    "resources/data/glist-hg19"
  conda: 
    "../envs/GenoFunc.yaml"
  shell:
    "wget -P resources/data/ https://www.cog-genomics.org/static/bin/plink/glist-hg19"

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
# Format the CMAP data for TWAS-GSEA
####

rule format_lincs_for_twas_gsea:
  input:
    rules.subset_lincs.output
  output:
    "resources/data/lincs/lincs_core_subset.txt.gz"
  conda: 
    "../envs/GenoFunc.yaml"
  shell:
    "Rscript scripts/format_lincs_for_twas_gsea.R"

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
    if not os.path.isdir("".join(["resources/data/fusion_snp_weights/",external_weights_id_list[x]])):
      os.system("".join(["mkdir resources/data/fusion_snp_weights/",external_weights_id_list[x]]))
      os.system("".join(["cp -r ",external_weights_path_list[x],"/* resources/data/fusion_snp_weights/", external_weights_id_list[x],"/"]))

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

# Modify panel column in .pos file
rule feature_pred:
  resources: 
    mem_mb=50000,
    cpus=5
  input:
    rules.install_feature_pred.output,
    rules.format_psychencode.output,
    rules.update_gtex_coord_all_panel.input,
    rules.insert_n_nongtex_all_panel.input,
    rules.install_pigz.output
  output:
    "resources/data/predicted_expression/{weight}/Reference_Expression/Reference_Expression_{weight}.txt.gz"
  conda:
    "../envs/GenoFunc.yaml"
  shell:
    "Rscript resources/software/Predicting-TWAS-features/FeaturePred.V2.0.R \
    	--PLINK_prefix_chr resources/data/1kg/1KG.Phase3.EUR.MAF_001.chr \
    	--weights resources/data/fusion_snp_weights/{wildcards.weight}/{wildcards.weight}.pos \
    	--weights_dir resources/data/fusion_snp_weights/{wildcards.weight} \
    	--ref_ld_chr resources/data/1kg/1KG.Phase3.EUR.MAF_001.chr \
    	--targ_pred F \
    	--save_ref_expr T \
    	--save_score F \
    	--plink plink \
    	--ref_maf resources/data/1kg/1KG.Phase3.EUR.MAF_001.chr \
    	--pigz resources/software/pigz/pigz/pigz \
    	--memory 40000 \
      --n_cores 5 \
    	--output resources/data/predicted_expression/{wildcards.weight}"

# Format expression data for TWAs-GSEA (i.e. remove PANEL from column names)
rule format_pred:
  input:
    "resources/data/predicted_expression/{weight}/Reference_Expression/Reference_Expression_{weight}.txt.gz"
  output:
    touch("resources/data/predicted_expression/format_pred_{weight}.done")
  conda:
    "../envs/GenoFunc.yaml"
  shell:
    "zcat resources/data/predicted_expression/{wildcards.weight}/Reference_Expression/Reference_Expression_{wildcards.weight}.txt.gz | sed -e s/{wildcards.weight}.//g | gzip > resources/data/predicted_expression/{wildcards.weight}/Reference_Expression/Reference_Expression_{wildcards.weight}_mod.txt.gz; mv resources/data/predicted_expression/{wildcards.weight}/Reference_Expression/Reference_Expression_{wildcards.weight}_mod.txt.gz resources/data/predicted_expression/{wildcards.weight}/Reference_Expression/Reference_Expression_{wildcards.weight}.txt.gz"

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

##########
# Analyse GWAS summary statistics
##########

###
# Run TWAS
###

# run twas
rule run_twas:
  resources:
    mem_mb=20000
  input:
    "resources/data/gwas_sumstat/{gwas}/{gwas}.cleaned.munged.sumstats.gz",
    "resources/data/gwas_sumstat/{gwas}/{gwas}.cleaned.munged.median_N.txt",
    rules.install_fusion.output,
    rules.install_plink2R.output,
    rules.prep_1kg.output,
    rules.format_psychencode.output,
    rules.update_gtex_coord_all_panel.input,
    rules.insert_n_nongtex_all_panel.input,
  output:
    "results/{gwas}/twas/{weights}/{gwas}_twas_{weights}_chr{chr}"
  conda: 
    "../envs/GenoFunc.yaml"
  shell:
    "N=$(cat resources/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}.cleaned.munged.median_N.txt); Rscript resources/software/fusion/FUSION.assoc_test.R \
    --sumstats resources/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}.cleaned.munged.sumstats.gz \
    --weights resources/data/fusion_snp_weights/{wildcards.weights}/{wildcards.weights}.pos \
    --weights_dir resources/data/fusion_snp_weights/{wildcards.weights} \
    --ref_ld_chr resources/data/1kg/1KG.Phase3.EUR.MAF_001.chr \
    --out {output} \
    --chr {wildcards.chr} \
    --coloc_P 1e-3 \
    --GWASN ${{N}}"

rule twas_all_chr:
    input: 
      lambda w: expand("results/{gwas}/twas/{weight}/{gwas}_twas_{weight}_chr{chr}", gwas=w.gwas, weight=w.weight, chr=range(1, 23))
    output: 
      touch("results/{gwas}/checks/twas_{weight}_all_chr.done")

rule twas_all_panel:
    input: 
      lambda w: expand("results/{gwas}/checks/twas_{weight}_all_chr.done", gwas=w.gwas, weight=weights)
    output: 
      touch("results/{gwas}/checks/twas_all_panel.done")

# Combine TWAS results
checkpoint combine_twas_res:
  input:
    "results/{gwas}/checks/twas_all_panel.done"
  output:
    "results/{gwas}/twas/{gwas}_twas_GW_clean.txt.gz"
  conda: 
    "../envs/GenoFunc.yaml"
  shell:
    "Rscript scripts/combine_twas.R \
      --gwas {wildcards.gwas}"

# Identify chromosomes with significant associations
from pathlib import Path

def sig_chr_munge(x):
    checkpoint_output = checkpoints.combine_twas_res.get(gwas=x).output[0]
    checkpoint_output = "results/" + x + "/twas/" + x + "_twas_sig_chr.txt"
    sig_chr_df = pd.read_table(checkpoint_output, sep=' ')
    return sig_chr_df['x'].tolist()

def get_mem_mb_cond(wildcards, attempt):
    return attempt * 50000

# Run conditional analysis
rule run_conditional:
  resources: 
    mem_mb=get_mem_mb_cond 
  input:
    "results/{gwas}/twas/{gwas}_twas_GW_clean.txt.gz",
    rules.download_glist.output
  output:
    touch("results/{gwas}/checks/run_conditional_{gwas}_{chr}.done")
  conda: 
    "../envs/GenoFunc.yaml"
  shell:
    "mkdir -p results/{wildcards.gwas}/twas/conditional; Rscript resources/software/fusion/FUSION.post_process.R \
      --input results/{wildcards.gwas}/twas/{wildcards.gwas}_twas_GW_clean_sig.txt \
      --sumstats resources/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}.cleaned.munged.sumstats.gz \
      --report \
      --ref_ld_chr resources/data/1kg/1KG.Phase3.EUR.MAF_001.chr \
      --out results/{wildcards.gwas}/twas/conditional/{wildcards.gwas}_twas_conditional_chr{wildcards.chr} \
      --chr {wildcards.chr} \
      --save_loci \
      --ldsc F \
      --locus_win 500000"

rule conditional:
    input: 
      lambda w: expand("results/{gwas}/checks/run_conditional_{gwas}_{chr}.done", gwas=w.gwas, chr=sig_chr_munge("{}".format(w.gwas)))
    output: 
      touch("results/{gwas}/checks/conditional_all_chr.done")

# Process conditional results
rule process_conditional:
  input: 
    "results/{gwas}/checks/conditional_all_chr.done"
  output:
    "results/{gwas}/twas/{gwas}_twas_novelty.csv"
  conda:
    "../envs/GenoFunc.yaml"
  shell:
    "Rscript scripts/process_conditional.R \
      --gwas {wildcards.gwas}"

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
    --coloc_P 5e-2 \
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
    --coloc_P 5e-2 \
    --GWASN ${{N}}"

rule banner_pwas_all_chr:
    input: 
      lambda w: expand("results/{gwas}/pwas/banner/{gwas}_pwas_banner_chr{chr}", gwas=w.gwas, chr=range(1, 23))
    output: 
      touch("results/{gwas}/checks/banner_pwas_all_chr.done")

####
# Run So et al. analysis
####

# To enable the analysis to be run in parallel for each GWAS, I have split the So et al code into two parts.
# The LINCS data has been slit into 10 batches.
# Step 1 runs the pearson, spearman and KS tests on each batch of pertubagens seperately.
# Step 2 combines the results across all batches, estimates the average ranks, and calculates the p-value.

rule so_cmap_step1:
  input:
    "results/{gwas}/twas/{gwas}_twas_GW_clean.txt.gz",
    rules.install_cmapR.output,
    rules.subset_lincs.output
  output:
    "results/{gwas}/twas/cmap/res_{weight}_batch_{batch}.RDS"
  conda: 
    "../envs/GenoFunc.yaml"
  shell:
    "Rscript scripts/so_cmap_step1.R \
    --twas {wildcards.gwas} \
    --panel {wildcards.weight} \
    --batch {wildcards.batch}"

rule so_cmap_step1_all_batch:
    input: 
      lambda w: expand("results/{gwas}/twas/cmap/res_{weight}_batch_{batch}.RDS", gwas=w.gwas, batch=range(1, 11), weight=w.weight)
    output: 
      touch("results/{gwas}/checks/so_cmap_step1_{weight}_all_batch.done")

rule so_cmap_step2:
  input:
    "results/{gwas}/checks/so_cmap_step1_{weight}_all_batch.done"
  output:
    "results/{gwas}/twas/cmap/So_res_{weight}.csv"
  conda: 
    "../envs/GenoFunc.yaml"
  shell:
    "Rscript scripts/so_cmap_step2.R \
    --twas {wildcards.gwas} \
    --panel {wildcards.weight}"

# Format the so et al results 
rule format_so_results:
  input:
    "results/{gwas}/twas/cmap/So_res_{weight}.csv",
    rules.download_atc.output
  output:
    "results/{gwas}/twas/cmap/atc_res_{weight}.csv"
  conda: 
    "../envs/GenoFunc.yaml"
  params:
    lincs_siginfo= config["lincs_siginfo"]
  shell:
    "Rscript scripts/format_so_results.R \
    --twas {wildcards.gwas} \
    --panel {wildcards.weight} \
    --lincs_siginfo_path {params.lincs_siginfo}"

rule format_so_results_all_panel:
    input: 
      lambda w: expand("results/{gwas}/twas/cmap/atc_res_{weight}.csv", gwas=w.gwas, batch=range(1, 11), weight=weights_nosplice)
    output: 
      touch("results/{gwas}/checks/format_so_results_all_panel.done")


###
# TWAS-GSEA
###

# Run TWAS-GSEA
rule run_twas_gsea:
  resources: 
    mem_mb=50000,
    cpus=5
  input:
    rules.install_twas_gsea.output,
    "results/{gwas}/twas/{gwas}_twas_GW_clean.txt.gz",
    rules.format_lincs_for_twas_gsea.output,
    "resources/data/predicted_expression/format_pred_{weight}.done",
    rules.install_lme4qtl.output
  output:
    touch("results/{gwas}/twas/cmap/twas_gsea_cmap_{weight}.done")
  conda:
    "../envs/GenoFunc.yaml"
  shell:
    "Rscript resources/software/TWAS-GSEA/TWAS-GSEA.V1.2.R \
      --twas_results results/{wildcards.gwas}/twas/{wildcards.gwas}_twas_{wildcards.weight}_GW_clean.txt.gz \
      --pos resources/data/fusion_snp_weights/{wildcards.weight}/{wildcards.weight}.pos \
    	--prop_file resources/data/lincs/lincs_core_subset.txt.gz \
    	--expression_ref resources/data/predicted_expression/{wildcards.weight}/Reference_Expression/Reference_Expression_{wildcards.weight}.txt.gz \
    	--n_cores 5 \
    	--covar GeneLength,NSNP \
    	--use_alt_id ensembl_gene_id \
    	--linear_p_thresh 1 \
    	--min_Ngenes 2 \
    	--save_CorMat F \
    	--min_r2 0.01 \
    	--qqplot F \
    	--directional T \
    	--output results/{wildcards.gwas}/twas/cmap/twas_gsea_cmap_{wildcards.weight}"

# Format TWAS-GSEA results
rule format_twas_gsea_results:
  input:
    "results/{gwas}/twas/cmap/twas_gsea_cmap_{weight}.done",
    rules.download_atc.output
  output:
    "results/{gwas}/twas/cmap/twas_gsea_{weight}_res_atc_res.csv"
  conda: 
    "../envs/GenoFunc.yaml"
  params:
    lincs_siginfo= config["lincs_siginfo"]
  shell:
    "Rscript scripts/format_twas_gsea_results.R \
    --twas {wildcards.gwas} \
    --panel {wildcards.weight} \
    --lincs_siginfo_path {params.lincs_siginfo}"
    
rule format_twas_gsea_results_all_panel:
    input: 
      lambda w: expand("results/{gwas}/twas/cmap/twas_gsea_{weight}_res_atc_res.csv", gwas=w.gwas, batch=range(1, 11), weight=weights_nosplice)
    output: 
      touch("results/{gwas}/checks/format_twas_gsea_results_all_panel.done")
      
#######
# Run TWAS-GSEA using DrugTargetor sets
#######

# Format drugtargetor database
rule format_drug_targetor_for_twas_gsea:
  input:
    rules.download_drug_targetor.output,
    rules.download_magma_gene_loc.output
  output:
    "resources/data/drug_targetor/wholedatabase_for_targetor_directional.prop"
  conda: 
    "../envs/GenoFunc.yaml"
  shell:
    "Rscript scripts/format_drug_targetor_for_twas_gsea.R"

# Run TWAS-GSEA
rule run_twas_gsea_drug_targetor:
  resources: 
    mem_mb=50000,
    cpus=5
  input:
    rules.install_twas_gsea.output,
    "results/{gwas}/twas/{gwas}_twas_GW_clean.txt.gz",
    rules.format_drug_targetor_for_twas_gsea.output,
    "resources/data/predicted_expression/format_pred_{weight}.done",
    rules.install_lme4qtl.output
  output:
    touch("results/{gwas}/twas/drugtargetor/twas_gsea_drugtargetor_{weight}.done")
  conda:
    "../envs/GenoFunc.yaml"
  shell:
    "Rscript resources/software/TWAS-GSEA/TWAS-GSEA.V1.2.R \
      --twas_results results/{wildcards.gwas}/twas/{wildcards.gwas}_twas_{wildcards.weight}_GW_clean.txt.gz \
      --pos resources/data/fusion_snp_weights/{wildcards.weight}/{wildcards.weight}.pos \
    	--prop_file resources/data/drug_targetor/wholedatabase_for_targetor_directional.prop \
    	--expression_ref resources/data/predicted_expression/{wildcards.weight}/Reference_Expression/Reference_Expression_{wildcards.weight}.txt.gz \
    	--n_cores 5 \
    	--covar GeneLength,NSNP \
    	--use_alt_id ID \
    	--linear_p_thresh 1 \
    	--min_Ngenes 2 \
    	--save_CorMat F \
    	--min_r2 0.01 \
    	--qqplot F \
    	--directional T \
    	--output results/{wildcards.gwas}/twas/drugtargetor/twas_gsea_drugtargetor_{wildcards.weight}"
    	
# Format the output
rule format_twas_gsea_drugtargetor_results:
  input:
    "results/{gwas}/twas/drugtargetor/twas_gsea_drugtargetor_{weight}.done",
    rules.download_atc.output
  output:
    "results/{gwas}/twas/drugtargetor/twas_gsea_{weight}_res_atc_res.csv"
  conda: 
    "../envs/GenoFunc.yaml"
  shell:
    "Rscript scripts/format_twas_gsea_drugtargetor_results.R \
    --twas {wildcards.gwas} \
    --panel {wildcards.weight}"
    
rule format_twas_gsea_drugtargetor_results_all_panel:
    input: 
      lambda w: expand("results/{gwas}/twas/drugtargetor/twas_gsea_{weight}_res_atc_res.csv", gwas=w.gwas, batch=range(1, 11), weight=weights_nosplice)
    output: 
      touch("results/{gwas}/checks/format_twas_gsea_drugtargetor_results_all_panel.done")

