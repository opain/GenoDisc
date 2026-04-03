####
# Download software required for TWAS-related analysis
####

# Install fusion
rule install_fusion:
  output:
    directory(f"{resdir}/software/fusion/")
  conda:
    "../envs/main.yaml"
  shell:
    "git clone https://github.com/gusevlab/fusion_twas.git {output}; \
    cd {output}; \
    git reset --hard e1ba5f7f3907e6f586f7fb5bb115b35cc0d3c0c2"

# Download plink2R
rule download_plink2R:
  input:
    rules.install_fusion.output
  output:
    f"{resdir}/software/plink2R/plink2R-master/data.bed"
  conda:
    "../envs/main.yaml"
  params:
    resdir=resdir
  shell:
    "rm -rf {params.resdir}/software/plink2R; \
    mkdir -p {params.resdir}/software/plink2R; \
    wget -O {params.resdir}/software/plink2R/master.zip https://github.com/gabraham/plink2R/archive/master.zip; \
    unzip {params.resdir}/software/plink2R/master.zip -d {params.resdir}/software/plink2R"

# Install plink2R
rule install_plink2R:
  input:
    rules.download_plink2R.output,
    "envs/main.yaml"
  output:
    touch(f"{resdir}/software/install_plink2R")
  conda:
    "../envs/main.yaml"
  params:
    resdir=resdir
  shell:
    "Rscript -e 'install.packages(\"{params.resdir}/software/plink2R/plink2R-master/plink2R/\",repos=NULL)'"

# Install SNP-weights pipeline repo
rule install_snp_weight_pipe:
  output:
    directory(f"{resdir}/software/Calculating-FUSION-TWAS-weights-pipeline/")
  conda:
    "../envs/main.yaml"
  shell:
    "git clone https://github.com/opain/Calculating-FUSION-TWAS-weights-pipeline.git {output}; \
    cd {output}; \
    git reset --hard ab15a41e4568107f29bc5a538ea016a554d58589"

####
# Dowload data for TWAS related analysis
####

# We want to download weights and harmonise into a single repository
# First download the weights, then reformat accordingly

# Download PsychENCODE SNP-weights
rule download_psychENCODE_weights:
  output:
    touch(f"{resdir}/data/download_psychENCODE_weights.done")
  conda:
    "../envs/main.yaml"
  params:
    resdir=resdir
  shell:
    "mkdir -p {params.resdir}/data/fusion_snp_weights/psychencode; \
    wget -O {params.resdir}/data/fusion_snp_weights/psychencode/PEC_TWAS_weights.tar.gz http://resource.psychencode.org/Datasets/Derived/PEC_TWAS_weights.tar.gz; \
    mkdir -p {params.resdir}/data/fusion_snp_weights/psychencode/psychencode; \
    tar xvzf {params.resdir}/data/fusion_snp_weights/psychencode/PEC_TWAS_weights.tar.gz -C {params.resdir}/data/fusion_snp_weights/psychencode/psychencode; \
    rm {params.resdir}/data/fusion_snp_weights/psychencode/PEC_TWAS_weights.tar.gz"

# Format PsychENCODE SNP-weights
rule format_psychencode:
  input:
    psychencode_data=rules.download_psychENCODE_weights.output,
    weights_pipe=rules.install_snp_weight_pipe.output,
    biomart=rules.download_biomart.output
  output:
    f"{resdir}/data/format_psychencode.done"
  conda:
    "../envs/main.yaml"
  shell:
    "Rscript scripts/format_psychENCODE.R"

# Download FUSION GTEx v8 EUR SNP-weights
# I am using EUR instead of full sample to avoid LD mismatch
gtex_weights=config["gtex_weights"]

rule download_gtex_weights:
  output:
    touch(f"{resdir}/data/download_fusion_gtex_{{weight}}_weights.done")
  conda:
    "../envs/main.yaml"
  params:
    resdir=resdir
  shell:
    "mkdir -p {params.resdir}/data/fusion_snp_weights/{wildcards.weight}; wget -O {params.resdir}/data/fusion_snp_weights/GTExv8.EUR.{wildcards.weight}.tar.gz https://s3.us-west-1.amazonaws.com/gtex.v8.fusion/EUR/GTExv8.EUR.{wildcards.weight}.tar.gz; tar xf {params.resdir}/data/fusion_snp_weights/GTExv8.EUR.{wildcards.weight}.tar.gz -C {params.resdir}/data/fusion_snp_weights/{wildcards.weight}; rm {params.resdir}/data/fusion_snp_weights/GTExv8.EUR.{wildcards.weight}.tar.gz; mv {params.resdir}/data/fusion_snp_weights/{wildcards.weight}/GTExv8.EUR.{wildcards.weight}.pos {params.resdir}/data/fusion_snp_weights/{wildcards.weight}/{wildcards.weight}.pos; mv {params.resdir}/data/fusion_snp_weights/{wildcards.weight}/GTExv8.EUR.{wildcards.weight} {params.resdir}/data/fusion_snp_weights/{wildcards.weight}/{wildcards.weight}"

# Update GTEx v8 P0 and P1 to build GRCh 37
rule update_gtex_coord:
  input:
    f"{resdir}/data/download_fusion_gtex_{{weight}}_weights.done",
    rules.download_biomart.output
  output:
    touch(f"{resdir}/data/update_gtex_coord_{{weight}}.done")
  conda:
    "../envs/main.yaml"
  shell:
    "Rscript scripts/update_gtex_coord.R \
      --panel {wildcards.weight}"

rule update_gtex_coord_all_panel:
    input: expand(f"{resdir}/data/update_gtex_coord_{{weight}}.done", weight=gtex_weights)

# Download FUSION non-GTEx SNP-weights
non_gtex_weights=config["non_gtex_weights"]

rule download_non_gtex_weights:
  output:
    touch(f"{resdir}/data/download_non_gtex_{{weight}}_weights.done")
  conda:
    "../envs/main.yaml"
  params:
    resdir=resdir
  shell:
    "mkdir -p {params.resdir}/data/fusion_snp_weights/{wildcards.weight}; wget --no-check-certificate -O {params.resdir}/data/fusion_snp_weights/{wildcards.weight}.tar.bz2 https://data.broadinstitute.org/alkesgroup/FUSION/WGT/{wildcards.weight}.tar.bz2; tar xvjf {params.resdir}/data/fusion_snp_weights/{wildcards.weight}.tar.bz2 -C {params.resdir}/data/fusion_snp_weights/{wildcards.weight}; rm {params.resdir}/data/fusion_snp_weights/{wildcards.weight}.tar.bz2"

# Insert N into non-GTEX SNP-weights
rule insert_n_nongtex:
  input:
    f"{resdir}/data/download_non_gtex_{{weight}}_weights.done",
    rules.download_biomart.output
  output:
    touch(f"{resdir}/data/insert_n_nongtex_{{weight}}.done")
  conda:
    "../envs/main.yaml"
  shell:
    "Rscript scripts/insert_n_nongtex.R \
      --panel {wildcards.weight}"

rule insert_n_nongtex_all_panel:
    input: expand(f"{resdir}/data/insert_n_nongtex_{{weight}}.done", weight=non_gtex_weights)

# Download glist file
rule download_glist:
  output:
    f"{resdir}/data/glist-hg19"
  conda:
    "../envs/main.yaml"
  params:
    resdir=resdir
  shell:
    "wget -P {params.resdir}/data/ https://www.cog-genomics.org/static/bin/plink/glist-hg19"

####
# Download TWAS-GSEA
####

rule install_twas_gsea:
  output:
    directory(f"{resdir}/software/TWAS-GSEA/")
  conda:
    "../envs/main.yaml"
  shell:
    "git clone https://github.com/opain/TWAS-GSEA.git {output}; \
    cd {output}; \
    git reset --hard d9b98a670121bcf686448b5d65c8d3bc443ba494"

####
# Download FeaturePred
####

rule install_feature_pred:
  output:
    directory(f"{resdir}/software/Predicting-TWAS-features/")
  conda:
    "../envs/main.yaml"
  shell:
    "git clone https://github.com/opain/Predicting-TWAS-features.git {output}; \
    cd {output}; \
    git reset --hard b9defcf3c96145ab86f605629c48e0d29daebe0c"

####
# Download pigz
####

rule install_pigz:
  output:
    f"{resdir}/software/pigz/pigz/pigz"
  conda:
    "../envs/main.yaml"
  params:
    resdir=resdir
  shell:
    "wget -O {params.resdir}/software/pigz.tar.gz https://zlib.net/pigz/pigz.tar.gz; mkdir -p {params.resdir}/software/pigz; tar xvzf {params.resdir}/software/pigz.tar.gz -C {params.resdir}/software/pigz; rm {params.resdir}/software/pigz.tar.gz; cd {params.resdir}/software/pigz/pigz; make"

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
  conda:
    "../envs/main.yaml"
  params:
    resdir=resdir
  shell:
    "Rscript {params.resdir}/software/Predicting-TWAS-features/FeaturePred.V2.0.R \
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
    	--output {params.resdir}/data/predicted_expression/{wildcards.weight}"

# Format expression data for TWAS-GSEA (i.e. remove PANEL from column names)
rule format_pred:
  input:
    f"{resdir}/data/predicted_expression/{{weight}}/Reference_Expression/Reference_Expression_{{weight}}.txt.gz"
  output:
    touch(f"{resdir}/data/predicted_expression/format_pred_{{weight}}.done")
  conda:
    "../envs/main.yaml"
  params:
    resdir=resdir
  shell:
    "zcat {params.resdir}/data/predicted_expression/{wildcards.weight}/Reference_Expression/Reference_Expression_{wildcards.weight}.txt.gz | sed -e s/{wildcards.weight}.//g | gzip > {params.resdir}/data/predicted_expression/{wildcards.weight}/Reference_Expression/Reference_Expression_{wildcards.weight}_mod.txt.gz; mv {params.resdir}/data/predicted_expression/{wildcards.weight}/Reference_Expression/Reference_Expression_{wildcards.weight}_mod.txt.gz {params.resdir}/data/predicted_expression/{wildcards.weight}/Reference_Expression/Reference_Expression_{wildcards.weight}.txt.gz"

####
# Install lme4qtl
####
# Note the version of conda was not working in R 4.0.2

rule install_lme4qtl:
  input:
    "envs/main.yaml"
  output:
    touch(f"{resdir}/software/install_lme4qtl.done")
  conda:
    "../envs/main.yaml"
  shell:
    "Rscript -e 'devtools::install_github(\"variani/lme4qtl\", ref = \"0.1.10\")'"

####
# Format ROSMAP and Banner PWAS data
####

rule format_pwas_data:
  output:
    f"{resdir}/data/banner_twas/Banner.n152.fusion.WEIGHTS/train_weights_withN.pos"
  conda:
    "../envs/main.yaml"
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
    "{outdir}/data/gwas_sumstat/{gwas}/{gwas}.cleaned.munged.sumstats.gz",
    "{outdir}/data/gwas_sumstat/{gwas}/{gwas}.cleaned.munged.median_N.txt",
    rules.install_fusion.output,
    rules.install_plink2R.output,
    rules.prep_1kg.output,
    rules.format_psychencode.output,
    rules.update_gtex_coord_all_panel.input,
    rules.insert_n_nongtex_all_panel.input,
  output:
    "{outdir}/results/{gwas}/twas/{weights}/{gwas}_twas_{weights}_chr{chr}"
  conda:
    "../envs/main.yaml"
  params:
    resdir=resdir
  shell:
    "N=$(cat {outdir}/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}.cleaned.munged.median_N.txt); Rscript {params.resdir}/software/fusion/FUSION.assoc_test.R \
    --sumstats {outdir}/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}.cleaned.munged.sumstats.gz \
    --weights {params.resdir}/data/fusion_snp_weights/{wildcards.weights}/{wildcards.weights}.pos \
    --weights_dir {params.resdir}/data/fusion_snp_weights/{wildcards.weights} \
    --ref_ld_chr {params.resdir}/data/1kg/1KG.Phase3.EUR.MAF_001.chr \
    --out {output} \
    --chr {wildcards.chr} \
    --coloc_P 1e-3 \
    --GWASN ${{N}}"

rule twas_all_chr:
    input:
      lambda w: expand("{outdir}/results/{gwas}/twas/{weight}/{gwas}_twas_{weight}_chr{chr}", gwas=w.gwas, weight=w.weight, chr=chromosomes, outdir={outdir})
    output:
      touch("{outdir}/results/{gwas}/checks/twas_{weight}_all_chr.done")

rule twas_all_panel:
    input:
      lambda w: expand("{outdir}/results/{gwas}/checks/twas_{weight}_all_chr.done", gwas=w.gwas, weight=weights, outdir={outdir})
    output:
      touch("{outdir}/results/{gwas}/checks/twas_all_panel.done")

# Combine TWAS results
# Delete conditional results folder to avoid conflicts during reruns
checkpoint combine_twas_res:
  input:
    "{outdir}/results/{gwas}/checks/twas_all_panel.done",
    rules.download_biomart.output
  output:
    "{outdir}/results/{gwas}/twas/{gwas}_twas_GW_clean.txt.gz"
  conda:
    "../envs/main.yaml"
  params:
    config_file= config["config_file"]
  shell:
    "Rscript scripts/combine_twas.R \
      --gwas {wildcards.gwas} \
      --config_file {params.config_file}; \
      rm -rf {outdir}/results/{wildcards.gwas}/twas/conditional"

# Identify chromosomes with significant associations
from pathlib import Path

def sig_chr_munge(x):
    checkpoint_output = checkpoints.combine_twas_res.get(gwas=x, outdir=outdir).output[0]
    checkpoint_output = outdir + "/results/" + x + "/twas/" + x + "_twas_sig_chr.txt"
    sig_chr_df = pd.read_table(checkpoint_output, sep=' ')
    return sig_chr_df['x'].tolist()

def get_mem_mb_cond(wildcards, attempt):
    return attempt * 50000

# Run conditional analysis
rule run_conditional:
  resources:
    mem_mb=get_mem_mb_cond
  input:
    "{outdir}/results/{gwas}/twas/{gwas}_twas_GW_clean.txt.gz",
    rules.download_glist.output
  output:
    touch("{outdir}/results/{gwas}/checks/run_conditional_{gwas}_{chr}.done")
  conda:
    "../envs/main.yaml"
  params:
    resdir=resdir
  shell:
    "mkdir -p {outdir}/results/{wildcards.gwas}/twas/conditional; Rscript {params.resdir}/software/fusion/FUSION.post_process.R \
      --input {outdir}/results/{wildcards.gwas}/twas/{wildcards.gwas}_twas_GW_clean_sig.txt \
      --sumstats {outdir}/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}.cleaned.munged.sumstats.gz \
      --report \
      --ref_ld_chr {params.resdir}/data/1kg/1KG.Phase3.EUR.MAF_001.chr \
      --out {outdir}/results/{wildcards.gwas}/twas/conditional/{wildcards.gwas}_twas_conditional_chr{wildcards.chr} \
      --chr {wildcards.chr} \
      --save_loci \
      --ldsc F \
      --locus_win 500000"

rule conditional:
    input:
      lambda w: expand("{outdir}/results/{gwas}/checks/run_conditional_{gwas}_{chr}.done", gwas=w.gwas, chr=sig_chr_munge("{}".format(w.gwas)), outdir={outdir})
    output:
      touch("{outdir}/results/{gwas}/checks/conditional_all_chr.done")

# Process conditional results
rule process_conditional:
  input:
    "{outdir}/results/{gwas}/checks/conditional_all_chr.done"
  output:
    "{outdir}/results/{gwas}/twas/{gwas}_twas_novelty.csv"
  conda:
    "../envs/main.yaml"
  params:
    config_file=config['config_file']
  shell:
    "Rscript scripts/process_conditional.R \
      --gwas {wildcards.gwas} \
      --config_file {params.config_file}"

###
# Run PWAS
###

# Run twas using ROSMAP SNP-weights
rule run_rosmap_pwas:
  resources: mem_mb=20000
  input:
    sumstats="{outdir}/data/gwas_sumstat/{gwas}/{gwas}.cleaned.munged.sumstats.gz",
    neff_txt="{outdir}/data/gwas_sumstat/{gwas}/{gwas}.cleaned.munged.median_N.txt",
    fusion=rules.install_fusion.output,
    plink2R=rules.install_plink2R.output,
    format_psychencode=rules.format_pwas_data.output,
    prep_1kg=rules.prep_1kg.output
  output:
    "{outdir}/results/{gwas}/pwas/rosmap/{gwas}_pwas_rosmap_chr{chr}"
  conda:
    "../envs/main.yaml"
  params:
    resdir=resdir
  shell:
    "N=$(cat {outdir}/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}.cleaned.munged.median_N.txt); Rscript {params.resdir}/software/fusion/FUSION.assoc_test.R \
    --sumstats {input.sumstats} \
    --weights {params.resdir}/data/rosmap_twas/ROSMAP.n376.fusion.WEIGHTS/train_weights_withN.pos \
    --weights_dir {params.resdir}/data/rosmap_twas/ROSMAP.n376.fusion.WEIGHTS \
    --ref_ld_chr {params.resdir}/data/1kg/1KG.Phase3.EUR.MAF_001.chr \
    --out {output} \
    --chr {wildcards.chr} \
    --coloc_P 5e-2 \
    --GWASN ${{N}}"

rule rosmap_pwas_all_chr:
    input:
      lambda w: expand("{outdir}/results/{gwas}/pwas/rosmap/{gwas}_pwas_rosmap_chr{chr}", gwas=w.gwas, chr=chromosomes, outdir={outdir})
    output:
      touch("{outdir}/results/{gwas}/checks/rosmap_pwas_all_chr.done")

# Run twas using Banner SNP-weights
rule run_banner_pwas:
  resources: mem_mb=20000
  input:
    sumstats="{outdir}/data/gwas_sumstat/{gwas}/{gwas}.cleaned.munged.sumstats.gz",
    neff_txt="{outdir}/data/gwas_sumstat/{gwas}/{gwas}.cleaned.munged.median_N.txt",
    fusion=rules.install_fusion.output,
    plink2R=rules.install_plink2R.output,
    format_psychencode=rules.format_pwas_data.output,
    prep_1kg=rules.prep_1kg.output
  output:
    "{outdir}/results/{gwas}/pwas/banner/{gwas}_pwas_banner_chr{chr}"
  conda:
    "../envs/main.yaml"
  params:
    resdir=resdir
  shell:
    "N=$(cat {outdir}/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}.cleaned.munged.median_N.txt); Rscript {params.resdir}/software/fusion/FUSION.assoc_test.R \
    --sumstats {input.sumstats} \
    --weights {params.resdir}/data/banner_twas/Banner.n152.fusion.WEIGHTS/train_weights_withN.pos \
    --weights_dir {params.resdir}/data/banner_twas/Banner.n152.fusion.WEIGHTS \
    --ref_ld_chr {params.resdir}/data/1kg/1KG.Phase3.EUR.MAF_001.chr \
    --out {output} \
    --chr {wildcards.chr} \
    --coloc_P 5e-2 \
    --GWASN ${{N}}"

rule banner_pwas_all_chr:
    input:
      lambda w: expand("{outdir}/results/{gwas}/pwas/banner/{gwas}_pwas_banner_chr{chr}", gwas=w.gwas, chr=chromosomes, outdir={outdir})
    output:
      touch("{outdir}/results/{gwas}/checks/banner_pwas_all_chr.done")

#######
# Run TWAS-GSEA using DrugTargetor sets
#######

# Format drugtargetor database
rule format_drug_targetor_for_twas_gsea:
  input:
    rules.download_drug_targetor.output,
    rules.download_magma_gene_loc.output
  output:
    f"{resdir}/data/drug_targetor/wholedatabase_for_targetor_directional.prop"
  conda:
    "../envs/main.yaml"
  shell:
    "Rscript scripts/format_drug_targetor_for_twas_gsea.R"

# Run TWAS-GSEA
rule run_twas_gsea_drug_targetor:
  resources:
    mem_mb=50000,
    cpus=5
  input:
    rules.install_twas_gsea.output,
    "{outdir}/results/{gwas}/twas/{gwas}_twas_GW_clean.txt.gz",
    rules.format_drug_targetor_for_twas_gsea.output,
    f"{resdir}/data/predicted_expression/format_pred_{{weight}}.done",
    rules.install_lme4qtl.output
  output:
    touch("{outdir}/results/{gwas}/twas/drugtargetor/twas_gsea_drugtargetor_{weight}.done")
  conda:
    "../envs/main.yaml"
  params:
    resdir=resdir
  shell:
    "Rscript {params.resdir}/software/TWAS-GSEA/TWAS-GSEA.V1.2.R \
      --twas_results {outdir}/results/{wildcards.gwas}/twas/{wildcards.gwas}_twas_{wildcards.weight}_GW_clean.txt.gz \
      --pos {params.resdir}/data/fusion_snp_weights/{wildcards.weight}/{wildcards.weight}.pos \
    	--prop_file {params.resdir}/data/drug_targetor/wholedatabase_for_targetor_directional.prop \
    	--expression_ref {params.resdir}/data/predicted_expression/{wildcards.weight}/Reference_Expression/Reference_Expression_{wildcards.weight}.txt.gz \
    	--n_cores 5 \
    	--covar GeneLength,NSNP \
    	--use_alt_id ID \
    	--linear_p_thresh 1 \
    	--min_Ngenes 2 \
    	--save_CorMat F \
    	--min_r2 0.01 \
    	--qqplot F \
    	--directional T \
    	--output {outdir}/results/{wildcards.gwas}/twas/drugtargetor/twas_gsea_drugtargetor_{wildcards.weight}"

# Format the output
rule format_twas_gsea_drugtargetor_results:
  input:
    "{outdir}/results/{gwas}/twas/drugtargetor/twas_gsea_drugtargetor_{weight}.done",
    rules.download_atc.output
  output:
    "{outdir}/results/{gwas}/twas/drugtargetor/twas_gsea_{weight}_res_atc_res.csv"
  conda:
    "../envs/main.yaml"
  params:
    config_file=config['config_file']
  shell:
    "Rscript scripts/format_twas_gsea_drugtargetor_results.R \
    --twas {wildcards.gwas} \
    --panel {wildcards.weight} \
    --config_file {params.config_file}"

rule format_twas_gsea_drugtargetor_results_all_panel:
    input:
      lambda w: expand("{outdir}/results/{gwas}/twas/drugtargetor/twas_gsea_{weight}_res_atc_res.csv", gwas=w.gwas, batch=range(1, 11), weight=weights_nosplice, outdir={outdir})
    output:
      touch("{outdir}/results/{gwas}/checks/format_twas_gsea_drugtargetor_results_all_panel.done")

