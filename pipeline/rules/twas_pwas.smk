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
    "{outdir}/results/{gwas}/gwas_sumstat/{gwas}.cleaned.munged.sumstats.gz",
    "{outdir}/results/{gwas}/gwas_sumstat/{gwas}.cleaned.munged.median_N.txt",
    rules.install_fusion.output,
    rules.install_plink2R.output,
    rules.prep_1kg.output,
    rules.format_psychencode.output,
    rules.update_gtex_coord_all_panel.input,
    rules.insert_n_nongtex_all_panel.input,
  output:
    "{outdir}/results/{gwas}/twas/{weights}/{gwas}_twas_{weights}_chr{chr}"
  benchmark:
    "{outdir}/benchmarks/run_twas_{gwas}_{weights}_chr{chr}.tsv"
  conda:
    "../envs/main.yaml"
  params:
    resdir=resdir
  log:
    "{outdir}/logs/run_twas-{gwas}-{weights}-chr{chr}.log"
  shell:
    "(N=$(cat {outdir}/results/{wildcards.gwas}/gwas_sumstat/{wildcards.gwas}.cleaned.munged.median_N.txt); Rscript --vanilla {params.resdir}/software/fusion/FUSION.assoc_test.R \
    --sumstats {outdir}/results/{wildcards.gwas}/gwas_sumstat/{wildcards.gwas}.cleaned.munged.sumstats.gz \
    --weights {params.resdir}/data/fusion_snp_weights/{wildcards.weights}/{wildcards.weights}.pos \
    --weights_dir {params.resdir}/data/fusion_snp_weights/{wildcards.weights} \
    --ref_ld_chr {params.resdir}/data/1kg/1KG.Phase3.EUR.MAF_001.chr \
    --out {output} \
    --chr {wildcards.chr} \
    --coloc_P 1e-3 \
    --GWASN ${{N}}) > {log} 2>&1"

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
    "{outdir}/results/{gwas}/twas/{gwas}_twas_GW_clean.txt.gz",
    expand("{{outdir}}/results/{{gwas}}/twas/{{gwas}}_twas_{weight}_GW_clean.txt.gz",
           weight=weights_nosplice)
  benchmark:
    "{outdir}/benchmarks/combine_twas_res_{gwas}.tsv"
  conda:
    "../envs/main.yaml"
  params:
    config_file= config["config_file"]
  log:
    "{outdir}/logs/combine_twas_res-{gwas}.log"
  shell:
    "(Rscript --vanilla scripts/combine_twas.R \
      --gwas {wildcards.gwas} \
      --config_file {params.config_file}; \
      rm -rf {outdir}/results/{wildcards.gwas}/twas/conditional) > {log} 2>&1"

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
  benchmark:
    "{outdir}/benchmarks/run_conditional_{gwas}_chr{chr}.tsv"
  conda:
    "../envs/main.yaml"
  params:
    resdir=resdir
  log:
    "{outdir}/logs/run_conditional-{gwas}-chr{chr}.log"
  shell:
    "(mkdir -p {outdir}/results/{wildcards.gwas}/twas/conditional; Rscript --vanilla {params.resdir}/software/fusion/FUSION.post_process.R \
      --input {outdir}/results/{wildcards.gwas}/twas/{wildcards.gwas}_twas_GW_clean_sig.txt \
      --sumstats {outdir}/results/{wildcards.gwas}/gwas_sumstat/{wildcards.gwas}.cleaned.munged.sumstats.gz \
      --report \
      --ref_ld_chr {params.resdir}/data/1kg/1KG.Phase3.EUR.MAF_001.chr \
      --out {outdir}/results/{wildcards.gwas}/twas/conditional/{wildcards.gwas}_twas_conditional_chr{wildcards.chr} \
      --chr {wildcards.chr} \
      --save_loci \
      --ldsc F \
      --locus_win 500000) > {log} 2>&1"

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
  benchmark:
    "{outdir}/benchmarks/process_conditional_{gwas}.tsv"
  conda:
    "../envs/main.yaml"
  params:
    config_file=config['config_file']
  log:
    "{outdir}/logs/process_conditional-{gwas}.log"
  shell:
    "Rscript --vanilla scripts/process_conditional.R \
      --gwas {wildcards.gwas} \
      --config_file {params.config_file} > {log} 2>&1"

###
# Run PWAS
###

# Run twas using ROSMAP SNP-weights
rule run_rosmap_pwas:
  resources: mem_mb=20000
  input:
    sumstats="{outdir}/results/{gwas}/gwas_sumstat/{gwas}.cleaned.munged.sumstats.gz",
    neff_txt="{outdir}/results/{gwas}/gwas_sumstat/{gwas}.cleaned.munged.median_N.txt",
    fusion=rules.install_fusion.output,
    plink2R=rules.install_plink2R.output,
    format_psychencode=rules.format_pwas_data.output,
    prep_1kg=rules.prep_1kg.output
  output:
    "{outdir}/results/{gwas}/pwas/rosmap/{gwas}_pwas_rosmap_chr{chr}"
  benchmark:
    "{outdir}/benchmarks/run_rosmap_pwas_{gwas}_chr{chr}.tsv"
  conda:
    "../envs/main.yaml"
  params:
    resdir=resdir
  log:
    "{outdir}/logs/run_rosmap_pwas-{gwas}-chr{chr}.log"
  shell:
    "(N=$(cat {outdir}/results/{wildcards.gwas}/gwas_sumstat/{wildcards.gwas}.cleaned.munged.median_N.txt); Rscript --vanilla {params.resdir}/software/fusion/FUSION.assoc_test.R \
    --sumstats {input.sumstats} \
    --weights {params.resdir}/data/rosmap_twas/ROSMAP.n376.fusion.WEIGHTS/train_weights_withN.pos \
    --weights_dir {params.resdir}/data/rosmap_twas/ROSMAP.n376.fusion.WEIGHTS \
    --ref_ld_chr {params.resdir}/data/1kg/1KG.Phase3.EUR.MAF_001.chr \
    --out {output} \
    --chr {wildcards.chr} \
    --coloc_P 5e-2 \
    --GWASN ${{N}}) > {log} 2>&1"

rule rosmap_pwas_all_chr:
    input:
      lambda w: expand("{outdir}/results/{gwas}/pwas/rosmap/{gwas}_pwas_rosmap_chr{chr}", gwas=w.gwas, chr=chromosomes, outdir={outdir})
    output:
      touch("{outdir}/results/{gwas}/checks/rosmap_pwas_all_chr.done")

# Run twas using Banner SNP-weights
rule run_banner_pwas:
  resources: mem_mb=20000
  input:
    sumstats="{outdir}/results/{gwas}/gwas_sumstat/{gwas}.cleaned.munged.sumstats.gz",
    neff_txt="{outdir}/results/{gwas}/gwas_sumstat/{gwas}.cleaned.munged.median_N.txt",
    fusion=rules.install_fusion.output,
    plink2R=rules.install_plink2R.output,
    format_psychencode=rules.format_pwas_data.output,
    prep_1kg=rules.prep_1kg.output
  output:
    "{outdir}/results/{gwas}/pwas/banner/{gwas}_pwas_banner_chr{chr}"
  benchmark:
    "{outdir}/benchmarks/run_banner_pwas_{gwas}_chr{chr}.tsv"
  conda:
    "../envs/main.yaml"
  params:
    resdir=resdir
  log:
    "{outdir}/logs/run_banner_pwas-{gwas}-chr{chr}.log"
  shell:
    "(N=$(cat {outdir}/results/{wildcards.gwas}/gwas_sumstat/{wildcards.gwas}.cleaned.munged.median_N.txt); Rscript --vanilla {params.resdir}/software/fusion/FUSION.assoc_test.R \
    --sumstats {input.sumstats} \
    --weights {params.resdir}/data/banner_twas/Banner.n152.fusion.WEIGHTS/train_weights_withN.pos \
    --weights_dir {params.resdir}/data/banner_twas/Banner.n152.fusion.WEIGHTS \
    --ref_ld_chr {params.resdir}/data/1kg/1KG.Phase3.EUR.MAF_001.chr \
    --out {output} \
    --chr {wildcards.chr} \
    --coloc_P 5e-2 \
    --GWASN ${{N}}) > {log} 2>&1"

rule banner_pwas_all_chr:
    input:
      lambda w: expand("{outdir}/results/{gwas}/pwas/banner/{gwas}_pwas_banner_chr{chr}", gwas=w.gwas, chr=chromosomes, outdir={outdir})
    output:
      touch("{outdir}/results/{gwas}/checks/banner_pwas_all_chr.done")

#######
# Run TWAS-GSEA using DrugTargetor sets
#######

# Pre-compute the gene-gene predicted-expression correlation matrix once per
# weight panel. The result is reused across every (gwas, gene-set) call to
# TWAS-GSEA-fast.R for that panel, so this rule has no {gwas} wildcard.
rule build_twas_gsea_cormat:
  resources:
    mem_mb=50000,
    cpus=5
  input:
    rules.install_twas_gsea.output,
    f"{resdir}/data/fusion_snp_weights/{{weight}}/{{weight}}.pos",
    f"{resdir}/data/predicted_expression/format_pred_{{weight}}.done"
  output:
    f"{resdir}/data/predicted_expression/{{weight}}/Reference_Expression/{{weight}}.CorMat.RDS"
  benchmark:
    f"{resdir}/benchmarks/build_twas_gsea_cormat_{{weight}}.tsv"
  conda:
    "../envs/main.yaml"
  params:
    resdir=resdir
  log:
    f"{resdir}/logs/build_twas_gsea_cormat-{{weight}}.log"
  shell:
    "Rscript --vanilla {params.resdir}/software/TWAS-GSEA/build_cor_matrix.R \
      --expression_ref {params.resdir}/data/predicted_expression/{wildcards.weight}/Reference_Expression/Reference_Expression_{wildcards.weight}.txt.gz \
      --pos {params.resdir}/data/fusion_snp_weights/{wildcards.weight}/{wildcards.weight}.pos \
      --min_r2 0.01 \
      --n_cores 5 \
      --output {params.resdir}/data/predicted_expression/{wildcards.weight}/Reference_Expression/{wildcards.weight} > {log} 2>&1"

# Run TWAS-GSEA-fast.R against the precomputed cor matrix.
rule run_twas_gsea_drug_targetor:
  wildcard_constraints:
    weight="(?!nondir_).+"
  resources:
    mem_mb=50000,
    cpus=5
  input:
    rules.install_twas_gsea.output,
    "{outdir}/results/{gwas}/twas/{gwas}_twas_{weight}_GW_clean.txt.gz",
    rules.format_drug_targetor_for_twas_gsea.output,
    f"{resdir}/data/predicted_expression/{{weight}}/Reference_Expression/{{weight}}.CorMat.RDS"
  output:
    touch("{outdir}/results/{gwas}/twas/drugtargetor/twas_gsea_drugtargetor_{weight}.done")
  benchmark:
    "{outdir}/benchmarks/run_twas_gsea_drug_targetor_{gwas}_{weight}.tsv"
  conda:
    "../envs/main.yaml"
  params:
    resdir=resdir
  log:
    "{outdir}/logs/run_twas_gsea_drug_targetor-{gwas}-{weight}.log"
  shell:
    "Rscript --vanilla {params.resdir}/software/TWAS-GSEA/TWAS-GSEA-fast.R \
      --twas_results {outdir}/results/{wildcards.gwas}/twas/{wildcards.gwas}_twas_{wildcards.weight}_GW_clean.txt.gz \
      --pos {params.resdir}/data/fusion_snp_weights/{wildcards.weight}/{wildcards.weight}.pos \
      --input_CorMat {params.resdir}/data/predicted_expression/{wildcards.weight}/Reference_Expression/{wildcards.weight}.CorMat.RDS \
      --prop_file {params.resdir}/data/drug_targetor/wholedatabase_for_targetor_directional.prop \
      --n_cores 5 \
      --covar GeneLength,NSNP \
      --use_alt_id ID \
      --min_Ngenes 2 \
      --directional T \
      --output {outdir}/results/{wildcards.gwas}/twas/drugtargetor/twas_gsea_drugtargetor_{wildcards.weight} > {log} 2>&1"

# Format the output
rule format_twas_gsea_drugtargetor_results:
  wildcard_constraints:
    weight="(?!nondir_).+"
  input:
    "{outdir}/results/{gwas}/twas/drugtargetor/twas_gsea_drugtargetor_{weight}.done",
    rules.download_atc.output
  output:
    "{outdir}/results/{gwas}/twas/drugtargetor/twas_gsea_{weight}_res_atc_res.csv"
  benchmark:
    "{outdir}/benchmarks/format_twas_gsea_drugtargetor_results_{gwas}_{weight}.tsv"
  conda:
    "../envs/main.yaml"
  params:
    config_file=config['config_file']
  log:
    "{outdir}/logs/format_twas_gsea_drugtargetor_results-{gwas}-{weight}.log"
  shell:
    "Rscript --vanilla scripts/format_twas_gsea_drugtargetor_results.R \
    --twas {wildcards.gwas} \
    --panel {wildcards.weight} \
    --config_file {params.config_file} > {log} 2>&1"

rule format_twas_gsea_drugtargetor_results_all_panel:
    input:
      lambda w: expand("{outdir}/results/{gwas}/twas/drugtargetor/twas_gsea_{weight}_res_atc_res.csv", gwas=w.gwas, batch=range(1, 11), weight=weights_nosplice, outdir={outdir})
    output:
      touch("{outdir}/results/{gwas}/checks/format_twas_gsea_drugtargetor_results_all_panel.done")

#######
# Run TWAS-GSEA using DrugTargetor sets, non-directional (full .gmt, comparable to MAGMA)
#######

# Non-directional variant of run_twas_gsea_drug_targetor: uses the full .gmt
# (no signed drug effect) and disables --directional so the outcome is
# probit(1 - TWAS.P) instead of signed TWAS.Z. Reuses the same precomputed
# CorMat and TWAS results.
rule run_twas_gsea_drug_targetor_nondirectional:
  resources:
    mem_mb=50000,
    cpus=5
  input:
    rules.install_twas_gsea.output,
    "{outdir}/results/{gwas}/twas/{gwas}_twas_{weight}_GW_clean.txt.gz",
    rules.format_drug_targetor.output,
    f"{resdir}/data/predicted_expression/{{weight}}/Reference_Expression/{{weight}}.CorMat.RDS"
  output:
    touch("{outdir}/results/{gwas}/twas/drugtargetor/twas_gsea_drugtargetor_nondir_{weight}.done")
  benchmark:
    "{outdir}/benchmarks/run_twas_gsea_drug_targetor_nondirectional_{gwas}_{weight}.tsv"
  conda:
    "../envs/main.yaml"
  params:
    resdir=resdir
  log:
    "{outdir}/logs/run_twas_gsea_drug_targetor_nondirectional-{gwas}-{weight}.log"
  shell:
    "Rscript --vanilla {params.resdir}/software/TWAS-GSEA/TWAS-GSEA-fast.R \
      --twas_results {outdir}/results/{wildcards.gwas}/twas/{wildcards.gwas}_twas_{wildcards.weight}_GW_clean.txt.gz \
      --pos {params.resdir}/data/fusion_snp_weights/{wildcards.weight}/{wildcards.weight}.pos \
      --input_CorMat {params.resdir}/data/predicted_expression/{wildcards.weight}/Reference_Expression/{wildcards.weight}.CorMat.RDS \
      --gmt_file {params.resdir}/data/drug_targetor/wholedatabase_for_targetor_symbols.gmt \
      --n_cores 5 \
      --covar GeneLength,NSNP \
      --use_alt_id ID \
      --min_Ngenes 2 \
      --directional F \
      --output {outdir}/results/{wildcards.gwas}/twas/drugtargetor/twas_gsea_drugtargetor_nondir_{wildcards.weight} > {log} 2>&1"

# Format the non-directional output (parallel to format_twas_gsea_drugtargetor_results)
rule format_twas_gsea_drugtargetor_nondirectional_results:
  input:
    "{outdir}/results/{gwas}/twas/drugtargetor/twas_gsea_drugtargetor_nondir_{weight}.done",
    rules.download_atc.output
  output:
    "{outdir}/results/{gwas}/twas/drugtargetor/twas_gsea_nondir_{weight}_res_atc_res.csv"
  benchmark:
    "{outdir}/benchmarks/format_twas_gsea_drugtargetor_nondirectional_results_{gwas}_{weight}.tsv"
  conda:
    "../envs/main.yaml"
  params:
    config_file=config['config_file']
  log:
    "{outdir}/logs/format_twas_gsea_drugtargetor_nondirectional_results-{gwas}-{weight}.log"
  shell:
    "Rscript --vanilla scripts/format_twas_gsea_drugtargetor_results.R \
    --twas {wildcards.gwas} \
    --panel {wildcards.weight} \
    --mode nondirectional \
    --config_file {params.config_file} > {log} 2>&1"

rule format_twas_gsea_drugtargetor_nondirectional_results_all_panel:
    input:
      lambda w: expand("{outdir}/results/{gwas}/twas/drugtargetor/twas_gsea_nondir_{weight}_res_atc_res.csv", gwas=w.gwas, weight=weights_nosplice, outdir={outdir})
    output:
      touch("{outdir}/results/{gwas}/checks/format_twas_gsea_drugtargetor_nondirectional_results_all_panel.done")

#######
# Run TWAS-GSEA against reprocessed CMAP level5 drug signatures (directional)
#######

# Same shape as run_twas_gsea_drug_targetor, but using a user-supplied prop file
# generated from CMAP level5 (see DrugRepurposing/Code/Published/make_twas_gsea_prop.R).
# The prop file is keyed on HGNC gene symbol to match the ID column produced by
# scripts/combine_twas.R.
rule run_twas_gsea_cmap:
  wildcard_constraints:
    weight="(?!nondir_).+"
  resources:
    mem_mb=100000,
    cpus=8
  input:
    rules.install_twas_gsea.output,
    "{outdir}/results/{gwas}/twas/{gwas}_twas_{weight}_GW_clean.txt.gz",
    f"{resdir}/data/predicted_expression/{{weight}}/Reference_Expression/{{weight}}.CorMat.RDS",
    config["cmap_level5_prop"] if config["twas_gsea_cmap"] == "T" else []
  output:
    touch("{outdir}/results/{gwas}/twas/cmap/twas_gsea_cmap_{weight}.done")
  benchmark:
    "{outdir}/benchmarks/run_twas_gsea_cmap_{gwas}_{weight}.tsv"
  conda:
    "../envs/main.yaml"
  params:
    resdir=resdir,
    prop=config["cmap_level5_prop"]
  log:
    "{outdir}/logs/run_twas_gsea_cmap-{gwas}-{weight}.log"
  shell:
    "Rscript --vanilla {params.resdir}/software/TWAS-GSEA/TWAS-GSEA-fast.R \
      --twas_results {outdir}/results/{wildcards.gwas}/twas/{wildcards.gwas}_twas_{wildcards.weight}_GW_clean.txt.gz \
      --pos {params.resdir}/data/fusion_snp_weights/{wildcards.weight}/{wildcards.weight}.pos \
      --input_CorMat {params.resdir}/data/predicted_expression/{wildcards.weight}/Reference_Expression/{wildcards.weight}.CorMat.RDS \
      --prop_file {params.prop} \
      --n_cores 8 \
      --covar GeneLength,NSNP \
      --use_alt_id ID \
      --min_Ngenes 2 \
      --directional T \
      --output {outdir}/results/{wildcards.gwas}/twas/cmap/twas_gsea_cmap_{wildcards.weight} > {log} 2>&1"

rule run_twas_gsea_cmap_all_panel:
    input:
      lambda w: expand("{outdir}/results/{gwas}/twas/cmap/twas_gsea_cmap_{weight}.done", gwas=w.gwas, weight=weights_nosplice, outdir={outdir})
    output:
      touch("{outdir}/results/{gwas}/checks/run_twas_gsea_cmap_all_panel.done")

# Format the per-(gwas, weight) CMAP TWAS-GSEA output: parse signature names,
# attach MOA from compoundinfo_beta.txt, and write per-signature + per-MOA CSVs.
rule format_twas_gsea_cmap_results:
  input:
    "{outdir}/results/{gwas}/twas/cmap/twas_gsea_cmap_{weight}.done"
  output:
    "{outdir}/results/{gwas}/twas/cmap/twas_gsea_cmap_{weight}_drug_res.csv",
    "{outdir}/results/{gwas}/twas/cmap/twas_gsea_cmap_{weight}_moa_res.csv"
  benchmark:
    "{outdir}/benchmarks/format_twas_gsea_cmap_results_{gwas}_{weight}.tsv"
  conda:
    "../envs/main.yaml"
  params:
    config_file=config['config_file']
  log:
    "{outdir}/logs/format_twas_gsea_cmap_results-{gwas}-{weight}.log"
  shell:
    "Rscript --vanilla scripts/format_twas_gsea_cmap_results.R \
    --twas {wildcards.gwas} \
    --panel {wildcards.weight} \
    --config_file {params.config_file} > {log} 2>&1"

rule format_twas_gsea_cmap_results_all_panel:
    input:
      lambda w: expand("{outdir}/results/{gwas}/twas/cmap/twas_gsea_cmap_{weight}_drug_res.csv", gwas=w.gwas, weight=weights_nosplice, outdir={outdir})
    output:
      touch("{outdir}/results/{gwas}/checks/format_twas_gsea_cmap_results_all_panel.done")
