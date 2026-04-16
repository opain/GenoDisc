##########
# Analyse GWAS summary statistics
##########

##
# QC and format GWAS summary statistics
##

# Read in GWAS list
gwas_list_df = pd.read_table(config["gwas_list"], sep=' ')

# Subset to EUR-based GWAS
gwas_list_df_eur = gwas_list_df.loc[gwas_list_df['population'] == 'EUR']

rule sumstat_prep_i:
  resources:
    mem_mb=lambda wildcards, input: max(
      8000,
      int(3500 + (16 * os.path.getsize(input[2]) / 1024**2))
    )
  input:
    rules.prep_1kg.output,
    rules.install_genoutils.output,
    lambda w: gwas_list_df.loc[gwas_list_df['name'] == "{}".format(w.gwas), 'path'].iloc[0]
  output:
    f"{outdir}/results/{{gwas}}/gwas_sumstat/{{gwas}}.cleaned.gz"
  benchmark:
    f"{outdir}/benchmarks/sumstat_prep_i_{{gwas}}.tsv"
  conda:
    "../envs/main.yaml"
  params:
    outdir=config["outdir"],
    config_file = config["config_file"],
    population= lambda w: gwas_list_df.loc[gwas_list_df['name'] == "{}".format(w.gwas), 'population'].iloc[0],
    n= lambda w: gwas_list_df.loc[gwas_list_df['name'] == "{}".format(w.gwas), 'n'].iloc[0],
    path= lambda w: gwas_list_df.loc[gwas_list_df['name'] == "{}".format(w.gwas), 'path'].iloc[0],
    resdir=resdir
  log:
    f"{outdir}/logs/sumstat_prep_i-{{gwas}}.log"
  shell:
    """
    (sumstat_cleaner_script=$(Rscript -e 'cat(system.file("scripts", "sumstat_cleaner.R", package = "GenoUtils"))')
    Rscript $sumstat_cleaner_script \
      --sumstats {params.path} \
      --n {params.n} \
      --ref_chr {params.resdir}/data/1kg/1KG.Phase3.MAF_001.chr \
      --population {params.population} \
      --output {outdir}/results/{wildcards.gwas}/gwas_sumstat/{wildcards.gwas}.cleaned) > {log} 2>&1
    """

rule sumstat_prep:
  input: expand(f"{outdir}/results/{{gwas}}/gwas_sumstat/{{gwas}}.cleaned.gz", gwas=gwas_list_df_eur['name'])

###
# Munge sumstats
###

# munge sumstats using FOCUS munge function
rule focus_munge:
  input:
    premunged="{outdir}/results/{gwas}/gwas_sumstat/{gwas}.cleaned.gz"
  output:
    "{outdir}/results/{gwas}/gwas_sumstat/{gwas}.cleaned.munged.sumstats.gz"
  benchmark:
    "{outdir}/benchmarks/focus_munge_{gwas}.tsv"
  conda:
    "../envs/focus.yaml"
  log:
    "{outdir}/logs/focus_munge-{gwas}.log"
  shell:
    "focus munge {input.premunged} --output {outdir}/results/{wildcards.gwas}/gwas_sumstat/{wildcards.gwas}.cleaned.munged > {log} 2>&1"

# Calculate median effective sample size
# FUSION requires this parameter to be specified despite having the N column in the sumstats
rule retrieve_N:
  input:
    "{outdir}/results/{gwas}/gwas_sumstat/{gwas}.cleaned.munged.sumstats.gz"
  output:
    "{outdir}/results/{gwas}/gwas_sumstat/{gwas}.cleaned.munged.median_N.txt"
  benchmark:
    "{outdir}/benchmarks/retrieve_N_{gwas}.tsv"
  conda:
    "../envs/main.yaml"
  log:
    "{outdir}/logs/retrieve_N-{gwas}.log"
  shell:
    "Rscript scripts/median_n.R --munged {input} --out {output} > {log} 2>&1"

###
# Run LDSC
###

rule ldsc:
  input:
    "{outdir}/results/{gwas}/gwas_sumstat/{gwas}.cleaned.munged.sumstats.gz",
    f"{resdir}/software/ldsc/",
    f"{resdir}/data/ldsc/eur_w_ld_chr/10.l2.ldscore.gz",
    f"{resdir}/data/ldsc/w_hm3.snplist"
  output:
    "{outdir}/results/{gwas}/ldsc/{gwas}_ldsc_res.log"
  benchmark:
    "{outdir}/benchmarks/ldsc_{gwas}.tsv"
  conda:
    "../envs/ldsc.yaml"
  params:
    resdir=resdir
  log:
    "{outdir}/logs/ldsc-{gwas}.log"
  shell:
    "(mkdir -p {outdir}/results/{wildcards.gwas}/ldsc/; python2.7 {params.resdir}/software/ldsc/ldsc.py \
      --h2 {outdir}/results/{wildcards.gwas}/gwas_sumstat/{wildcards.gwas}.cleaned.munged.sumstats.gz \
      --ref-ld-chr {params.resdir}/data/ldsc/eur_w_ld_chr/ \
      --w-ld-chr {params.resdir}/data/ldsc/eur_w_ld_chr/ \
      --out {outdir}/results/{wildcards.gwas}/ldsc/{wildcards.gwas}_ldsc_res) > {log} 2>&1"

###
# Run LD clumping
###

rule clump:
  input:
    "{outdir}/results/{gwas}/gwas_sumstat/{gwas}.cleaned.gz"
  output:
    touch("{outdir}/results/{gwas}/checks/{gwas}_chr{chr}.clumped.done")
  benchmark:
    "{outdir}/benchmarks/clump_{gwas}_chr{chr}.tsv"
  conda:
    "../envs/main.yaml"
  params:
    population= lambda w: gwas_list_df_eur.loc[gwas_list_df_eur['name'] == "{}".format(w.gwas), 'population'].iloc[0],
    resdir=resdir
  log:
    "{outdir}/logs/clump-{gwas}-chr{chr}.log"
  shell:
    "(mkdir -p {outdir}/results/{wildcards.gwas}/clump; plink \
      --bfile {params.resdir}/data/1kg/1KG.Phase3.{params.population}.MAF_001.chr{wildcards.chr} \
      --chr {wildcards.chr} \
      --maf 0.01 \
      --clump {outdir}/results/{wildcards.gwas}/gwas_sumstat/{wildcards.gwas}.cleaned.gz \
      --clump-p1 1e-5 \
      --clump-r2 0.1 \
      --clump-kb 500 \
      --out {outdir}/results/{wildcards.gwas}/clump/{wildcards.gwas}_chr{wildcards.chr}) > {log} 2>&1"

rule clump_all_chr:
    input:
      lambda w: expand("{outdir}/results/{gwas}/checks/{gwas}_chr{chr}.clumped.done", gwas=w.gwas, chr=chromosomes, outdir={outdir})
    output:
      touch("{outdir}/results/{gwas}/checks/clump_all_chr.done")

###
# Process clumping results
###

rule process_clump:
  input:
    "{outdir}/results/{gwas}/checks/clump_all_chr.done",
    rules.download_biomart.output
  output:
    "{outdir}/results/{gwas}/clump/{gwas}.GW.clump.clean.csv"
  benchmark:
    "{outdir}/benchmarks/process_clump_{gwas}.tsv"
  conda:
    "../envs/main.yaml"
  params:
    config_file=config['config_file']
  log:
    "{outdir}/logs/process_clump-{gwas}.log"
  shell:
    "Rscript scripts/process_clump.R \
      --gwas {wildcards.gwas} \
      --config_file {params.config_file} > {log} 2>&1"

###
# Run COJO
###
# Note this is not ideal as we are using the EUR subset of 1KG
# The developers suggest using a reference that matches a large sample in the GWAS
# Or at at least a reference that is >4000 individuals
# Should allow the user to specify a reference of their own

rule cojo:
  input:
    rules.download_gcta.output,
    "{outdir}/results/{gwas}/gwas_sumstat/{gwas}.cleaned.cojo"
  output:
    touch("{outdir}/results/{gwas}/checks/{gwas}_cojo_chr{chr}.done")
  benchmark:
    "{outdir}/benchmarks/cojo_{gwas}_chr{chr}.tsv"
  conda:
    "../envs/ldsc.yaml"
  params:
    population= lambda w: gwas_list_df_eur.loc[gwas_list_df_eur['name'] == "{}".format(w.gwas), 'population'].iloc[0],
    resdir=resdir
  log:
    "{outdir}/logs/cojo-{gwas}-chr{chr}.log"
  shell:
    "(mkdir -p {outdir}/results/{wildcards.gwas}/cojo; {params.resdir}/software/gcta/gcta_v1.94.0Beta_linux_kernel_3_x86_64/gcta_v1.94.0Beta_linux_kernel_3_x86_64_static \
      --bfile {params.resdir}/data/1kg/1KG.Phase3.{params.population}.MAF_001.chr{wildcards.chr} \
      --chr {wildcards.chr} \
      --maf 0.01 \
      --cojo-file {outdir}/results/{wildcards.gwas}/gwas_sumstat/{wildcards.gwas}.cleaned.cojo \
      --cojo-slct \
      --cojo-p 1e-5 \
      --out {outdir}/results/{wildcards.gwas}/cojo/{wildcards.gwas}_chr{wildcards.chr}) > {log} 2>&1"

rule cojo_all_chr:
    input:
      lambda w: expand("{outdir}/results/{gwas}/checks/{gwas}_cojo_chr{chr}.done", gwas=w.gwas, chr=chromosomes, outdir={outdir})
    output:
      touch("{outdir}/results/{gwas}/checks/cojo_all_chr.done")

###
# Process COJO results
###

rule process_cojo:
  input:
    "{outdir}/results/{gwas}/checks/cojo_all_chr.done",
    rules.download_biomart.output
  output:
    "{outdir}/results/{gwas}/cojo/{gwas}.GW.cojo.clean.csv"
  benchmark:
    "{outdir}/benchmarks/process_cojo_{gwas}.tsv"
  conda:
    "../envs/main.yaml"
  params:
    config_file=config['config_file']
  log:
    "{outdir}/logs/process_cojo-{gwas}.log"
  shell:
    "Rscript scripts/process_cojo.R \
      --gwas {wildcards.gwas} \
      --config_file {params.config_file} > {log} 2>&1"
