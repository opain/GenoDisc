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
    "../envs/main.yaml"
  shell:
    "mkdir -p resources/software/liftover/; wget --no-check-certificate -O resources/software/liftover/liftover https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver"

# Download liftover track
rule download_liftover_track:
  output:
    touch("resources/software/download_liftover_track.done")
  conda:
    "../envs/main.yaml"
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
    "../envs/main.yaml"
  shell:
    "Rscript scripts/prep_1kg.R"

####
# Install FOCUS
####

rule install_focus:
  output:
    directory("resources/software/ma-focus/")
  conda:
    "../envs/focus.yaml"
  shell:
    "git clone https://github.com/mancusolab/ma-focus.git {output}; \
    cd {output}; \
    git reset --hard 8af424a2d38222f76bf7a0422cce8acf274dc610; \
    python3 -m pip install ."

####
# Download LDSC
####

# Install LDSC
rule install_ldsc:
  output:
    directory("resources/software/ldsc/")
  conda:
    "../envs/main.yaml"
  shell:
    "git clone https://github.com/bulik/ldsc.git {output}; \
    cd {output}; \
    git reset --hard aa33296abac9569a6422ee6ba7eb4b902422cc74"

# Download LDSC reference data
rule download_ldsc_scores:
  output:
    "resources/data/ldsc/eur_w_ld_chr/10.l2.ldscore.gz"
  conda:
    "../envs/main.yaml"
  shell:
    "mkdir -p resources/data/ldsc; wget --no-check-certificate -O resources/data/ldsc/eur_w_ld_chr.tar.gz https://zenodo.org/record/8182036/files/eur_w_ld_chr.tar.gz?download=1; tar -xf resources/data/ldsc/eur_w_ld_chr.tar.gz -C resources/data/ldsc; rm resources/data/ldsc/eur_w_ld_chr.tar.gz"

rule download_ldsc_hm3:
  output:
    "resources/data/ldsc/w_hm3.snplist"
  conda:
    "../envs/main.yaml"
  shell:
    "mkdir -p resources/data/ldsc; wget --no-check-certificate -O resources/data/ldsc/w_hm3.snplist.gz https://zenodo.org/record/7773502/files/w_hm3.snplist.gz?download=1; gzip -d resources/data/ldsc/w_hm3.snplist.gz"

# Download GCTA
rule download_gcta:
  output:
    directory('resources/software/gcta/gcta_v1.94.0Beta_linux_kernel_3_x86_64/')
  conda:
    "../envs/main.yaml"
  shell:
    "mkdir -p resources/software/gcta; wget -O resources/software/gcta/gcta_v1.94.0Beta_linux_kernel_3_x86_64.zip https://yanglab.westlake.edu.cn/software/gcta/bin/gcta_v1.94.0Beta_linux_kernel_3_x86_64.zip; unzip resources/software/gcta/gcta_v1.94.0Beta_linux_kernel_3_x86_64.zip -d resources/software/gcta; rm resources/software/gcta/gcta_v1.94.0Beta_linux_kernel_3_x86_64.zip"

# Install GenoUtils
rule install_genoutils:
  input:
    "envs/main.yaml"
  output:
    touch("resources/software/install_genoutils.done")
  conda:
    "../envs/main.yaml"
  shell:
    "Rscript -e 'devtools::install_github(\"opain/GenoUtils@4beb75620f3291b633598acd06febb22298418c8\")'"

##########
# Analyse GWAS summary statistics
##########

# For the time being, assume the GWAS sumstats are in Rosalind format
##
# QC and format GWAS summary statistics
##

# Read in GWAS list
gwas_list_df = pd.read_table(config["gwas_list"], sep=' ')

# Subset to EUR-based GWAS
gwas_list_df_eur = gwas_list_df.loc[gwas_list_df['population'] == 'EUR']

rule sumstat_prep_i:
  input:
    rules.prep_1kg.output,
    rules.install_genoutils.output,
    lambda w: gwas_list_df.loc[gwas_list_df['name'] == "{}".format(w.gwas), 'path'].iloc[0]
  output:
    f"{outdir}/data/gwas_sumstat/{{gwas}}/{{gwas}}.cleaned.gz"
  conda:
    "../envs/main.yaml"
  params:
    outdir=config["outdir"],
    config_file = config["config_file"],
    population= lambda w: gwas_list_df.loc[gwas_list_df['name'] == "{}".format(w.gwas), 'population'].iloc[0],
    n= lambda w: gwas_list_df.loc[gwas_list_df['name'] == "{}".format(w.gwas), 'n'].iloc[0],
    path= lambda w: gwas_list_df.loc[gwas_list_df['name'] == "{}".format(w.gwas), 'path'].iloc[0]
  shell:
    """
    sumstat_cleaner_script=$(Rscript -e 'cat(system.file("scripts", "sumstat_cleaner.R", package = "GenoUtils"))')
    Rscript $sumstat_cleaner_script \
      --sumstats {params.path} \
      --n {params.n} \
      --ref_chr resources/data/1kg/1KG.Phase3.MAF_001.chr \
      --population {params.population} \
      --output {outdir}/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}.cleaned
    """

rule sumstat_prep:
  input: expand(f"{outdir}/data/gwas_sumstat/{{gwas}}/{{gwas}}.cleaned.gz", gwas=gwas_list_df_eur['name'])

###
# Munge sumstats
###

# munge sumstats using FOCUS munge function
rule focus_munge:
  input:
    premunged="{outdir}/data/gwas_sumstat/{gwas}/{gwas}.cleaned.gz",
    focus=rules.install_focus.output
  output:
    "{outdir}/data/gwas_sumstat/{gwas}/{gwas}.cleaned.munged.sumstats.gz"
  conda:
    "../envs/focus.yaml"
  shell:
    "focus munge {input.premunged} --output {outdir}/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}.cleaned.munged"

# Calculate median effective sample size
# FUSION requires this parameter to be specified despite having the N column in the sumstats
rule retrieve_N:
  input:
    "{outdir}/data/gwas_sumstat/{gwas}/{gwas}.cleaned.munged.sumstats.gz"
  output:
    "{outdir}/data/gwas_sumstat/{gwas}/{gwas}.cleaned.munged.median_N.txt"
  conda:
    "../envs/main.yaml"
  shell:
    "Rscript scripts/median_n.R --munged {input} --out {output}"

###
# Run LDSC
###

rule ldsc:
  input:
    "{outdir}/data/gwas_sumstat/{gwas}/{gwas}.cleaned.munged.sumstats.gz",
    "resources/software/ldsc/",
    "resources/data/ldsc/eur_w_ld_chr/10.l2.ldscore.gz",
    "resources/data/ldsc/w_hm3.snplist"
  output:
    "{outdir}/results/{gwas}/ldsc/{gwas}_ldsc_res.log"
  conda:
    "../envs/ldsc.yaml"
  shell:
    "mkdir -p {outdir}/results/{wildcards.gwas}/ldsc/; python2.7 resources/software/ldsc/ldsc.py \
      --h2 {outdir}/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}.cleaned.munged.sumstats.gz \
      --ref-ld-chr resources/data/ldsc/eur_w_ld_chr/ \
      --w-ld-chr resources/data/ldsc/eur_w_ld_chr/ \
      --out {outdir}/results/{wildcards.gwas}/ldsc/{wildcards.gwas}_ldsc_res"

###
# Run LD clumping
###

rule clump:
  input:
    "{outdir}/data/gwas_sumstat/{gwas}/{gwas}.cleaned.gz"
  output:
    touch("{outdir}/results/{gwas}/checks/{gwas}_chr{chr}.clumped.done")
  conda:
    "../envs/main.yaml"
  params:
    population= lambda w: gwas_list_df_eur.loc[gwas_list_df_eur['name'] == "{}".format(w.gwas), 'population'].iloc[0]
  shell:
    "mkdir -p {outdir}/results/{wildcards.gwas}/clump; plink \
      --bfile resources/data/1kg/1KG.Phase3.{params.population}.MAF_001.chr{wildcards.chr} \
      --chr {wildcards.chr} \
      --maf 0.01 \
      --clump {outdir}/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}.cleaned.gz \
      --clump-p1 1e-5 \
      --clump-r2 0.1 \
      --clump-kb 500 \
      --out {outdir}/results/{wildcards.gwas}/clump/{wildcards.gwas}_chr{wildcards.chr}"

rule clump_all_chr:
    input:
      lambda w: expand("{outdir}/results/{gwas}/checks/{gwas}_chr{chr}.clumped.done", gwas=w.gwas, chr=range(1, 23), outdir={outdir})
    output:
      touch("{outdir}/results/{gwas}/checks/clump_all_chr.done")

###
# Process clumping results
###

rule process_clump:
  input:
    "{outdir}/results/{gwas}/checks/clump_all_chr.done"
  output:
    "{outdir}/results/{gwas}/clump/{gwas}.GW.clump.clean.csv"
  conda:
    "../envs/main.yaml"
  params:
    config_file=config['config_file']
  shell:
    "Rscript scripts/process_clump.R \
      --gwas {wildcards.gwas} \
      --config_file {params.config_file}"

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
    "{outdir}/data/gwas_sumstat/{gwas}/{gwas}.cleaned.cojo"
  output:
    touch("{outdir}/results/{gwas}/checks/{gwas}_cojo_chr{chr}.done")
  conda:
    "../envs/ldsc.yaml"
  params:
    population= lambda w: gwas_list_df_eur.loc[gwas_list_df_eur['name'] == "{}".format(w.gwas), 'population'].iloc[0]
  shell:
    "mkdir -p {outdir}/results/{wildcards.gwas}/cojo; resources/software/gcta/gcta_v1.94.0Beta_linux_kernel_3_x86_64/gcta_v1.94.0Beta_linux_kernel_3_x86_64_static \
      --bfile resources/data/1kg/1KG.Phase3.{params.population}.MAF_001.chr{wildcards.chr} \
      --chr {wildcards.chr} \
      --maf 0.01 \
      --cojo-file {outdir}/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}.cleaned.cojo \
      --cojo-slct \
      --cojo-p 1e-5 \
      --out {outdir}/results/{wildcards.gwas}/cojo/{wildcards.gwas}_chr{wildcards.chr}"

rule cojo_all_chr:
    input:
      lambda w: expand("{outdir}/results/{gwas}/checks/{gwas}_cojo_chr{chr}.done", gwas=w.gwas, chr=range(1, 23), outdir={outdir})
    output:
      touch("{outdir}/results/{gwas}/checks/cojo_all_chr.done")

###
# Process COJO results
###

rule process_cojo:
  input:
    "{outdir}/results/{gwas}/checks/cojo_all_chr.done"
  output:
    "{outdir}/results/{gwas}/cojo/{gwas}.GW.cojo.clean.csv"
  conda:
    "../envs/main.yaml"
  params:
    config_file=config['config_file']
  shell:
    "Rscript scripts/process_cojo.R \
      --gwas {wildcards.gwas} \
      --config_file {params.config_file}"



