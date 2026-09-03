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

# Path to bivariate LDSC secondary list (NA when not configured).
# Read only as a path for input-staleness tracking; per-pair file existence is
# verified at runtime so a missing secondary fails in isolation (NA row) rather
# than at parse time.
gencor_gwas_list_path = config.get("gencor_gwas_list", "NA")

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
    f"{outdir}/results/{{gwas}}/gwas_sumstat/{{gwas}}.cleaned.gz",
    f"{outdir}/results/{{gwas}}/gwas_sumstat/{{gwas}}.cleaned.munged.sumstats.gz"
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
    resdir=resdir,
    # When the run is restricted to a single chromosome (e.g. chromosomes: [22]
    # for test mode, or a chr-restricted input like the chr22 demo), pass it to
    # sumstat_cleaner's --test so ref_harmonise only harmonises that chromosome.
    # Without it, ref_harmonise defaults to chr 1:22 and merge()s an empty target
    # chunk for the absent chromosomes -> "'by' must specify a uniquely valid
    # column". Omitted for full genome-wide runs (unchanged behaviour).
    test_arg=("--test " + str(chromosomes[0])) if len(chromosomes) == 1 else ""
  log:
    f"{outdir}/logs/sumstat_prep_i-{{gwas}}.log"
  shell:
    """
    (sumstat_cleaner_script=$(Rscript --vanilla -e 'cat(system.file("scripts", "sumstat_cleaner.R", package = "GenoUtils"))')
    Rscript --vanilla $sumstat_cleaner_script \
      --sumstats {params.path} \
      --n {params.n} \
      --ref_chr {params.resdir}/data/1kg/1KG.Phase3.MAF_001.chr \
      --population {params.population} \
      --munged T \
      {params.test_arg} \
      --output {outdir}/results/{wildcards.gwas}/gwas_sumstat/{wildcards.gwas}.cleaned) > {log} 2>&1
    """

rule sumstat_prep:
  input: expand(f"{outdir}/results/{{gwas}}/gwas_sumstat/{{gwas}}.cleaned.gz", gwas=gwas_list_df_eur['name'])

###
# Calculate median effective sample size
###

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
    "Rscript --vanilla {workflow.basedir}/scripts/median_n.R --pipeline_dir {workflow.basedir} --munged {input} --out {output} > {log} 2>&1"

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
# Bivariate LDSC (genetic correlation)
###

rule ldsc_gencor:
  input:
    "{outdir}/results/{gwas}/gwas_sumstat/{gwas}.cleaned.munged.sumstats.gz",
    f"{resdir}/software/ldsc/",
    f"{resdir}/data/ldsc/eur_w_ld_chr/10.l2.ldscore.gz",
    f"{resdir}/data/ldsc/w_hm3.snplist",
    gencor_gwas_list_path
  output:
    touch("{outdir}/results/{gwas}/gencor/{gwas}_gencor_pairs.done")
  benchmark:
    "{outdir}/benchmarks/ldsc_gencor_{gwas}.tsv"
  conda:
    "../envs/ldsc.yaml"
  params:
    resdir=resdir,
    gencor_list=gencor_gwas_list_path
  log:
    "{outdir}/logs/ldsc_gencor-{gwas}.log"
  shell:
    """
    (mkdir -p {outdir}/results/{wildcards.gwas}/gencor/

    # Pre-check: every secondary file in the path column must exist. If any are
    # missing we abort the rule rather than write NA rows, so typos in the
    # gencor_gwas_list are surfaced immediately. Missing-secondary is treated
    # as a configuration error; other per-pair runtime failures (e.g. LDSC
    # crashes on malformed sumstats) still record NA via the || fallback below.
    missing=$(tail -n +2 {params.gencor_list} | awk '{{print $1, $2}}' | while read name path; do
      [ -f "${{path}}" ] || echo "  ${{name}} -> ${{path}}"
    done)
    if [ -n "$missing" ]; then
      echo "ERROR: secondary GWAS files listed in {params.gencor_list} are missing:"
      echo "$missing"
      exit 1
    fi

    tail -n +2 {params.gencor_list} | awk '{{print $1, $2}}' | while read name path; do
      out_prefix={outdir}/results/{wildcards.gwas}/gencor/{wildcards.gwas}__${{name}}
      python2.7 {params.resdir}/software/ldsc/ldsc.py \
        --rg {outdir}/results/{wildcards.gwas}/gwas_sumstat/{wildcards.gwas}.cleaned.munged.sumstats.gz,${{path}} \
        --ref-ld-chr {params.resdir}/data/ldsc/eur_w_ld_chr/ \
        --w-ld-chr {params.resdir}/data/ldsc/eur_w_ld_chr/ \
        --out ${{out_prefix}} \
        || echo "GENCOR_PAIR_FAILED: ${{name}}" >> ${{out_prefix}}.log
    done) > {log} 2>&1
    """

rule process_ldsc_gencor:
  input:
    "{outdir}/results/{gwas}/gencor/{gwas}_gencor_pairs.done"
  output:
    "{outdir}/results/{gwas}/gencor/{gwas}_gencor_res.csv"
  benchmark:
    "{outdir}/benchmarks/process_ldsc_gencor_{gwas}.tsv"
  conda:
    "../envs/main.yaml"
  params:
    config_file=config['config_file']
  log:
    "{outdir}/logs/process_ldsc_gencor-{gwas}.log"
  shell:
    "Rscript --vanilla {workflow.basedir}/scripts/process_ldsc_gencor.R --pipeline_dir {workflow.basedir} \
      --gwas {wildcards.gwas} \
      --config_file {params.config_file} > {log} 2>&1"

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
    "Rscript --vanilla {workflow.basedir}/scripts/process_clump.R --pipeline_dir {workflow.basedir} \
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
    """
    # GCTA-COJO aborts a chromosome with "Error: too many SNPs..." when the number of
    # independent genome-wide-significant signals exceeds the LD reference sample size
    # (1KG-EUR, ~503). Tolerate ONLY that error so the pipeline can continue with the
    # chromosomes that succeeded (recording a per-chr status); any other GCTA failure
    # still exits non-zero and fails the rule as before.
    mkdir -p {outdir}/results/{wildcards.gwas}/cojo
    sf={outdir}/results/{wildcards.gwas}/cojo/{wildcards.gwas}_chr{wildcards.chr}.cojo.status
    if {params.resdir}/software/gcta/gcta_v1.94.0Beta_linux_kernel_3_x86_64/gcta_v1.94.0Beta_linux_kernel_3_x86_64_static \
        --bfile {params.resdir}/data/1kg/1KG.Phase3.{params.population}.MAF_001.chr{wildcards.chr} \
        --chr {wildcards.chr} \
        --maf 0.01 \
        --cojo-file {outdir}/results/{wildcards.gwas}/gwas_sumstat/{wildcards.gwas}.cleaned.cojo \
        --cojo-slct \
        --cojo-p 5e-8 \
        --out {outdir}/results/{wildcards.gwas}/cojo/{wildcards.gwas}_chr{wildcards.chr} > {log} 2>&1; then
      echo ok > "$sf"
    elif grep -q 'too many SNPs' {log}; then
      echo reference_too_small > "$sf"
      echo 'NOTE: COJO skipped for chr{wildcards.chr} - independent signals exceed LD reference sample size.' >> {log}
    else
      exit 1
    fi
    """

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
    "Rscript --vanilla {workflow.basedir}/scripts/process_cojo.R --pipeline_dir {workflow.basedir} \
      --gwas {wildcards.gwas} \
      --config_file {params.config_file} > {log} 2>&1"

###
# QQ plot
###

rule qq_plot:
  resources:
    mem_mb=lambda wildcards, input: max(
      4000,
      int(2000 + (8 * os.path.getsize(input[0]) / 1024**2))
    )
  input:
    f"{outdir}/results/{{gwas}}/gwas_sumstat/{{gwas}}.cleaned.gz"
  output:
    f"{outdir}/results/{{gwas}}/gwas_sumstat/{{gwas}}.qq_plot.png"
  benchmark:
    f"{outdir}/benchmarks/qq_plot_{{gwas}}.tsv"
  conda:
    "../envs/main.yaml"
  params:
    config_file=config['config_file']
  log:
    f"{outdir}/logs/qq_plot-{{gwas}}.log"
  shell:
    "Rscript --vanilla {workflow.basedir}/scripts/plot_qq.R --pipeline_dir {workflow.basedir} \
      --gwas {wildcards.gwas} \
      --config_file {params.config_file} > {log} 2>&1"

###
# Manhattan plot
###

def _manhattan_inputs(wildcards):
    inputs = {"sumstats": f"{outdir}/results/{wildcards.gwas}/gwas_sumstat/{wildcards.gwas}.cleaned.gz"}
    if config.get("clump", "F") == "T":
        inputs["clump_csv"] = f"{outdir}/results/{wildcards.gwas}/clump/{wildcards.gwas}.GW.clump.clean.csv"
    return inputs

rule manhattan_plot:
  resources:
    mem_mb=lambda wildcards, input: max(
      6000,
      int(3000 + (12 * os.path.getsize(input.sumstats) / 1024**2))
    )
  input:
    unpack(_manhattan_inputs)
  output:
    unlabelled=f"{outdir}/results/{{gwas}}/gwas_sumstat/{{gwas}}.manhattan_plot.unlabelled.png",
    labelled=f"{outdir}/results/{{gwas}}/gwas_sumstat/{{gwas}}.manhattan_plot.labelled.png",
    data=f"{outdir}/results/{{gwas}}/gwas_sumstat/{{gwas}}.manhattan_data.rds"
  benchmark:
    f"{outdir}/benchmarks/manhattan_plot_{{gwas}}.tsv"
  conda:
    "../envs/main.yaml"
  params:
    config_file=config['config_file']
  log:
    f"{outdir}/logs/manhattan_plot-{{gwas}}.log"
  shell:
    "Rscript --vanilla {workflow.basedir}/scripts/plot_manhattan.R --pipeline_dir {workflow.basedir} \
      --gwas {wildcards.gwas} \
      --config_file {params.config_file} > {log} 2>&1"

###
# Per-locus zoom plots
###

rule locus_plots:
  resources:
    mem_mb=lambda wildcards, input: max(
      6000,
      int(3000 + (12 * os.path.getsize(input.sumstats) / 1024**2))
    )
  input:
    sumstats=f"{outdir}/results/{{gwas}}/gwas_sumstat/{{gwas}}.cleaned.gz",
    clump_csv=f"{outdir}/results/{{gwas}}/clump/{{gwas}}.GW.clump.clean.csv",
    recomb_maps=expand(f"{resdir}/data/recomb_maps/chr{{chr}}.interpolated_genetic_map.gz",
                        chr=chromosomes)
  output:
    touch(f"{outdir}/results/{{gwas}}/locus_plots/{{gwas}}.locus_plots.done")
  benchmark:
    f"{outdir}/benchmarks/locus_plots_{{gwas}}.tsv"
  conda:
    "../envs/main.yaml"
  params:
    config_file=config['config_file']
  log:
    f"{outdir}/logs/locus_plots-{{gwas}}.log"
  shell:
    "Rscript --vanilla {workflow.basedir}/scripts/plot_locus_zoom.R --pipeline_dir {workflow.basedir} \
      --gwas {wildcards.gwas} \
      --config_file {params.config_file} > {log} 2>&1"
