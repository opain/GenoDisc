####
# Perform TWAS using GTEx v7 weights
####

# run twas
rule run_twas_gcsc:
  resources:
    mem_mb=20000
  input:
    "{outdir}/results/{gwas}/gwas_sumstat/{gwas}.cleaned.munged.sumstats.gz",
    "{outdir}/results/{gwas}/gwas_sumstat/{gwas}.cleaned.munged.median_N.txt",
    rules.install_fusion.output,
    rules.install_plink2R.output,
    rules.prep_1kg.output,
    f"{resdir}/data/GCSC/coreg/{{gcsc_tissue}}_geneNames.txt",
    f"{resdir}/data/GCSC/twas_weights/GTEx.{{gcsc_tissue}}.P01"
  output:
    "{outdir}/results/{gwas}/gcsc/twas/{gcsc_tissue}/{gwas}_twas_{gcsc_tissue}_chr{chr}.dat"
  benchmark:
    "{outdir}/benchmarks/run_twas_gcsc_{gwas}_{gcsc_tissue}_chr{chr}.tsv"
  conda:
    "../envs/main.yaml"
  params:
    resdir=resdir
  log:
    "{outdir}/logs/run_twas_gcsc-{gwas}-{gcsc_tissue}-chr{chr}.log"
  shell:
    "(mkdir -p {outdir}/results/{wildcards.gwas}/gcsc/twas/{wildcards.gcsc_tissue}; N=$(cat {outdir}/results/{wildcards.gwas}/gwas_sumstat/{wildcards.gwas}.cleaned.munged.median_N.txt); Rscript {params.resdir}/software/fusion/FUSION.assoc_test.R \
    --sumstats {outdir}/results/{wildcards.gwas}/gwas_sumstat/{wildcards.gwas}.cleaned.munged.sumstats.gz \
    --weights {params.resdir}/data/GCSC/twas_weights/GTEx.{wildcards.gcsc_tissue}.P01/{wildcards.gcsc_tissue}.P01.pos \
    --weights_dir {params.resdir}/data/GCSC/twas_weights/GTEx.{wildcards.gcsc_tissue}.P01 \
    --ref_ld_chr {params.resdir}/data/1kg/1KG.Phase3.EUR.MAF_001.chr \
    --out {output} \
    --chr {wildcards.chr}) > {log} 2>&1"

rule twas_gcsc_all_chr:
    input:
      lambda w: expand("{outdir}/results/{gwas}/gcsc/twas/{gcsc_tissue}/{gwas}_twas_{gcsc_tissue}_chr{chr}.dat", gwas=w.gwas, gcsc_tissue=w.gcsc_tissue, chr=chromosomes, outdir={outdir})
    output:
      touch("{outdir}/results/{gwas}/checks/gcsc_twas_{gcsc_tissue}_all_chr.done")

rule twas_gcsc_all_panel:
    input:
      lambda w: expand("{outdir}/results/{gwas}/checks/gcsc_twas_{gcsc_tissue}_all_chr.done", gwas=w.gwas, gcsc_tissue=gcsc_tissues, outdir={outdir})
    output:
      touch("{outdir}/results/{gwas}/checks/gcsc_twas_all_panel.done")

####
# Prepare drug-gene interaction data
####

checkpoint prep_set_gcsc:
  input:
    "{outdir}/results/{gwas}/checks/gcsc_twas_all_panel.done",
    rules.download_drug_targetor.output,
    rules.install_gcsc.output,
    rules.download_biomart.output
  output:
    "{outdir}/results/{gwas}/gcsc/drugtargetor_gcsc_sets.nset.txt"
  benchmark:
    "{outdir}/benchmarks/prep_set_gcsc_{gwas}.tsv"
  conda:
    "../envs/main.yaml"
  params:
    config_file= config["config_file"]
  log:
    "{outdir}/logs/prep_set_gcsc-{gwas}.log"
  shell:
    "Rscript scripts/prep_set_gcsc.R \
      --gwas {wildcards.gwas} \
      --config_file {params.config_file} > {log} 2>&1"

def n_chunk_gcsc(x):
    checkpoint_output = checkpoints.prep_set_gcsc.get(gwas=x, outdir=outdir).output[0]
    checkpoint_output = outdir + "/results/" + x + "/gcsc/drugtargetor_gcsc_sets.nset.txt"
    n_chunk_gcsc_df = pd.read_table(checkpoint_output, sep=' ')
    return n_chunk_gcsc_df['x'].tolist()

####
# Run GCSC
####

# run gcsc with drugtargetor sets
rule run_gcsc_drugtargetor:
  resources:
    mem_mb=10000
  input:
    "{outdir}/results/{gwas}/gcsc/drugtargetor_gcsc_sets.nset.txt"
  output:
    "{outdir}/results/{gwas}/gcsc/drugtargetor/{chunk}/GCSCresults.txt"
  benchmark:
    "{outdir}/benchmarks/run_gcsc_drugtargetor_{gwas}_{chunk}.tsv"
  conda:
    "../envs/gcsc.yaml"
  params:
    gcsc_tissues= config["gcsc_tissues"],
    resdir=resdir
  log:
    "{outdir}/logs/run_gcsc_drugtargetor-{gwas}-{chunk}.log"
  shell:
    """
      (
      mkdir -p {wildcards.outdir}/results/{wildcards.gwas}/gcsc/drugtargetor/{wildcards.chunk};

      # Read the value and round to the nearest integer using printf
      RAW_N=$(cat {wildcards.outdir}/results/{wildcards.gwas}/gwas_sumstat/{wildcards.gwas}.cleaned.munged.median_N.txt);
      N=$(printf "%.0f" "$RAW_N");

      python {params.resdir}/software/GCSC/gcsc.py \
          --geneSets {wildcards.outdir}/results/{wildcards.gwas}/gcsc/drugtargetor_gcsc_sets_{wildcards.chunk}.csv \
          --TWASdir {wildcards.outdir}/results/{wildcards.gwas}/gcsc/twas/tissue \
          --N ${{N}} \
          --tissues {params.gcsc_tissues} \
          --coreg {params.resdir}/data/GCSC/coreg \
          --out {wildcards.outdir}/results/{wildcards.gwas}/gcsc/drugtargetor/{wildcards.chunk}
      ) > {log} 2>&1
    """

rule run_gcsc_all_chunk:
    input:
      lambda w: expand("{outdir}/results/{gwas}/gcsc/drugtargetor/{chunk}/GCSCresults.txt", gwas=w.gwas, chunk=n_chunk_gcsc("{}".format(w.gwas)), outdir={outdir})
    output:
      touch("{outdir}/results/{gwas}/checks/run_gcsc_all_chunk.done")

####
# Combine GCSC results
####

rule combine_gcsc:
  input:
    "{outdir}/results/{gwas}/checks/run_gcsc_all_chunk.done"
  output:
    "{outdir}/results/{gwas}/gcsc/{gwas}_drugtargetor_gcsc_res_atc.csv"
  benchmark:
    "{outdir}/benchmarks/combine_gcsc_{gwas}.tsv"
  conda:
    "../envs/main.yaml"
  params:
    config_file= config["config_file"]
  log:
    "{outdir}/logs/combine_gcsc-{gwas}.log"
  shell:
    "Rscript scripts/combine_gcsc.R \
      --gwas {wildcards.gwas} \
      --config_file {params.config_file} > {log} 2>&1"
