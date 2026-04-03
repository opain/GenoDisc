####
# Download GCSC
####

rule install_gcsc:
  output:
    directory(f"{resdir}/software/GCSC/")
  conda:
    "../envs/main.yaml"
  shell:
    "git clone https://github.com/ksiewert/GCSC.git {output}; \
     cd {output}; \
     git reset --hard b10ea77b9a43399801b46ef70c80516599264123"

####
# Download GCSC gene co-regulation scores
####

gcsc_tissues=config["gcsc_tissues"]

rule download_gcsc_coreg:
  output:
    f"{resdir}/data/GCSC/coreg/{{gcsc_tissue}}_geneNames.txt"
  conda:
    "../envs/main.yaml"
  params:
    resdir=resdir
  shell:
    "wget -O {params.resdir}/data/GCSC/coreg/{wildcards.gcsc_tissue}_coregscores.npz https://storage.googleapis.com/broad-alkesgroup-public/GCSC/Coreg_scores/{wildcards.gcsc_tissue}_coregscores.npz; \
    wget -O {params.resdir}/data/GCSC/coreg/{wildcards.gcsc_tissue}_geneNames.txt https://storage.googleapis.com/broad-alkesgroup-public/GCSC/Coreg_scores/{wildcards.gcsc_tissue}_geneNames.txt"

rule download_gcsc_coreg_all_tissue:
    input: expand(f"{resdir}/data/GCSC/coreg/{{gcsc_tissue}}_geneNames.txt", gcsc_tissue=gcsc_tissues)

####
# Download corresponding GTEx v7 TWAS weights
####

gcsc_tissues=config["gcsc_tissues"]

rule download_gcsc_twas_weights:
  output:
    directory(f"{resdir}/data/GCSC/twas_weights/GTEx.{{gcsc_tissue}}.P01")
  conda:
    "../envs/main.yaml"
  params:
    resdir=resdir
  shell:
    "mkdir {params.resdir}/data/GCSC/twas_weights/GTEx.{wildcards.gcsc_tissue}.P01; \
    wget -O {params.resdir}/data/GCSC/twas_weights/GTEx.{wildcards.gcsc_tissue}.P01.tar.bz2 http://gusevlab.org/projects/fusion/weights/GTEx.{wildcards.gcsc_tissue}.P01.tar.bz2; \
    tar xjvf {params.resdir}/data/GCSC/twas_weights/GTEx.{wildcards.gcsc_tissue}.P01.tar.bz2 -C {params.resdir}/data/GCSC/twas_weights/GTEx.{wildcards.gcsc_tissue}.P01; \
    rm {params.resdir}/data/GCSC/twas_weights/GTEx.{wildcards.gcsc_tissue}.P01.tar.bz2"

rule download_gcsc_twas_weights_all_tissue:
    input: expand(f"{resdir}/data/GCSC/twas_weights/GTEx.{{gcsc_tissue}}.P01", gcsc_tissue=gcsc_tissues)

####
# Perform TWAS using GTEx v7 weights
####

# run twas
rule run_twas_gcsc:
  resources:
    mem_mb=20000
  input:
    "{outdir}/data/gwas_sumstat/{gwas}/{gwas}.cleaned.munged.sumstats.gz",
    "{outdir}/data/gwas_sumstat/{gwas}/{gwas}.cleaned.munged.median_N.txt",
    rules.install_fusion.output,
    rules.install_plink2R.output,
    rules.prep_1kg.output,
    f"{resdir}/data/GCSC/coreg/{{gcsc_tissue}}_geneNames.txt",
    f"{resdir}/data/GCSC/twas_weights/GTEx.{{gcsc_tissue}}.P01"
  output:
    "{outdir}/results/{gwas}/gcsc/twas/{gcsc_tissue}/{gwas}_twas_{gcsc_tissue}_chr{chr}.dat"
  conda:
    "../envs/main.yaml"
  params:
    resdir=resdir
  shell:
    "mkdir -p {outdir}/results/{wildcards.gwas}/gcsc/twas/{wildcards.gcsc_tissue}; N=$(cat {outdir}/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}.cleaned.munged.median_N.txt); Rscript {params.resdir}/software/fusion/FUSION.assoc_test.R \
    --sumstats {outdir}/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}.cleaned.munged.sumstats.gz \
    --weights {params.resdir}/data/GCSC/twas_weights/GTEx.{wildcards.gcsc_tissue}.P01/{wildcards.gcsc_tissue}.P01.pos \
    --weights_dir {params.resdir}/data/GCSC/twas_weights/GTEx.{wildcards.gcsc_tissue}.P01 \
    --ref_ld_chr {params.resdir}/data/1kg/1KG.Phase3.EUR.MAF_001.chr \
    --out {output} \
    --chr {wildcards.chr}"

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
  conda: 
    "../envs/main.yaml"
  params:
    config_file= config["config_file"]
  shell:
    "Rscript scripts/prep_set_gcsc.R \
      --gwas {wildcards.gwas} \
      --config_file {params.config_file}"

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
  conda:
    "../envs/gcsc.yaml"
  params:
    gcsc_tissues= config["gcsc_tissues"],
    resdir=resdir
  shell:
    "mkdir -p {outdir}/results/{wildcards.gwas}/gcsc/drugtargetor/{wildcards.chunk}; N=$(cat {outdir}/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}.cleaned.munged.median_N.txt); python {params.resdir}/software/GCSC/gcsc.py \
--geneSets {outdir}/results/{wildcards.gwas}/gcsc/drugtargetor_gcsc_sets_{wildcards.chunk}.csv \
--TWASdir {outdir}/results/{wildcards.gwas}/gcsc/twas/tissue \
--N ${{N}} \
--tissues {params.gcsc_tissues} \
--coreg {params.resdir}/data/GCSC/coreg \
--out {outdir}/results/{wildcards.gwas}/gcsc/drugtargetor/{wildcards.chunk}"

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
  conda: 
    "../envs/main.yaml"
  params:
    config_file= config["config_file"]
  shell:
    "Rscript scripts/combine_gcsc.R \
      --gwas {wildcards.gwas} \
      --config_file {params.config_file}"


