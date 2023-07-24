####
# Download GCSC
####

rule install_gcsc:
  output:
    directory("resources/software/GCSC/")
  conda:
    "../envs/GenoFunc.yaml"
  shell:
    "git clone git@github.com:ksiewert/GCSC.git {output}"

####
# Download GCSC gene co-regulation scores
####

gcsc_tissues=config["gcsc_tissues"]

rule download_gcsc_coreg:
  output:
    "resources/data/GCSC/coreg/{gcsc_tissue}_geneNames.txt"
  conda:
    "../envs/GenoFunc.yaml"
  shell:
    "wget -O resources/data/GCSC/coreg/{wildcards.gcsc_tissue}_coregscores.npz https://storage.googleapis.com/broad-alkesgroup-public/GCSC/Coreg_scores/{wildcards.gcsc_tissue}_coregscores.npz; wget -O resources/data/GCSC/coreg/{wildcards.gcsc_tissue}_geneNames.txt https://storage.googleapis.com/broad-alkesgroup-public/GCSC/Coreg_scores/{wildcards.gcsc_tissue}_geneNames.txt"

rule download_gcsc_coreg_all_tissue:
    input: expand("resources/data/GCSC/coreg/{gcsc_tissue}_geneNames.txt", gcsc_tissue=gcsc_tissues)

####
# Download corresponding GTEx v7 TWAS weights
####

gcsc_tissues=config["gcsc_tissues"]

rule download_gcsc_twas_weights:
  output:
    directory("resources/data/GCSC/twas_weights/GTEx.{gcsc_tissue}.P01")
  conda:
    "../envs/GenoFunc.yaml"
  shell:
    "mkdir resources/data/GCSC/twas_weights/GTEx.{wildcards.gcsc_tissue}.P01; wget -O resources/data/GCSC/twas_weights/GTEx.{wildcards.gcsc_tissue}.P01.tar.bz2 http://gusevlab.org/projects/fusion/weights/GTEx.{wildcards.gcsc_tissue}.P01.tar.bz2; tar xjvf resources/data/GCSC/twas_weights/GTEx.{wildcards.gcsc_tissue}.P01.tar.bz2 -C resources/data/GCSC/twas_weights/GTEx.{wildcards.gcsc_tissue}.P01; rm resources/data/GCSC/twas_weights/GTEx.{wildcards.gcsc_tissue}.P01.tar.bz2"

rule download_gcsc_twas_weights_all_tissue:
    input: expand("resources/data/GCSC/twas_weights/GTEx.{gcsc_tissue}.P01", gcsc_tissue=gcsc_tissues)

####
# Perform TWAS using GTEx v7 weights
####

# run twas
rule run_twas_gcsc:
  resources:
    mem_mb=20000
  input:
    "resources/data/gwas_sumstat/{gwas}/{gwas}.cleaned.munged.sumstats.gz",
    "resources/data/gwas_sumstat/{gwas}/{gwas}.cleaned.munged.median_N.txt",
    rules.install_fusion.output,
    rules.install_plink2R.output,
    rules.prep_1kg.output,
    "resources/data/GCSC/coreg/{gcsc_tissue}_geneNames.txt",
    "resources/data/GCSC/twas_weights/GTEx.{gcsc_tissue}.P01"
  output:
    "results/{gwas}/gcsc/twas/{gcsc_tissue}/{gwas}_twas_{gcsc_tissue}_chr{chr}"
  conda:
    "../envs/GenoFunc.yaml"
  shell:
    "mkdir -p results/{wildcards.gwas}/gcsc/twas/{wildcards.gcsc_tissue}; N=$(cat resources/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}.cleaned.munged.median_N.txt); Rscript resources/software/fusion/FUSION.assoc_test.R \
    --sumstats resources/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}.cleaned.munged.sumstats.gz \
    --weights resources/data/GCSC/twas_weights/GTEx.{wildcards.gcsc_tissue}.P01/{wildcards.gcsc_tissue}.P01.pos \
    --weights_dir resources/data/GCSC/twas_weights/GTEx.{wildcards.gcsc_tissue}.P01 \
    --ref_ld_chr resources/data/1kg/1KG.Phase3.EUR.MAF_001.chr \
    --out {output}.dat \
    --chr {wildcards.chr}"

rule twas_gcsc_all_chr:
    input: 
      lambda w: expand("results/{gwas}/gcsc/twas/{gcsc_tissue}/{gwas}_twas_{gcsc_tissue}_chr{chr}", gwas=w.gwas, gcsc_tissue=w.gcsc_tissue, chr=range(1, 23))
    output: 
      touch("results/{gwas}/checks/gcsc_twas_{gcsc_tissue}_all_chr.done")

rule twas_gcsc_all_panel:
    input: 
      lambda w: expand("results/{gwas}/checks/gcsc_twas_{gcsc_tissue}_all_chr.done", gwas=w.gwas, gcsc_tissue=gcsc_tissues)
    output: 
      touch("results/{gwas}/checks/gcsc_twas_all_panel.done")
      
####
# Prepare drug-gene interaction data
####

checkpoint prep_set_gcsc:
  input:
    "results/{gwas}/checks/gcsc_twas_all_panel.done",
    rules.download_drug_targetor.output
  output:
    "results/{gwas}/gcsc/drugtargetor_gcsc_sets.nset.txt"
  conda: 
    "../envs/GenoFunc.yaml"
  params:
    config_file= config["config_file"]
  shell:
    "Rscript scripts/prep_set_gcsc.R \
      --gwas {wildcards.gwas} \
      --config {params.config_file}"

def n_chunk_gcsc(x):
    checkpoint_output = checkpoints.prep_set_gcsc.get(gwas=x).output[0]
    checkpoint_output = "results/" + x + "/gcsc/drugtargetor_gcsc_sets.nset.txt"
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
    "results/{gwas}/gcsc/drugtargetor_gcsc_sets.nset.txt"
  output:
    "results/{gwas}/gcsc/drugtargetor/{chunk}/GCSCresults.txt"
  conda:
    "../envs/gcsc.yaml"
  params:
    gcsc_tissues= config["gcsc_tissues"]
  shell:
    "mkdir -p results/ALS_only/gcsc/drugtargetor/{wildcards.chunk}; N=$(cat resources/data/gwas_sumstat/{wildcards.gwas}/{wildcards.gwas}.cleaned.munged.median_N.txt); python resources/software/GCSC/gcsc.py \
--geneSets results/{wildcards.gwas}/gcsc/drugtargetor_gcsc_sets_{wildcards.chunk}.csv \
--TWASdir results/{wildcards.gwas}/gcsc/twas/tissue \
--N ${{N}} \
--tissues {params.gcsc_tissues} \
--coreg resources/data/GCSC/coreg \
--out results/{wildcards.gwas}/gcsc/drugtargetor/{wildcards.chunk}"

rule run_gcsc_all_chunk:
    input: 
      lambda w: expand("results/{gwas}/gcsc/drugtargetor/{chunk}/GCSCresults.txt", gwas=w.gwas, chunk=n_chunk_gcsc("{}".format(w.gwas)))
    output: 
      touch("results/{gwas}/checks/run_gcsc_all_chunk.done")
      
####
# Combine GCSC results
####

rule combine_gcsc:
  input:
    "results/{gwas}/checks/run_gcsc_all_chunk.done"
  output:
    "results/{gwas}/gcsc/{gwas}_drugtargetor_gcsc_res_atc.txt"
  conda: 
    "../envs/GenoFunc.yaml"
  params:
    config_file= config["config_file"]
  shell:
    "Rscript scripts/combine_gcsc.R \
      --gwas {wildcards.gwas}"


