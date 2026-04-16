##################
# Run SuSiE
##################

def get_mem_mb_fine(wildcards, attempt):
    return attempt * 20000
    
rule finemap:
  resources:
    mem_mb=get_mem_mb_fine
  input:
    "{outdir}/results/{gwas}/gwas_sumstat/{gwas}.cleaned.gz",
    "{outdir}/results/{gwas}/clump/{gwas}.GW.clump.clean.csv"
  output:
    "{outdir}/results/{gwas}/checks/{gwas}.chr{chr}.finemap.done"
  benchmark:
    "{outdir}/benchmarks/finemap_{gwas}_chr{chr}.tsv"
  conda:
    "../envs/main.yaml"
  params:
    population= lambda w: gwas_list_df_eur.loc[gwas_list_df_eur['name'] == "{}".format(w.gwas), 'population'].iloc[0],
    config_file=config['config_file']
  log:
    "{outdir}/logs/finemap-{gwas}-chr{chr}.log"
  shell:
    "Rscript scripts/finemap.R \
      --gwas {wildcards.gwas} \
      --config_file {params.config_file} \
      --chr {wildcards.chr} > {log} 2>&1"

rule finemap_all_chr:
    input: 
      lambda w: expand("{outdir}/results/{gwas}/checks/{gwas}.chr{chr}.finemap.done", gwas=w.gwas, chr=chromosomes, outdir={outdir})
    output: 
      touch("{outdir}/results/{gwas}/checks/finemap_all_chr.done")

#################
# Process SuSiE results
#################

rule process_finemap:
  input:
    "{outdir}/results/{gwas}/checks/finemap_all_chr.done",
    rules.download_biomart.output
  output:
    "{outdir}/results/{gwas}/finemap/{gwas}.GW.finemap.L1.csv"
  benchmark:
    "{outdir}/benchmarks/process_finemap_{gwas}.tsv"
  conda:
    "../envs/main.yaml"
  params:
    config_file=config['config_file']
  log:
    "{outdir}/logs/process_finemap-{gwas}.log"
  shell:
    "Rscript scripts/process_finemap.R \
      --gwas {wildcards.gwas} \
      --config_file {params.config_file} > {log} 2>&1"


















