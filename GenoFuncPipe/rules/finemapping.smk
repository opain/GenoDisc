##################
# Run SuSiE
##################

rule finemap:
  input:
    "resources/data/gwas_sumstat/{gwas}/{gwas}.cleaned.gz",
    "results/{gwas}/clump/{gwas}.GW.clump.clean.csv"
  output:
    "results/{gwas}/checks/{gwas}.chr{chr}.finemap.done"
  conda:
    "../envs/GenoFunc.yaml"
  params:
    population= lambda w: gwas_list_df_eur.loc[gwas_list_df_eur['name'] == "{}".format(w.gwas), 'population'].iloc[0]
  shell:
    "Rscript scripts/finemap.R \
      --gwas {wildcards.gwas} \
      --chr {wildcards.chr}"

rule finemap_all_chr:
    input: 
      lambda w: expand("results/{gwas}/checks/{gwas}.chr{chr}.finemap.done", gwas=w.gwas, chr=range(1, 23))
    output: 
      touch("results/{gwas}/checks/finemap_all_chr.done")

#################
# Process SuSiE results
#################

rule process_finemap:
  input:
    "results/{gwas}/checks/finemap_all_chr.done"
  output:
    "results/{gwas}/finemap/{gwas}.GW.finemap.L1.csv"
  conda:
    "../envs/GenoFunc.yaml"
  shell:
    "Rscript scripts/process_finemap.R \
      --gwas {wildcards.gwas}"


















