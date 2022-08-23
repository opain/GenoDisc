##################
# Query DGIdb
##################

myoutput = list()

if config["twas_panel_fusion"] == "T" or config["twas_panel_psychencode"] == "T":
    myoutput.append("results/{gwas}/twas/{gwas}_twas_GW_clean.txt.gz")

if config["pwas_panel_rosmap"] == "T":
    myoutput.append("results/{gwas}/checks/rosmap_pwas_all_chr.done")

if config["pwas_panel_banner"] == "T":
    myoutput.append("results/{gwas}/checks/banner_pwas_all_chr.done")

if config["smr_expression_panel_psychencode"] == "T":
    myoutput.append("results/{gwas}/checks/psychencode_smr_all_chr.done")

if config["smr_expression_panel_metabrain_basalganglia"] == "T":
    myoutput.append("results/{gwas}/smr/metabrain/{gwas}_smr_metabrain_GW.txt.gz")

if config["smr_expression_panel_metabrain_cerebellum"] == "T":
    myoutput.append("results/{gwas}/smr/metabrain/{gwas}_smr_metabrain_GW.txt.gz")

if config["smr_expression_panel_metabrain_cortex"] == "T":
    myoutput.append("results/{gwas}/smr/metabrain/{gwas}_smr_metabrain_GW.txt.gz")

if config["smr_expression_panel_metabrain_hippocampus"] == "T":
    myoutput.append("results/{gwas}/smr/metabrain/{gwas}_smr_metabrain_GW.txt.gz")

if config["smr_expression_panel_metabrain_spinalcord"] == "T":
    myoutput.append("results/{gwas}/smr/metabrain/{gwas}_smr_metabrain_GW.txt.gz")

if config["smr_protein_panel_rosmap"] == "T":
    myoutput.append("results/{gwas}/smr/rosmap/{gwas}_smr_rosmap_GW.txt.gz")

rule query_dgidb:
  input:
    myoutput
  output:
    "results/{gwas}/DGIdb/DGIdb_opposing_clean.csv"
  conda:
    "../envs/GenoFunc.yaml"
  params:
    config_file=config["config_file"]
  shell:
    "Rscript scripts/DGIdb_query.R \
      --gwas {wildcards.gwas} \
      --config_file {params.config_file}"















