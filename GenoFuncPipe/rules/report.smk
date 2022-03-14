####
# Create results report
####
# Create a report for each GWAS

myoutput = list()

if config["ldsc"] == "T":
    myoutput.append("results/{gwas}/ldsc/{gwas}_ldsc_res.log")
    
if config["magma_gene"] == "T":
    myoutput.append("results/{gwas}/magma/magma_gene_level.genes.raw")
    
if config["magma_drugtargetor"] == "T":
    myoutput.append("results/{gwas}/magma/magma_drug_targetor_atc_res.csv")

if config["twas_panel_psychencode"] == "T":
    myoutput.append("results/{gwas}/checks/psychencode_twas_all_chr.done")

if config["twas_gsea_lincs"] == "T":
    myoutput.append("results/{gwas}/twas/cmap/twas_gsea_res_atc_res.csv")

if config["twas_so_lincs"] == "T":
    myoutput.append("results/{gwas}/twas/cmap/So_res_atc_res.csv")

if config["pwas_panel_rosmap"] == "T":
    myoutput.append("results/{gwas}/checks/rosmap_pwas_all_chr.done")

if config["pwas_panel_banner"] == "T":
    myoutput.append("results/{gwas}/checks/banner_pwas_all_chr.done")

if config["smr_expression_panel_psychencode"] == "T":
    myoutput.append("results/{gwas}/checks/psychencode_smr_all_chr.done")

if config["smr_panel_metabrain_basalganglia"] == "T":
    myoutput.append("results/{gwas}/checks/metabrain_smr_basalganglia_all_chr.done")

if config["smr_panel_metabrain_cerebellum"] == "T":
    myoutput.append("results/{gwas}/checks/metabrain_smr_cerebellum_all_chr.done")

if config["smr_panel_metabrain_cortex"] == "T":
    myoutput.append("results/{gwas}/checks/metabrain_smr_cortex_all_chr.done")

if config["smr_panel_metabrain_hippocampus"] == "T":
    myoutput.append("results/{gwas}/checks/metabrain_smr_hippocampus_all_chr.done")

if config["smr_panel_metabrain_spinalcord"] == "T":
    myoutput.append("results/{gwas}/checks/metabrain_smr_spinalcord_all_chr.done")

if config["smr_protein_panel_rosmap"] == "T":
    myoutput.append("results/{gwas}/checks/rosmap_smr_all_chr.done")

rule create_report:
  input:
    "results/{gwas}/cojo/{gwas}.GW.cojo.clean.csv",
    myoutput
  output:
    "results/{gwas}/reports/{gwas}_report.html"
  conda:
    "../envs/GenoFunc.yaml"
  shell:
    "mkdir -p results/{wildcards.gwas}/reports; Rscript -e \"rmarkdown::render(\'scripts/create_report.Rmd\', \
     output_file = \'../results/{wildcards.gwas}/reports/{wildcards.gwas}_report.html\', \
     params = list(gwas = \'{wildcards.gwas}\'))\""

rule run_create_report:
  input: expand("results/{gwas}/reports/{gwas}_report.html", gwas=gwas_list_df_eur['name'])

