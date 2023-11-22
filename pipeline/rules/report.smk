####
# Create results package
####
    
myoutput = list()

if config["clump"] == "T":
    myoutput.append(expand("{outdir}/results/{gwas}/clump/{gwas}.GW.clump.clean.csv", gwas=gwas_list_df_eur['name'], outdir={outdir}))

if config["cojo"] == "T":
    myoutput.append(expand("{outdir}/results/{gwas}/cojo/{gwas}.GW.cojo.clean.csv", gwas=gwas_list_df_eur['name'], outdir={outdir}))

if config["finemap"] == "T":
    myoutput.append(expand("{outdir}/results/{gwas}/finemap/{gwas}.GW.finemap.L1.csv", gwas=gwas_list_df_eur['name'], outdir={outdir}))

if config["ldsc"] == "T":
    myoutput.append(expand("{outdir}/results/{gwas}/ldsc/{gwas}_ldsc_res.log", gwas=gwas_list_df_eur['name'], outdir={outdir}))
    
if config["magma_gene"] == "T":
    myoutput.append(expand("{outdir}/results/{gwas}/magma/magma_gene_level.clean.csv", gwas=gwas_list_df_eur['name'], outdir={outdir}))
    
if config["magma_drugtargetor"] == "T":
    myoutput.append(expand("{outdir}/results/{gwas}/magma/magma_drug_targetor_atc_res.csv", gwas=gwas_list_df_eur['name'], outdir={outdir}))

if config["twas_panel_fusion"] == "T" or config["twas_panel_psychencode"] == "T":
    myoutput.append(expand("{outdir}/results/{gwas}/twas/{gwas}_twas_GW_clean.txt.gz", gwas=gwas_list_df_eur['name'], outdir={outdir}))

if config["twas_conditional"] == "T" and (config["twas_panel_fusion"] == "T" or config["twas_panel_psychencode"] == "T"):
    myoutput.append(expand("{outdir}/results/{gwas}/twas/{gwas}_twas_novelty.csv", gwas=gwas_list_df_eur['name'], outdir={outdir}))

if config["pwas_panel_rosmap"] == "T":
    myoutput.append(expand("{outdir}/results/{gwas}/checks/rosmap_pwas_all_chr.done", gwas=gwas_list_df_eur['name'], outdir={outdir}))

if config["pwas_panel_banner"] == "T":
    myoutput.append(expand("{outdir}/results/{gwas}/checks/banner_pwas_all_chr.done", gwas=gwas_list_df_eur['name'], outdir={outdir}))

if config["smr_expression_panel_psychencode"] == "T":
    myoutput.append(expand("{outdir}/results/{gwas}/checks/psychencode_smr_all_chr.done", gwas=gwas_list_df_eur['name'], outdir={outdir}))

if config["smr_expression_panel_metabrain_basalganglia"] == "T":
    myoutput.append(expand("{outdir}/results/{gwas}/smr/metabrain/{gwas}_smr_metabrain_GW.txt.gz", gwas=gwas_list_df_eur['name'], outdir={outdir}))

if config["smr_expression_panel_metabrain_cerebellum"] == "T":
    myoutput.append(expand("{outdir}/results/{gwas}/smr/metabrain/{gwas}_smr_metabrain_GW.txt.gz", gwas=gwas_list_df_eur['name'], outdir={outdir}))

if config["smr_expression_panel_metabrain_cortex"] == "T":
    myoutput.append(expand("{outdir}/results/{gwas}/smr/metabrain/{gwas}_smr_metabrain_GW.txt.gz", gwas=gwas_list_df_eur['name'], outdir={outdir}))

if config["smr_expression_panel_metabrain_hippocampus"] == "T":
    myoutput.append(expand("{outdir}/results/{gwas}/smr/metabrain/{gwas}_smr_metabrain_GW.txt.gz", gwas=gwas_list_df_eur['name'], outdir={outdir}))

if config["smr_expression_panel_metabrain_spinalcord"] == "T":
    myoutput.append(expand("{outdir}/results/{gwas}/smr/metabrain/{gwas}_smr_metabrain_GW.txt.gz", gwas=gwas_list_df_eur['name'], outdir={outdir}))

if config["smr_expression_panel_eqtlgen"] == "T":
    myoutput.append(expand("{outdir}/results/{gwas}/smr/eqtlgen/{gwas}_smr_eqtlgen_GW.txt.gz", gwas=gwas_list_df_eur['name'], outdir={outdir}))

if config["smr_protein_panel_rosmap"] == "T":
    myoutput.append(expand("{outdir}/results/{gwas}/smr/rosmap/{gwas}_smr_rosmap_GW.txt.gz", gwas=gwas_list_df_eur['name'], outdir={outdir}))

if config["twas_gsea_drugtargetor"] == "T":
    myoutput.append(expand("{outdir}/results/{gwas}/checks/format_twas_gsea_drugtargetor_results_all_panel.done", gwas=gwas_list_df_eur['name'], outdir={outdir}))

if config["gcsc"] == "T":
    myoutput.append(expand("{outdir}/results/{gwas}/gcsc/{gwas}_drugtargetor_gcsc_res_atc.csv", gwas=gwas_list_df_eur['name'], outdir={outdir}))

rule package_results:
  input:
    myoutput
  output:
    "{outdir}/results/results_package.rds"
  conda:
    "../envs/main.yaml"
  params:
    config_file= config["config_file"]
  shell:
    "Rscript scripts/package_results.R \
    --config {params.config_file}"

rule run_package_results:
  input: expand("{outdir}/results/results_package.rds", outdir={outdir})