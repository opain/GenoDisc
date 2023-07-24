####
# Create results report
####
# Create a report for each GWAS

rule install_poolr:
  output:
    touch("resources/software/install_poolr.done")
  conda:
    "../envs/GenoFunc.yaml"
  shell:
    "Rscript -e 'install.packages(\"poolr\", repos = \"http://cran.us.r-project.org\")'"
    
myoutput = list()

if config["clump"] == "T":
    myoutput.append("results/{gwas}/clump/{gwas}.GW.clump.clean.csv")

if config["cojo"] == "T":
    myoutput.append("results/{gwas}/cojo/{gwas}.GW.cojo.clean.csv")

if config["finemap"] == "T":
    myoutput.append("results/{gwas}/finemap/{gwas}.GW.finemap.L1.csv")

if config["ldsc"] == "T":
    myoutput.append("results/{gwas}/ldsc/{gwas}_ldsc_res.log")
    
if config["magma_gene"] == "T":
    myoutput.append("results/{gwas}/magma/magma_gene_level.clean.csv")
    
if config["magma_drugtargetor"] == "T":
    myoutput.append("results/{gwas}/magma/magma_drug_targetor_atc_res.csv")

if config["twas_panel_fusion"] == "T" or config["twas_panel_psychencode"] == "T":
    myoutput.append("results/{gwas}/twas/{gwas}_twas_GW_clean.txt.gz")

if config["twas_conditional"] == "T" and (config["twas_panel_fusion"] == "T" or config["twas_panel_psychencode"] == "T"):
    myoutput.append("results/{gwas}/twas/{gwas}_twas_novelty.csv")

if config["twas_gsea_lincs"] == "T":
    myoutput.append("results/{gwas}/checks/format_twas_gsea_results_all_panel.done")

if config["twas_so_lincs"] == "T":
    myoutput.append("results/{gwas}/checks/format_so_results_all_panel.done")

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

if config["smr_expression_panel_eqtlgen"] == "T":
    myoutput.append("results/{gwas}/smr/eqtlgen/{gwas}_smr_eqtlgen_GW.txt.gz")

if config["smr_protein_panel_rosmap"] == "T":
    myoutput.append("results/{gwas}/smr/rosmap/{gwas}_smr_rosmap_GW.txt.gz")

if config["dgi_db_comp"] == "T":
    myoutput.append("results/{gwas}/DGIdb/DGIdb_opposing_clean.csv")

if config["twas_gsea_drugtargetor"] == "T":
    myoutput.append("results/{gwas}/checks/format_twas_gsea_drugtargetor_results_all_panel.done")
    
rule create_report:
  input:
    myoutput,
    "resources/software/install_poolr.done"
  output:
    "results/{gwas}/reports/{gwas}_report.html"
  conda:
    "../envs/GenoFunc.yaml"
  params:
    config_file= config["config_file"]
  shell:
    "mkdir -p results/{wildcards.gwas}/reports; Rscript -e \"rmarkdown::render(\'scripts/create_report.Rmd\', \
     output_file = \'../results/{wildcards.gwas}/reports/{wildcards.gwas}_report.html\', \
     params = list(gwas = \'{wildcards.gwas}\', config_file = \'{params.config_file}\'))\""

rule run_create_report:
  input: expand("results/{gwas}/reports/{gwas}_report.html", gwas=gwas_list_df_eur['name'])

