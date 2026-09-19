##########
# MAGMA gene-set x tissue interaction (post-hoc analysis of de Leeuw et al.
# 2018, Nat Commun 9:3768). One rule per GWAS; the R wrapper iterates every
# pathway GMT and every gene set >= min_set_size, invokes MAGMA once per set
# with condition-interaction=Average,{S}, applies BH FDR across the full
# (tissue x set) grid within the GWAS, and runs both spec-mandated follow-up
# tests on rows passing FDR. Reuses the existing magma_gene_level.genes.raw
# and GTEx_v8_tissue.tsv - no new gene analysis is performed here.
##########

_mi_n_cores = int(magma_interaction_cfg.get('n_cores', 1)) if magma_interaction_enabled else 1

rule magma_gene_set_by_tissue_interaction:
  input:
    genes_raw = "{outdir}/results/{gwas}/magma/magma_gene_level.genes.raw",
    tissue    = f"{resdir}/data/gtex/GTEx_v8_tissue.tsv",
    magma     = f"{resdir}/software/magma/magma",
    gmts      = [os.path.join(pathway_gmt_dir_val, f"{g}.gmt") for g in pathway_gmts] if pathway_gmts else []
  output:
    results     = "{outdir}/results/{gwas}/magma/interaction/interaction_results.tsv",
    calibration = "{outdir}/results/{gwas}/magma/interaction/interaction_calibration.tsv",
    commands    = "{outdir}/results/{gwas}/magma/interaction/interaction_commands.json"
  benchmark:
    "{outdir}/benchmarks/magma_gene_set_by_tissue_interaction_{gwas}.tsv"
  conda:
    "../envs/main.yaml"
  resources:
    cpus=_mi_n_cores
  params:
    resdir=resdir,
    config_file=config['config_file'],
    n_cores=_mi_n_cores
  log:
    "{outdir}/logs/magma_gene_set_by_tissue_interaction-{gwas}.log"
  shell:
    "Rscript --vanilla {workflow.basedir}/scripts/magma_gene_set_by_tissue_interaction.R \
       --pipeline_dir {workflow.basedir} \
       --gwas {wildcards.gwas} \
       --config_file {params.config_file} \
       --resdir {params.resdir} \
       --outdir {outdir} \
       --n_cores {params.n_cores} > {log} 2>&1"

rule magma_gene_set_by_tissue_interaction_all:
  input:
    expand("{outdir}/results/{gwas}/magma/interaction/interaction_results.tsv",
           gwas=gwas_list_df_eur['name'], outdir=[outdir])
