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
# 'auto' -> R script prefers /dev/shm if writable, else R's tempdir(); on shared-FS
# clusters (CephFS/Lustre) this is ~9x faster than the default when running with
# n_cores > 1. Users can override to any node-local path (e.g. /dev/shm or a
# local scratch SSD) via magma_interaction.tmpdir in the config.
_mi_tmpdir = (magma_interaction_cfg.get('tmpdir', 'auto') if magma_interaction_enabled else 'auto')

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
  # `threads:` is what slurm profiles universally translate into
  # `--cpus-per-task`, so declaring it here makes sure each job actually
  # gets the requested cores. `resources.cpus` is kept alongside as
  # supplementary metadata for profiles that read it (or for callers
  # setting it via --set-resources).
  threads: _mi_n_cores
  resources:
    cpus=_mi_n_cores
  params:
    resdir=resdir,
    config_file=config['config_file'],
    n_cores=_mi_n_cores,
    tmpdir=_mi_tmpdir
  log:
    "{outdir}/logs/magma_gene_set_by_tissue_interaction-{gwas}.log"
  shell:
    "Rscript --vanilla {workflow.basedir}/scripts/magma_gene_set_by_tissue_interaction.R \
       --pipeline_dir {workflow.basedir} \
       --gwas {wildcards.gwas} \
       --config_file {params.config_file} \
       --resdir {params.resdir} \
       --outdir {outdir} \
       --n_cores {params.n_cores} \
       --tmpdir {params.tmpdir} > {log} 2>&1"

rule magma_gene_set_by_tissue_interaction_all:
  input:
    expand("{outdir}/results/{gwas}/magma/interaction/interaction_results.tsv",
           gwas=gwas_list_df_eur['name'], outdir=[outdir])
