##################
# Run SBayesRC via the SBayesRC R package (SNP-h2 + polygenicity)
##################
#
# ALTERNATIVE SBayesRC implementation (flag: sbayesrc_r). Uses the SBayesRC R
# package (zhilizheng, v0.2.6) tidy -> impute -> sbayesrc, multi-threaded via
# OpenMP over a precomputed eigen-LD panel (the GenoPred route) -- typically
# faster than the from-scratch gctb-binary method in rules/sbayesrc.smk. Keeps
# only the model-level parameters (SNP-h2 + polygenicity); no PGS weights.
#
# References: the R package consumes the SBayesRC-author eigen-LD panel (shared
# with GenoPred) + the baseline2.2 annotation. config.yaml sbayesrc_r_ldref /
# sbayesrc_r_annot default to NA = use the copies staged in resdir
# (download_sbayesrc_r_ref + the shared sbayesrc_annot); set an explicit path to
# reuse an existing copy (e.g. GenoPred's staged refs on CREATE).
_sbrc_r_ldref_cfg = config.get("sbayesrc_r_ldref", "NA")
sbayesrc_r_ldref_base = (f"{resdir}/data/sbayesrc_r_ref"
                         if _sbrc_r_ldref_cfg in (None, "NA") else _sbrc_r_ldref_cfg)

_sbrc_r_annot_cfg = config.get("sbayesrc_r_annot", "NA")
sbayesrc_r_annot = (f"{resdir}/data/sbayesrc_annot/annot_baseline2.2.txt"
                    if _sbrc_r_annot_cfg in (None, "NA") else _sbrc_r_annot_cfg)

# The R package parallelises via OpenMP (OMP_NUM_THREADS); 10 matches GenoPred's
# default. threads / resources.cpus must agree (the SLURM profile reads cpus).
sbayesrc_r_threads = int(config.get("sbayesrc_r_threads", 10))
sbayesrc_r_seed    = int(config.get("sbayesrc_r_seed", 1))

def _sbrc_r_pop(w):
    return gwas_list_df.loc[gwas_list_df['name'] == "{}".format(w.gwas), 'population'].iloc[0]

def get_mem_mb_sbayesrc_r(wildcards, attempt):
    return attempt * 20000

rule sbayesrc_r:
  resources:
    mem_mb=get_mem_mb_sbayesrc_r,
    cpus=sbayesrc_r_threads
  threads: sbayesrc_r_threads
  input:
    sumstats="{outdir}/results/{gwas}/gwas_sumstat/{gwas}.cleaned.gz",
    sbrc_pkg=rules.install_sbayesrc_rpkg.output,
    ldref=lambda w: f"{sbayesrc_r_ldref_base}/{_sbrc_r_pop(w)}/ldm.info",
    annot=sbayesrc_r_annot
  output:
    "{outdir}/results/{gwas}/sbayesrc_r/{gwas}.GW.sbayesrc_r.clean.csv"
  benchmark:
    "{outdir}/benchmarks/sbayesrc_r_{gwas}.tsv"
  conda:
    "../envs/sbayesrc_r.yaml"
  params:
    ldref=lambda w: f"{sbayesrc_r_ldref_base}/{_sbrc_r_pop(w)}",
    annot=sbayesrc_r_annot,
    seed=sbayesrc_r_seed,
    config_file=config['config_file']
  log:
    "{outdir}/logs/sbayesrc_r-{gwas}.log"
  shell:
    "export OMP_NUM_THREADS={resources.cpus}; \
     Rscript --vanilla {workflow.basedir}/scripts/sbayesrc_r.R --pipeline_dir {workflow.basedir} \
      --gwas {wildcards.gwas} \
      --config_file {params.config_file} \
      --ldref {params.ldref} \
      --annot {params.annot} \
      --seed {params.seed} \
      --threads {resources.cpus} > {log} 2>&1"
