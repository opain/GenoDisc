##################
# Run GCTB SBayesRC (SNP-h2 + polygenicity)
##################
#
# Genome-wide SBayesRC (one run per GWAS, not per chromosome). We keep only the
# model-level parameters (SNP-h2 + polygenicity); the per-SNP effects (.snpRes)
# are not retained. The LD reference is chosen from the per-GWAS population
# column (EUR/EAS/AFR; see gwas_list_df_sbrc in sumstat_qc.smk).
#
# Reference data (eigen LD + baseline annotation) is staged into resdir by
# download_sbayesrc_ref / download_sbayesrc_annot (dependencies.smk) from the
# official GCTB site. config.yaml sbayesrc_ldref / sbayesrc_annot default to NA
# = use those staged copies; set an explicit path to point elsewhere.
_sbrc_ldref_cfg = config.get("sbayesrc_ldref", "NA")
sbayesrc_ldref_base = (f"{resdir}/data/sbayesrc_ref"
                       if _sbrc_ldref_cfg in (None, "NA") else _sbrc_ldref_cfg)

_sbrc_annot_cfg = config.get("sbayesrc_annot", "NA")
sbayesrc_annot = (f"{resdir}/data/sbayesrc_annot/annot_baseline2.2.txt"
                  if _sbrc_annot_cfg in (None, "NA") else _sbrc_annot_cfg)

# MCMC controls (explicit + reproducible); overridable per deployment/run.
sbayesrc_threads      = int(config.get("sbayesrc_threads", 4))
sbayesrc_seed         = int(config.get("sbayesrc_seed", 1))
sbayesrc_chain_length = int(config.get("sbayesrc_chain_length", 3000))
sbayesrc_burn_in      = int(config.get("sbayesrc_burn_in", 1000))

def _sbrc_pop(w):
    return gwas_list_df.loc[gwas_list_df['name'] == "{}".format(w.gwas), 'population'].iloc[0]

# Observed peak RSS for a whole-genome run is ~9.5 GB (EUR); 24 GB gives ample
# headroom for larger refs (AFR) and retries.
def get_mem_mb_sbayesrc(wildcards, attempt):
    return attempt * 24000

rule sbayesrc:
  resources:
    mem_mb=get_mem_mb_sbayesrc
  threads: sbayesrc_threads
  input:
    sumstats="{outdir}/results/{gwas}/gwas_sumstat/{gwas}.cleaned.gz",
    gctb=rules.download_gctb.output,
    ldref=lambda w: f"{sbayesrc_ldref_base}/{_sbrc_pop(w)}/ldm.info",
    annot=sbayesrc_annot
  output:
    "{outdir}/results/{gwas}/sbayesrc/{gwas}.GW.sbayesrc.clean.csv"
  benchmark:
    "{outdir}/benchmarks/sbayesrc_{gwas}.tsv"
  conda:
    "../envs/main.yaml"
  params:
    ldref=lambda w: f"{sbayesrc_ldref_base}/{_sbrc_pop(w)}",
    annot=sbayesrc_annot,
    seed=sbayesrc_seed,
    chain_length=sbayesrc_chain_length,
    burn_in=sbayesrc_burn_in,
    config_file=config['config_file']
  log:
    "{outdir}/logs/sbayesrc-{gwas}.log"
  shell:
    "Rscript --vanilla {workflow.basedir}/scripts/sbayesrc.R --pipeline_dir {workflow.basedir} \
      --gwas {wildcards.gwas} \
      --config_file {params.config_file} \
      --gctb {input.gctb}/gctb \
      --ldref {params.ldref} \
      --annot {params.annot} \
      --seed {params.seed} \
      --chain_length {params.chain_length} \
      --burn_in {params.burn_in} \
      --threads {threads} > {log} 2>&1"
