##################
# Run GCTB SBayesRC (SNP-h2 + polygenicity)
##################
#
# Genome-wide SBayesRC (one run per GWAS, not per chromosome). We keep only the
# model-level parameters (SNP-h2 + polygenicity) from the .parRes output; the
# per-SNP effects (.snpRes) are not retained. LD reference is chosen from the
# per-GWAS population column (EUR/EAS/AFR; see gwas_list_df_sbrc in sumstat_qc.smk).
#
# Reference data (eigen LD + baseline annotation) is shared with GenoPred.
# Set explicit paths in config.yaml (sbayesrc_ldref / sbayesrc_annot); NA falls
# back to GenoPred's staged copies. sbayesrc_ldref is a directory with
# per-population subdirs (<ldref>/<POP>/{ldm.info,snp.info,block*.eigen.bin}).
_GENOPRED_RES = "/scratch/prj/neurohackpain/GenoPred/pipeline_resources"

_sbrc_ldref_cfg = config.get("sbayesrc_ldref", "NA")
sbayesrc_ldref_base = (f"{_GENOPRED_RES}/data/sbayesrc_ref"
                       if _sbrc_ldref_cfg in (None, "NA") else _sbrc_ldref_cfg)

_sbrc_annot_cfg = config.get("sbayesrc_annot", "NA")
sbayesrc_annot = (f"{_GENOPRED_RES}/data/sbayesrc_annot/annot_baseline2.2.txt"
                  if _sbrc_annot_cfg in (None, "NA") else _sbrc_annot_cfg)

def get_mem_mb_sbayesrc(wildcards, attempt):
    return attempt * 32000

rule sbayesrc:
  resources:
    mem_mb=get_mem_mb_sbayesrc
  threads: 4
  input:
    sumstats="{outdir}/results/{gwas}/gwas_sumstat/{gwas}.cleaned.gz",
    gctb=rules.download_gctb.output
  output:
    "{outdir}/results/{gwas}/sbayesrc/{gwas}.GW.sbayesrc.clean.csv"
  benchmark:
    "{outdir}/benchmarks/sbayesrc_{gwas}.tsv"
  conda:
    "../envs/main.yaml"
  params:
    ldref=lambda w: f"{sbayesrc_ldref_base}/" + gwas_list_df.loc[gwas_list_df['name'] == "{}".format(w.gwas), 'population'].iloc[0],
    annot=sbayesrc_annot,
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
      --threads {threads} > {log} 2>&1"
