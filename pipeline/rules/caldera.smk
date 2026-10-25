##################
# CALDERA — effector-gene probabilities per credible set
# (flag: caldera; requires pops AND finemap, enforced at config parse)
##################
# Applies the published PRETRAINED CALDERA model (shipped in the repo) to each
# SuSiE credible set, using the PoPS gene scores (Phase 2) and the per-variant
# credible-set table (Phase 0.5). Also emits a nearest-gene baseline per set, and
# joins the Phase 1 coding summary when available. See scripts/run_caldera.R.

_caldera_path = config.get("caldera_repo", "NA")
_caldera_path = _caldera_path if _caldera_path not in (None, "NA", "") else f"{resdir}/software/caldera"

rule caldera_run:
  input:
    cv="{outdir}/results/{gwas}/finemap/{gwas}.GW.credset_variants.csv",
    preds="{outdir}/results/{gwas}/pops/{gwas}.GW.pops.preds",
    model=rules.download_caldera.output.model
  output:
    "{outdir}/results/{gwas}/caldera/{gwas}.GW.caldera.csv"
  benchmark:
    "{outdir}/benchmarks/caldera_run_{gwas}.tsv"
  conda:
    "../envs/main.yaml"
  params:
    config_file=config['config_file'],
    caldera_path=_caldera_path
  log:
    "{outdir}/logs/caldera_run-{gwas}.log"
  shell:
    "Rscript --vanilla {workflow.basedir}/scripts/run_caldera.R --pipeline_dir {workflow.basedir} \
      --gwas {wildcards.gwas} \
      --config_file {params.config_file} \
      --caldera_path {params.caldera_path} \
      --pops_preds {outdir}/results/{wildcards.gwas}/pops/{wildcards.gwas}.GW.pops.preds > {log} 2>&1"
