##################
# PoPS — Polygenic Priority Score (flag: pops; requires magma_gene)
##################
# Genome-wide gene prioritisation from MAGMA gene Z-scores + a published feature
# matrix. PoPS features are Ensembl-gene-ID keyed, so a DEDICATED MAGMA run against
# the ENSGID annotation (magma_ensgid_annot) is used - the existing Entrez
# magma_gene output is untouched. Three steps: per-GWAS ENSGID-MAGMA -> one-off
# feature munge (cached under resdir) -> per-GWAS pops.py.
#
# NOTE (flagged for the empirical benchmark in the plan): the one-off munge gunzips
# the ~7 GB feature matrix (large uncompressed footprint) and pops.py's
# --num_feature_chunks must equal the number of munged column chunks; both are
# handled dynamically below but should be confirmed on the first full-matrix run.

_pops_features = config.get("pops_features", "NA")
_pops_features = _pops_features if _pops_features not in (None, "NA", "") \
    else f"{resdir}/data/pops/features_raw/PoPS.features.txt.gz"
_pops_gene_annot = config.get("pops_gene_annot", "NA")
_pops_gene_annot = _pops_gene_annot if _pops_gene_annot not in (None, "NA", "") \
    else f"{resdir}/software/pops/example/data/utils/gene_annot_jun10.txt"
_pops_control = f"{resdir}/software/pops/example/data/utils/features_jul17_control.txt"
_pops_py = f"{resdir}/software/pops/pops.py"
_pops_munge_py = f"{resdir}/software/pops/munge_feature_directory.py"
_pops_munged_dir = f"{resdir}/data/pops/munged"

# ENSGID MAGMA gene analysis for PoPS (snp-wise=mean, as PoPS expects).
rule pops_magma:
  input:
    "{outdir}/results/{gwas}/gwas_sumstat/{gwas}.cleaned.gz",
    rules.magma_ensgid_annot.output.annot
  output:
    raw="{outdir}/results/{gwas}/pops/{gwas}.pops_magma.genes.raw",
    out="{outdir}/results/{gwas}/pops/{gwas}.pops_magma.genes.out"
  benchmark:
    "{outdir}/benchmarks/pops_magma_{gwas}.tsv"
  conda:
    "../envs/main.yaml"
  params:
    resdir=resdir,
    population= lambda w: gwas_list_df.loc[gwas_list_df['name'] == "{}".format(w.gwas), 'population'].iloc[0]
  log:
    "{outdir}/logs/pops_magma-{gwas}.log"
  shell:
    "(set -eo pipefail; mkdir -p {outdir}/results/{wildcards.gwas}/pops; "
    "gzip -f -d -c {outdir}/results/{wildcards.gwas}/gwas_sumstat/{wildcards.gwas}.cleaned.gz > {outdir}/results/{wildcards.gwas}/gwas_sumstat/{wildcards.gwas}.cleaned; "
    "{params.resdir}/software/magma/magma "
    "  --bfile {params.resdir}/data/magma_ref/1kg_hgdp_{params.population} "
    "  --pval {outdir}/results/{wildcards.gwas}/gwas_sumstat/{wildcards.gwas}.cleaned use=SNP,P ncol=N "
    "  --gene-model snp-wise=mean "
    "  --gene-annot {params.resdir}/data/magma/pops_ensgid.genes.annot "
    "  --out {outdir}/results/{wildcards.gwas}/pops/{wildcards.gwas}.pops_magma; "
    "rm -f {outdir}/results/{wildcards.gwas}/gwas_sumstat/{wildcards.gwas}.cleaned) > {log} 2>&1"

# One-off feature munge (cached under resdir; independent of any GWAS).
rule pops_munge:
  input:
    features=_pops_features,
    gene_annot=_pops_gene_annot,
    munge_py=_pops_munge_py
  output:
    rows=f"{_pops_munged_dir}/pops_features.rows.txt"
  benchmark:
    f"{resdir}/benchmarks/pops_munge.tsv"
  conda:
    "../envs/pops.yaml"
  params:
    munged=_pops_munged_dir, featdir=f"{resdir}/data/pops/features_raw",
    munge_py=_pops_munge_py, gene_annot=_pops_gene_annot
  resources:
    # munge reads the ~7 GB gzipped feature matrix straight into pandas (which
    # auto-decompresses the .gz, so no 21.5 GB temp file / page-cache blow-up);
    # measured peak ~39 GB RSS, so 64 GB is enough.
    mem_mb=64000
  log:
    f"{resdir}/logs/pops_munge.log"
  shell:
    "(set -eo pipefail; mkdir -p {params.munged}; "
    "python {params.munge_py} "
    "  --gene_annot_path {params.gene_annot} "
    "  --feature_dir {params.featdir} "
    "  --save_prefix {params.munged}/pops_features "
    "  --max_cols 500) > {log} 2>&1"

# Per-GWAS PoPS prediction. num_feature_chunks is derived from the munged outputs.
rule pops_run:
  input:
    rows=f"{_pops_munged_dir}/pops_features.rows.txt",
    magma_raw="{outdir}/results/{gwas}/pops/{gwas}.pops_magma.genes.raw",
    magma_out="{outdir}/results/{gwas}/pops/{gwas}.pops_magma.genes.out",
    control=_pops_control,
    pops_py=_pops_py,
    gene_annot=_pops_gene_annot
  output:
    preds="{outdir}/results/{gwas}/pops/{gwas}.GW.pops.preds"
  benchmark:
    "{outdir}/benchmarks/pops_run_{gwas}.tsv"
  conda:
    "../envs/pops.yaml"
  resources:
    mem_mb=16000
  params:
    munged=_pops_munged_dir, pops_py=_pops_py, control=_pops_control, gene_annot=_pops_gene_annot
  log:
    "{outdir}/logs/pops_run-{gwas}.log"
  shell:
    "(set -eo pipefail; "
    "N=$(ls {params.munged}/pops_features.cols.*.txt 2>/dev/null | wc -l); "
    "if [ \"$N\" -eq 0 ]; then echo 'no munged feature chunks found' >&2; exit 1; fi; "
    "python {params.pops_py} "
    "  --gene_annot_path {params.gene_annot} "
    "  --feature_mat_prefix {params.munged}/pops_features "
    "  --num_feature_chunks $N "
    "  --control_features_path {params.control} "
    "  --magma_prefix {outdir}/results/{wildcards.gwas}/pops/{wildcards.gwas}.pops_magma "
    "  --out_prefix {outdir}/results/{wildcards.gwas}/pops/{wildcards.gwas}.GW.pops) > {log} 2>&1"
