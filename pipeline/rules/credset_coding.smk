##################
# Credible-set coding annotation (flag: credset_coding_annot; requires finemap)
##################
# Annotate each SuSiE credible-set variant with its predicted coding consequence
# (AlphaMissense missense pathogenicity + bcftools csq for stop/frameshift/splice/
# inframe/start-lost), and summarise per credible set the PIP carried by
# protein-altering variants. GRCh37/hg19 throughout. See scripts/
# annotate_credset_coding.R and scripts/extract_credset_variants.R.

# Resource paths: use the config override if set, else the staged default.
def _resource_path(key, default):
    v = config.get(key, "NA")
    return v if v not in (None, "NA", "") else default

_am_path    = _resource_path("alphamissense_hg19", f"{resdir}/data/alphamissense/AlphaMissense_hg19.bgz")
_csq_gff3   = _resource_path("csq_gff3",  f"{resdir}/data/csq/Homo_sapiens.GRCh37.87.gff3.gz")
_csq_fasta  = _resource_path("csq_fasta", f"{resdir}/data/csq/Homo_sapiens.GRCh37.dna.primary_assembly.fa")
_coding_pip_thresh = config.get("credset_coding_pip_thresh", 0.1)

rule credset_coding_annot:
  input:
    cv="{outdir}/results/{gwas}/finemap/{gwas}.GW.credset_variants.csv",
    am=_am_path,
    amtbi=_am_path + ".tbi",
    gff3=_csq_gff3,
    fasta=_csq_fasta,
    fai=_csq_fasta + ".fai"
  output:
    summary="{outdir}/results/{gwas}/finemap/{gwas}.GW.credset_coding_summary.csv",
    long="{outdir}/results/{gwas}/finemap/{gwas}.GW.credset_coding.csv"
  benchmark:
    "{outdir}/benchmarks/credset_coding_annot_{gwas}.tsv"
  conda:
    "../envs/coding_annot.yaml"
  params:
    config_file=config['config_file'],
    am=_am_path,
    gff3=_csq_gff3,
    fasta=_csq_fasta,
    thr=_coding_pip_thresh
  log:
    "{outdir}/logs/credset_coding_annot-{gwas}.log"
  shell:
    "Rscript --vanilla {workflow.basedir}/scripts/annotate_credset_coding.R --pipeline_dir {workflow.basedir} \
      --gwas {wildcards.gwas} \
      --config_file {params.config_file} \
      --alphamissense {params.am} \
      --gff3 {params.gff3} \
      --fasta {params.fasta} \
      --pip_thresh {params.thr} > {log} 2>&1"
