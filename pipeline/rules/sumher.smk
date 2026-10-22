##################
# LDAK SumHer SNP-heritability (BLD-LDAK, pre-computed HapMap3 tagging)
##################
#
# A THIRD, independent SNP-h2 estimator (flag: sumher), alongside LDSC and the
# two SBayesRC routes. Runs LDAK's SumHer (--sum-hers) with the BLD-LDAK model
# over a pre-computed HapMap3 tagging file. Genome-wide; the tagging file is
# picked by the per-GWAS population (EUR->gbr, EAS->eas, SAS->sas, AFR->afr;
# AMR is skipped). Observed scale only (no prevalence/ascertainment). A different
# heritability model from LDSC/SBayesRC, so estimates are not expected to agree.
#
# References: the LDAK binary (download_ldak) + BLD-LDAK HapMap3 tagging files
# (download_sumher_tagging). config.yaml sumher_ldak / sumher_ref default to NA =
# use the copies staged in resdir; set an explicit path to reuse an existing copy.
_sumher_ldak_cfg = config.get("sumher_ldak", "NA")
sumher_ldak = (f"{resdir}/software/ldak/ldak7.01.linux"
               if _sumher_ldak_cfg in (None, "NA") else _sumher_ldak_cfg)

_sumher_ref_cfg = config.get("sumher_ref", "NA")
sumher_ref_base = (f"{resdir}/data/sumher_tagging/bld_ldak/hapmap"
                   if _sumher_ref_cfg in (None, "NA") else _sumher_ref_cfg)

sumher_cutoff      = config.get("sumher_cutoff", 0.01)
sumher_max_missing = config.get("sumher_max_missing", 0.2)
sumher_threads     = int(config.get("sumher_threads", 1))

# GWAS population (gnomAD super-pop code) -> LDAK BLD-LDAK tagging population code.
# CSA (central/south Asian) reuses the SAS tagging (closest available); AMR/MID
# have no tagging file and are excluded upstream (gwas_list_df_sumher).
_sumher_pop_map = {'EUR': 'gbr', 'EAS': 'eas', 'AFR': 'afr', 'CSA': 'sas'}

def _sumher_pop(w):
    return gwas_list_df.loc[gwas_list_df['name'] == "{}".format(w.gwas), 'population'].iloc[0]

def _sumher_tagpop(w):
    return _sumher_pop_map[_sumher_pop(w)]

# SumHer from a pre-computed tagfile is a light weighted regression; 8 GB is ample.
def get_mem_mb_sumher(wildcards, attempt):
    return attempt * 8000

rule sumher:
  # cpus (read by the SLURM profile) must equal threads; SumHer barely scales, so
  # sumher_threads defaults to 1.
  resources:
    mem_mb=get_mem_mb_sumher,
    cpus=sumher_threads
  threads: sumher_threads
  input:
    sumstats="{outdir}/results/{gwas}/gwas_sumstat/{gwas}.cleaned.gz",
    ldak=sumher_ldak,
    tagfile=lambda w: f"{sumher_ref_base}/{_sumher_tagpop(w)}/bld.ldak.hapmap.{_sumher_tagpop(w)}.tagging"
  output:
    "{outdir}/results/{gwas}/sumher/{gwas}.sumher.clean.csv"
  benchmark:
    "{outdir}/benchmarks/sumher_{gwas}.tsv"
  conda:
    "../envs/sumher.yaml"
  params:
    ldak=sumher_ldak,
    cutoff=sumher_cutoff,
    max_missing=sumher_max_missing,
    config_file=config['config_file']
  log:
    "{outdir}/logs/sumher-{gwas}.log"
  shell:
    "Rscript --vanilla {workflow.basedir}/scripts/sumher.R --pipeline_dir {workflow.basedir} \
      --gwas {wildcards.gwas} \
      --config_file {params.config_file} \
      --ldak {params.ldak} \
      --tagfile {input.tagfile} \
      --cutoff {params.cutoff} \
      --max_missing {params.max_missing} \
      --threads {resources.cpus} > {log} 2>&1"
