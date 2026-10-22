#!/usr/bin/Rscript
library(data.table)
library(optparse)

option_list = list(
  make_option("--resdir", type="character", default="resources")
)
opt = parse_args(OptionParser(option_list=option_list))

###################
# Format the dense 1KG+HGDP reference for GenoDisc
###################
# The reference tarball is downloaded and extracted by rule download_1kg_hgdp, leaving
# {resdir}/data/1kg_hgdp/ref/ with, per autosome: a plink2 pgen fileset (ref.chr{chr}.{pgen,pvar,psam})
# and a combined cross-population .rds (ref.chr{chr}.rds; the GenoDisc sumstat-cleaner schema
# CHR,SNP,BP_GRCh36,BP_GRCh37,BP_GRCh38,A1,A2,IUPAC,REF.FRQ.{POP}), plus per-population keep files
# (keep_files/{POP}.keep) and per-population allele frequencies (freq_files/...).
#
# Here we (1) convert the pgen to the per-population PLINK1 bfiles that downstream rules expect, and
# (2) stage the combined .rds under the name the sumstat cleaner reads via --ref_chr.

output_dir<-paste0(opt$resdir, '/data/1kg_hgdp')
ref_dir<-paste0(output_dir, '/ref')

# Super-populations present in the reference (driven by the shipped keep files: EUR/EAS/AFR/AMR/CSA/MID)
pops<-gsub('\\.keep$', '', list.files(paste0(ref_dir, '/keep_files'), pattern='\\.keep$'))

# Convert each population x chromosome from pgen to a PLINK1 bfile, under the name downstream rules
# (clump/cojo/smr/twas/gcsc/feature_pred) expect: 1KG.Phase3.{POP}.MAF_001.chr{chr}.
# No MAF filter is applied: the dense panel is already the dbSNP 'common' set (MAF > 1% in at least one
# 1KG super-population), and clump/COJO re-impose --maf 0.01 at runtime. The 'MAF_001' token is kept
# only to preserve the existing downstream path templates.
for(pop in pops){
  for(chr in 1:22){
    system(paste0('plink2 --pfile ',ref_dir,'/ref.chr',chr,
                  ' --keep ',ref_dir,'/keep_files/',pop,'.keep',
                  ' --make-bed --out ',output_dir,'/1KG.Phase3.',pop,'.MAF_001.chr',chr))
  }
}

# Stage the combined cross-population .rds under the name the sumstat cleaner reads (--ref_chr).
for(chr in 1:22){
  file.rename(paste0(ref_dir,'/ref.chr',chr,'.rds'),
              paste0(output_dir,'/1KG.Phase3.MAF_001.chr',chr,'.rds'))
}

# Remove the extracted reference directory (pgen/pvar/psam, keep_files, freq_files no longer needed)
unlink(ref_dir, recursive = TRUE)
