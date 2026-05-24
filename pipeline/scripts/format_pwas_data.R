#!/usr/bin/Rscript

suppressMessages(library("optparse"))

option_list = list(
  make_option("--rosmap", action="store", default=NA, type='character',
              help="Path to ROSMAP data [required]"),
  make_option("--banner", action="store", default=NA, type='character',
              help="Path to ROSMAP data [required]"),
  make_option("--resdir", type="character", default="resources")
)

option_list <- c(option_list, list(
  make_option("--pipeline_dir", action="store", default=NA, type="character",
              help="Path to the pipeline directory [required]")
))

opt = parse_args(OptionParser(option_list=option_list))
options(pipeline_dir = opt$pipeline_dir)

library(data.table)

########
# ROSMAP
########

dir.create(paste0(opt$resdir, '/data/rosmap_twas/'))
unzip(opt$rosmap, exdir=paste0(opt$resdir, '/data/rosmap_twas/'))

pos<-fread(paste0(opt$resdir, '/data/rosmap_twas/ROSMAP.n376.fusion.WEIGHTS/train_weights.pos'))
pos$N<-376
fwrite(pos, paste0(opt$resdir, '/data/rosmap_twas/ROSMAP.n376.fusion.WEIGHTS/train_weights_withN.pos'), quote=F, sep=' ', na='NA')

########
# Banner
########

dir.create(paste0(opt$resdir, '/data/banner_twas/'))
unzip(opt$banner, exdir=paste0(opt$resdir, '/data/banner_twas/'))

pos<-fread(paste0(opt$resdir, '/data/banner_twas/Banner.n152.fusion.WEIGHTS/train_weights.pos'))
pos$N<-152
fwrite(pos, paste0(opt$resdir, '/data/banner_twas/Banner.n152.fusion.WEIGHTS/train_weights_withN.pos'), quote=F, sep=' ', na='NA')


