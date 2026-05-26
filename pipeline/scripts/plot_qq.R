#!/usr/bin/Rscript
suppressMessages(library("optparse"))

option_list = list(
  make_option("--gwas", action="store", default=NA, type='character',
              help="Name of GWAS [required]"),
  make_option("--config_file", action="store", default=NA, type='character',
              help="Path to config file [required]")
)

option_list <- c(option_list, list(
  make_option("--pipeline_dir", action="store", default=NA, type="character",
              help="Path to the pipeline directory [required]")
))

opt = parse_args(OptionParser(option_list=option_list))
options(pipeline_dir = opt$pipeline_dir)

suppressMessages({
  library(data.table)
  library(ggplot2)
})
source(file.path(opt$pipeline_dir, 'scripts', 'functions', 'utils_functions.R'))

outdir <- read_param(config = opt$config_file, param = 'outdir', return_obj = F)

ss <- fread(paste0(outdir, '/results/', opt$gwas, '/gwas_sumstat/', opt$gwas, '.cleaned.gz'),
            select = "P")
ss <- ss[!is.na(P) & P > 0 & P <= 1]
setorder(ss, P)

n <- nrow(ss)
obs <- -log10(ss$P)
exp_q <- -log10(ppoints(n))

chi2 <- qchisq(1 - ss$P, df = 1)
lambda_gc <- median(chi2) / qchisq(0.5, 1)

p <- ggplot(data.frame(exp_q = exp_q, obs = obs), aes(x = exp_q, y = obs)) +
  geom_point(size = 0.6) +
  geom_abline(slope = 1, intercept = 0, colour = "red", linetype = "dashed") +
  labs(x = expression(Expected~-log[10](italic(P))),
       y = expression(Observed~-log[10](italic(P))),
       subtitle = bquote(lambda[GC] == .(round(lambda_gc, 3)))) +
  theme_bw(base_size = 10)

out_path <- paste0(outdir, '/results/', opt$gwas, '/gwas_sumstat/', opt$gwas, '.qq_plot.png')
ggsave(out_path, plot = p, width = 6, height = 6, dpi = 200)
