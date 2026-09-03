#!/usr/bin/Rscript
# THIS IS A TRACKED COPY. The web app actually runs the live copy at
# GenoDisc/bin/inspect_gwas_header.R, not this file directly. They're separate
# copies, NOT a symlink (symlinks across this CIFS mount don't resolve
# correctly when accessed from the HPC side). Edit one, then manually copy the
# change to the other. Whichever you edit first, do the copy in the same
# sitting - don't let them drift. Same split as run_genodisc.sh - see that
# file's header for the full explanation.
#
# Reads the header + first 1000 rows of a GWAS sumstats file and reports how
# each column was interpreted, reusing the exact logic the real pipeline runs
# on every job - scripts/sumstat_cleaner.R calls the same head_interp() (via
# scripts/functions/sumstat_cleaner_functions.R) during actual analysis. This
# is a read-only preview: it never writes/modifies the uploaded file.
suppressMessages(library("optparse"))
suppressMessages(library("data.table"))
suppressMessages(library("jsonlite"))

option_list = list(
  make_option("--sumstats", action="store", default=NA, type='character',
              help="Path to summary statistics file [required]"),
  make_option("--pipeline_dir", action="store", default=NA, type="character",
              help="Path to the pipeline directory [required]")
)
opt = parse_args(OptionParser(option_list=option_list))

source(file.path(opt$pipeline_dir, 'scripts', 'functions', 'sumstat_cleaner_functions.R'))

sub_ss <- fread(opt$sumstats, nrows = 1000)
header_interp <- head_interp(sub_ss)

# Same detection the real pipeline uses as its own fallback (scripts/sumstat_cleaner.R):
# N is usable if present directly, or derivable from N_CAS + N_CON.
n_detected <- ('N' %in% header_interp$Interpreted) ||
  all(c('N_CAS', 'N_CON') %in% header_interp$Interpreted)

# Same hard-requirement checks as format_header() in sumstat_cleaner_functions.R,
# called directly here (not via format_header) so this stays a read-only preview.
required_present <- list(
  a1_a2 = all(c('A1', 'A2') %in% header_interp$Interpreted),
  effect_size = any(c('OR', 'BETA', 'Z') %in% header_interp$Interpreted),
  variant_id = ('SNP' %in% header_interp$Interpreted) || all(c('CHR', 'BP') %in% header_interp$Interpreted)
)

cat(toJSON(list(
  columns = header_interp,
  n_detected = n_detected,
  required_present = required_present
), auto_unbox = TRUE, na = "null"))
