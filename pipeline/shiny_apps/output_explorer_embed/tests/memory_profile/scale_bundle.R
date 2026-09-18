#!/usr/bin/env Rscript
# Duplicate the single GWAS in a GenoDisc results bundle to synthesize an N-GWAS bundle.
# Usage: Rscript scale_bundle.R --src <in.tar.gz> --n <N> --out <out.tar.gz>

suppressPackageStartupMessages({
  library(jsonlite)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) stop("missing arg: ", flag)
  args[i + 1L]
}
SRC <- get_arg("--src")
N   <- as.integer(get_arg("--n"))
OUT <- get_arg("--out")
stopifnot(N >= 1L, file.exists(SRC))

work <- tempfile("scale_")
dir.create(work)
on.exit(unlink(work, recursive = TRUE), add = TRUE)

message(sprintf("extract %s -> %s", SRC, work))
ec <- utils::untar(SRC, exdir = work)
if (!identical(ec, 0L) && !is.null(ec)) stop("untar failed")

root <- list.files(work, pattern = "^manifest\\.json$", recursive = TRUE, full.names = TRUE)
if (length(root) != 1L) stop("expected exactly one manifest.json")
pkg_dir <- dirname(root)

mf <- jsonlite::fromJSON(root, simplifyVector = FALSE)
orig_names <- vapply(mf$gwas, function(g) g$name, character(1L))
if (length(orig_names) != 1L) {
  stop("scale_bundle only handles single-GWAS source bundles; got ", length(orig_names))
}
orig_name <- orig_names[1L]
orig_gwas_dir  <- file.path(pkg_dir, "gwas", orig_name)
orig_gwas_meta <- mf$gwas[[1L]]
orig_blocks    <- mf$blocks[[orig_name]]

new_names <- if (N == 1L) orig_name else {
  c(orig_name, sprintf("%s_v%02d", orig_name, 2:N))
}

if (N > 1L) {
  for (nm in new_names[-1L]) {
    dst <- file.path(pkg_dir, "gwas", nm)
    dir.create(dst, recursive = TRUE)
    fs <- list.files(orig_gwas_dir, recursive = TRUE, all.files = TRUE, full.names = FALSE)
    for (f in fs) {
      src_f <- file.path(orig_gwas_dir, f)
      dst_f <- file.path(dst, f)
      dir.create(dirname(dst_f), showWarnings = FALSE, recursive = TRUE)
      file.copy(src_f, dst_f, copy.mode = TRUE)
    }
  }
}

new_gwas_meta <- lapply(new_names, function(nm) {
  m <- orig_gwas_meta
  m$name  <- nm
  m$label <- paste0(orig_gwas_meta$label, if (nm == orig_name) "" else sprintf(" (dup %s)", sub(".*_v", "v", nm)))
  m
})
new_blocks <- setNames(rep(list(orig_blocks), length(new_names)), new_names)

mf$gwas   <- new_gwas_meta
mf$blocks <- new_blocks
jsonlite::write_json(mf, root, auto_unbox = TRUE, pretty = TRUE, null = "null")

gwas_list_path <- file.path(pkg_dir, "configuration", "gwas_list.rds")
if (file.exists(gwas_list_path)) {
  gl <- readRDS(gwas_list_path)
  gl <- as.data.table(gl)
  base_row <- gl[name == orig_name][1L]
  rows <- lapply(new_names, function(nm) {
    r <- copy(base_row)
    r[, name := nm]
    if ("label" %in% names(r)) r[, label := paste0(base_row$label, if (nm == orig_name) "" else sprintf(" (dup %s)", sub(".*_v", "v", nm)))]
    r
  })
  new_gl <- rbindlist(rows, use.names = TRUE, fill = TRUE)
  saveRDS(new_gl, gwas_list_path)
}

out_abs <- normalizePath(dirname(OUT), mustWork = TRUE)
out_file <- basename(OUT)
old_wd <- getwd()
setwd(work)
on.exit(setwd(old_wd), add = TRUE)
tar_bin <- Sys.which("tar")
if (tar_bin == "") stop("no tar in PATH")
cmd <- sprintf("%s czf %s package", shQuote(tar_bin),
               shQuote(file.path(out_abs, out_file)))
ec2 <- system(cmd)
if (ec2 != 0L) stop("tar failed")

size_b <- file.size(file.path(out_abs, out_file))
message(sprintf("wrote %s  (%d GWAS, %.1f MB)", file.path(out_abs, out_file), N, size_b / 1024 / 1024))
