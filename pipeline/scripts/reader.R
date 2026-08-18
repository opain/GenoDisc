# GenoDisc results reader
#
# Public API:
#   gd_open(path)             - open a bundle (tar.gz), a directory, or a legacy .rds
#   gd_gwas(gd)               - character vector of GWAS names
#   gd_blocks(gd, gwas)       - character vector of blocks present for a GWAS
#   gd_read(gd, gwas, block)  - read one block; cached on the gd object
#   gd_config(gd)             - list(flags_raw, gwas_list, pipeline_version, ...)
#   gd_manifest(gd)           - parsed manifest
#   gd_collect(gd)            - materialise the full nested list (old shape)
#
# Bundle layout (schema_version 1):
#   manifest.json, VERSION, reader.R
#   configuration/{config.rds, gwas_list.rds}
#   gwas/<GWAS>/<block>.rds        block IDs mirror file paths
#
# Legacy compat: passing a .rds file loads the monolithic object into
# memory and exposes the same API.

# The canonical set of block IDs, in the order gd_collect uses to rebuild
# the nested list. Each ID is a `/`-separated path whose components become
# nested list keys.
.gd_block_ids <- c(
  "gwas_qc",
  "snp_assoc",
  "mol_assoc/exp/fusion",
  "mol_assoc/exp/smr",
  "mol_assoc/protein/fusion",
  "mol_assoc/protein/smr",
  "mol_assoc/magma",
  "mol_assoc/nearest",
  "mol_assoc/finemap",
  "tx/drug",
  "tx/atc",
  "tx/cmap",
  "tissue"
)

.gd_new <- function(root_dir, manifest, legacy_data = NULL) {
  structure(
    list(
      root_dir    = root_dir,
      manifest    = manifest,
      legacy_data = legacy_data,
      cache       = new.env(parent = emptyenv())
    ),
    class = "gd_result"
  )
}

.gd_row_count <- function(x) {
  if (is.data.frame(x)) return(nrow(x))
  if (is.list(x) && !is.null(x$res) && is.data.frame(x$res)) return(nrow(x$res))
  NA_integer_
}

.gd_traverse <- function(x, keys) {
  for (k in keys) {
    if (is.null(x) || !is.list(x) || !(k %in% names(x))) return(NULL)
    x <- x[[k]]
  }
  x
}

.gd_assign_nested <- function(lst, keys, value) {
  if (length(keys) == 1L) {
    lst[[keys]] <- value
    return(lst)
  }
  head_key <- keys[1L]
  if (is.null(lst[[head_key]])) lst[[head_key]] <- list()
  lst[[head_key]] <- .gd_assign_nested(lst[[head_key]], keys[-1L], value)
  lst
}

.gd_synth_manifest_from_legacy <- function(rds) {
  gwas_names <- setdiff(names(rds), "configuration")
  gwas_list  <- rds$configuration$gwas_list
  gwas_meta <- lapply(gwas_names, function(g) {
    row <- if (!is.null(gwas_list)) gwas_list[gwas_list$name == g, , drop = FALSE] else NULL
    list(
      name  = g,
      label = if (!is.null(row) && nrow(row) > 0 && !is.null(row$label)) as.character(row$label) else g
    )
  })

  blocks <- setNames(vector("list", length(gwas_names)), gwas_names)
  for (g in gwas_names) {
    gdata <- rds[[g]]
    per <- list()
    for (bid in .gd_block_ids) {
      value <- .gd_traverse(gdata, strsplit(bid, "/", fixed = TRUE)[[1L]])
      if (is.null(value)) next
      per[[bid]] <- list(bytes = NA_integer_, rows = .gd_row_count(value))
    }
    blocks[[g]] <- per
  }

  list(
    schema_version   = 0L,   # legacy marker
    pipeline_version = if (!is.null(rds$configuration$repo$commit)) rds$configuration$repo$commit else NA_character_,
    created_at       = NA_character_,
    config_file      = NA_character_,
    gwas             = gwas_meta,
    blocks           = blocks
  )
}

gd_open <- function(path) {
  if (!file.exists(path)) stop("gd_open: path not found: ", path)

  # Legacy monolithic RDS
  if (!dir.exists(path) && grepl("\\.rds$", path, ignore.case = TRUE)) {
    rds <- readRDS(path)
    if (!is.list(rds) || !("configuration" %in% names(rds))) {
      stop("gd_open: file is not a GenoDisc results_package.rds")
    }
    manifest <- .gd_synth_manifest_from_legacy(rds)
    return(.gd_new(root_dir = NA_character_, manifest = manifest, legacy_data = rds))
  }

  # Tarball -> extract to tempdir
  if (!dir.exists(path) && grepl("\\.tar(\\.gz)?$", path, ignore.case = TRUE)) {
    extract_dir <- tempfile("gd_bundle_")
    dir.create(extract_dir)
    ec <- utils::untar(path, exdir = extract_dir)
    if (!identical(ec, 0L) && !is.null(ec)) stop("gd_open: failed to extract ", path)
    mf <- list.files(extract_dir, pattern = "^manifest\\.json$", recursive = TRUE, full.names = TRUE)
    if (length(mf) != 1L) stop("gd_open: expected exactly one manifest.json inside the bundle; found ", length(mf))
    root_dir <- dirname(mf[1L])
    manifest <- jsonlite::fromJSON(mf[1L], simplifyVector = FALSE)
    return(.gd_new(root_dir = root_dir, manifest = manifest))
  }

  # Directory containing manifest.json
  if (dir.exists(path)) {
    mf <- file.path(path, "manifest.json")
    if (!file.exists(mf)) stop("gd_open: no manifest.json in directory: ", path)
    manifest <- jsonlite::fromJSON(mf, simplifyVector = FALSE)
    return(.gd_new(root_dir = path, manifest = manifest))
  }

  stop("gd_open: unrecognised path (not .rds, .tar.gz, or directory): ", path)
}

gd_manifest <- function(gd) {
  stopifnot(inherits(gd, "gd_result"))
  gd$manifest
}

gd_gwas <- function(gd) {
  stopifnot(inherits(gd, "gd_result"))
  vapply(gd$manifest$gwas, function(g) g$name, character(1L))
}

gd_blocks <- function(gd, gwas) {
  stopifnot(inherits(gd, "gd_result"))
  entries <- gd$manifest$blocks[[gwas]]
  if (is.null(entries)) return(character(0L))
  names(entries)
}

.gd_is_legacy <- function(gd) !is.null(gd$legacy_data)

gd_read <- function(gd, gwas, block) {
  stopifnot(inherits(gd, "gd_result"))
  cache_key <- paste0(gwas, "|", block)
  hit <- gd$cache[[cache_key]]
  if (!is.null(hit)) return(hit$value)   # cached NULL still returns NULL

  if (.gd_is_legacy(gd)) {
    keys  <- strsplit(block, "/", fixed = TRUE)[[1L]]
    value <- .gd_traverse(gd$legacy_data[[gwas]], keys)
  } else {
    if (is.null(gd$manifest$blocks[[gwas]]) || is.null(gd$manifest$blocks[[gwas]][[block]])) {
      value <- NULL
    } else {
      file_path <- file.path(gd$root_dir, "gwas", gwas, paste0(block, ".rds"))
      value <- if (file.exists(file_path)) readRDS(file_path) else NULL
    }
  }

  gd$cache[[cache_key]] <- list(value = value)
  value
}

# Read a GWAS QC statistic (lambda_gc, max_chi2, n_sig_snp) from a gwas_qc block.
# These moved from the FOCUS munge log to the sumstat_cleaner log when the FOCUS
# munge step was dropped, so bundles packaged before that change carry them under
# focus_dat instead of cleaner_dat.
gd_qc_stat <- function(gwas_qc, field) {
  val <- gwas_qc$cleaner_dat$val[[field]]
  if (length(val) != 1 || is.na(val)) val <- gwas_qc$focus_dat$val[[field]]
  if (length(val) != 1) return(NA_real_)
  val
}

gd_config <- function(gd) {
  stopifnot(inherits(gd, "gd_result"))
  hit <- gd$cache[["__config__"]]
  if (!is.null(hit)) return(hit$value)

  if (.gd_is_legacy(gd)) {
    conf <- gd$legacy_data$configuration
    value <- list(
      flags_raw        = conf$config,
      gwas_list        = conf$gwas_list,
      pipeline_version = if (!is.null(conf$repo$commit)) conf$repo$commit else NA_character_,
      remote           = if (!is.null(conf$repo$remote)) conf$repo$remote else NA_character_,
      branch           = if (!is.null(conf$repo$branch)) conf$repo$branch else NA_character_,
      commit           = if (!is.null(conf$repo$commit)) conf$repo$commit else NA_character_
    )
  } else {
    config_rds    <- file.path(gd$root_dir, "configuration", "config.rds")
    gwas_list_rds <- file.path(gd$root_dir, "configuration", "gwas_list.rds")
    value <- list(
      flags_raw        = if (file.exists(config_rds))    readRDS(config_rds)    else character(0L),
      gwas_list        = if (file.exists(gwas_list_rds)) readRDS(gwas_list_rds) else NULL,
      pipeline_version = gd$manifest$pipeline_version,
      remote           = NA_character_,
      branch           = NA_character_,
      commit           = NA_character_
    )
  }

  gd$cache[["__config__"]] <- list(value = value)
  value
}

gd_collect <- function(gd) {
  stopifnot(inherits(gd, "gd_result"))
  out <- list()
  for (g in gd_gwas(gd)) {
    per <- list()
    for (bid in .gd_block_ids) {
      value <- gd_read(gd, g, bid)
      keys <- strsplit(bid, "/", fixed = TRUE)[[1L]]
      per <- .gd_assign_nested(per, keys, value)
    }
    out[[g]] <- per
  }

  conf  <- gd_config(gd)
  out$configuration <- list(
    repo = list(
      remote = conf$remote,
      branch = conf$branch,
      commit = if (!is.na(conf$commit)) conf$commit else conf$pipeline_version
    ),
    config    = conf$flags_raw,
    gwas_list = conf$gwas_list
  )
  out
}

print.gd_result <- function(x, ...) {
  n <- length(gd_gwas(x))
  cat(sprintf("<gd_result> %s bundle, %d GWAS, pipeline v%s\n",
              if (.gd_is_legacy(x)) "legacy" else "split",
              n,
              if (!is.na(x$manifest$pipeline_version)) x$manifest$pipeline_version else "?"))
  invisible(x)
}
