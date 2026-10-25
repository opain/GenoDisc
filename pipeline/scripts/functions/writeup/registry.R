#!/usr/bin/Rscript

# Registry + citation handling for the write-up generator.
#
# Loads writeup/registry.yaml and writeup/references.bib, resolves tool versions
# (live from envs/*.yaml where possible), and formats inline citations while
# recording which BibTeX keys were actually used so the reference list can be
# subset to the run. Any tool/dataset used without a registry entry, or whose
# bibkey is absent from references.bib, is FLAGGED (never silently dropped).

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || (length(a) == 1 && is.na(a))) b else a

# Version of the generator itself (bumped when the templates/logic change). Recorded
# into the manifest and the document headers.
WRITEUP_GENERATOR_VERSION <- "0.1.0"

# ---- resolved-config / params accessors (read from the run) -----------------

# Raw value of a resolved-config key (parsed merged_config list). NA if absent.
cfg_val <- function(ctx, key, default = NA) {
  v <- ctx$config[[key]]
  if (is.null(v)) return(default)
  v
}

# TRUE if a T/F config flag is on.
cfg_on <- function(ctx, key) {
  v <- cfg_val(ctx, key, "F")
  identical(toupper(as.character(v)[1]), "T")
}

# A path/string config param is "set" if not NA / NULL / "" / "NA".
cfg_set <- function(ctx, key) {
  v <- cfg_val(ctx, key, NA)
  !(is.null(v) || length(v) == 0 || (length(v) == 1 && (is.na(v) || v %in% c("", "NA"))))
}

# A list-valued config param (e.g. gtex_weights) as a clean character vector.
cfg_list <- function(ctx, key) {
  v <- cfg_val(ctx, key, NULL)
  if (is.null(v)) return(character(0))
  v <- as.character(unlist(v))
  v[nzchar(v) & v != "NA"]
}

# A fixed pipeline-version parameter from writeup/params.yaml, e.g. param(ctx,"clump","r2").
param <- function(ctx, ...) {
  x <- ctx$params
  for (k in c(...)) {
    if (is.null(x)) return(NULL)
    x <- x[[k]]
  }
  x
}

# ---- loading ----------------------------------------------------------------

writeup_load_registry <- function(pipeline_dir) {
  path <- file.path(pipeline_dir, "writeup", "registry.yaml")
  reg <- yaml::read_yaml(path)
  # Index software + datasets by key for O(1) lookup; keep `kind` for messages.
  idx <- list()
  for (e in reg$software) { e$kind <- "software"; idx[[e$key]] <- e }
  for (e in reg$datasets) { e$kind <- "dataset";  idx[[e$key]] <- e }
  list(raw = reg, idx = idx)
}

# Parse the @type{key, ...} openers of references.bib into a character vector of keys.
writeup_load_bibkeys <- function(pipeline_dir) {
  path <- file.path(pipeline_dir, "writeup", "references.bib")
  lines <- readLines(path, warn = FALSE)
  m <- regmatches(lines, regexpr("^@[A-Za-z]+\\{[^,]+,", lines))
  m <- m[nzchar(m)]
  keys <- sub("^@[A-Za-z]+\\{", "", m)
  keys <- sub(",$", "", keys)
  trimws(keys)
}

# Read a pinned conda/pip version out of an env YAML. Matches "  - plink=1.90b6.21",
# "  - r-susier=0.12.16" and pip "    - scipy==0.18". Returns NA if not found.
# Uses fixed-string prefix matching (no regex built from the package name).
writeup_env_version <- function(pipeline_dir, env_file, pkg) {
  path <- file.path(pipeline_dir, "envs", env_file)
  if (!file.exists(path)) return(NA_character_)
  lines <- trimws(readLines(path, warn = FALSE))
  # Check "==" (pip) before "=" (conda) since "=" is a prefix of "==".
  for (sep in c("==", "=")) {
    prefix <- paste0("- ", pkg, sep)
    hit <- lines[startsWith(lines, prefix)]
    if (length(hit) > 0) {
      v <- substring(hit[1], nchar(prefix) + 1L)
      v <- sub("#.*$", "", v)          # strip trailing comment
      v <- sub("\\s.*$", "", v)        # strip trailing whitespace/build
      v <- trimws(v)
      if (nzchar(v)) return(v)
    }
  }
  NA_character_
}

# ---- context ----------------------------------------------------------------

# A mutable context (environment) threaded through all block builders. Holds the
# registry/bib, formatting options, and accumulators for cited keys + flags/warnings.
writeup_new_ctx <- function(pipeline_dir, config, opt = writeup_fmt_opts(),
                            schema_version = NA_integer_) {
  ctx <- new.env(parent = emptyenv())
  ctx$pipeline_dir   <- pipeline_dir
  ctx$config         <- config            # parsed merged_config (named list)
  ctx$opt            <- opt
  ctx$registry       <- writeup_load_registry(pipeline_dir)
  ctx$params         <- yaml::read_yaml(file.path(pipeline_dir, "writeup", "params.yaml"))
  ctx$bib_keys       <- writeup_load_bibkeys(pipeline_dir)
  ctx$bibkeys_used   <- character(0)
  ctx$warnings       <- character(0)      # reconciliation / citation flags
  ctx$ran_blocks     <- character(0)      # block ids for which a paragraph was emitted
  ctx$schema_version <- schema_version
  ctx$generator_version <- WRITEUP_GENERATOR_VERSION
  ctx
}

writeup_warn <- function(ctx, msg) {
  if (!msg %in% ctx$warnings) ctx$warnings <- c(ctx$warnings, msg)
  invisible(NULL)
}

# ---- citation API -----------------------------------------------------------

# Resolve a registry entry's display version string.
writeup_resolve_version <- function(entry, pipeline_dir) {
  if (isTRUE(entry$unpinned)) return("version not pinned upstream")
  if (!is.null(entry$conda_env) && !is.null(entry$conda_pkg)) {
    v <- writeup_env_version(pipeline_dir, entry$conda_env, entry$conda_pkg)
    if (!is.na(v)) return(paste0("v", v))
    return(entry$version %||% "version unknown")
  }
  entry$version %||% "version unknown"
}

# Look up a registry key, record its bibkeys as used, flag problems, and return a
# small list(name, version, cite). `cite` is a pandoc-style "[@key1; @key2]".
writeup_ref <- function(ctx, key) {
  entry <- ctx$registry$idx[[key]]
  if (is.null(entry)) {
    writeup_warn(ctx, sprintf("No registry entry for '%s' — citation and version missing from the text.", key))
    return(list(name = key, version = "version unknown", cite = sprintf("[CITATION MISSING: %s]", key)))
  }
  bibkeys <- unlist(entry$bibkey)
  if (length(bibkeys) == 0) {
    writeup_warn(ctx, sprintf("Registry entry '%s' has no bibkey.", key))
  }
  missing <- setdiff(bibkeys, ctx$bib_keys)
  if (length(missing) > 0) {
    writeup_warn(ctx, sprintf("Registry entry '%s' references bibkey(s) absent from references.bib: %s.",
                              key, paste(missing, collapse = ", ")))
  }
  ctx$bibkeys_used <- union(ctx$bibkeys_used, bibkeys)
  ver <- writeup_resolve_version(entry, ctx$pipeline_dir)
  if (identical(ver, "version not pinned upstream")) {
    writeup_warn(ctx, sprintf("%s version is not recorded anywhere in the run (not pinned upstream); stated as unknown.", entry$name))
  }
  notes <- entry$notes %||% ""
  if (grepl("VERIFY", notes, fixed = TRUE)) {
    writeup_warn(ctx, sprintf("Citation for '%s' is not from the curated in-app table — verify before publication.", entry$name))
  }
  cite <- paste0("[", paste0("@", bibkeys, collapse = "; "), "]")
  list(name = entry$name, version = ver, cite = cite)
}

# "Name (version) [@key]" — the usual inline tool phrase for methods text.
writeup_tool <- function(ctx, key) {
  r <- writeup_ref(ctx, key)
  paste0(r$name, " (", r$version, ") ", r$cite)
}

# ---- coverage checks (used by tests/testthat/test-writeup-references.R) ------

# Every bibkey named by any registry entry must exist in references.bib.
writeup_missing_bibkeys <- function(pipeline_dir) {
  reg  <- writeup_load_registry(pipeline_dir)
  keys <- writeup_load_bibkeys(pipeline_dir)
  used <- unique(unlist(lapply(reg$idx, function(e) unlist(e$bibkey))))
  setdiff(used, keys)
}

# Every block-enabling config flag must be covered by >= 1 registry `used_by`.
# `flags` defaults to the canonical block flags enumerated below.
writeup_block_flags <- function() {
  c("clump", "cojo", "manc_cojo", "finemap", "sbayesrc", "sbayesrc_r", "sumher",
    "ldsc", "gencor_gwas_list", "gencor_within_gwas_list",
    "magma_gene", "magma_drugtargetor", "magma_drugtargetor_gls", "magma_pathway",
    "tissue_magma",
    "twas_panel_psychencode", "twas_panel_fusion", "twas_panel_fetal_brain",
    "twas_panel_herv", "external_weights",
    "pwas_panel_rosmap", "pwas_panel_banner",
    "smr_expression_panel_psychencode", "smr_expression_panel_metabrain_basalganglia",
    "smr_expression_panel_metabrain_cerebellum", "smr_expression_panel_metabrain_cortex",
    "smr_expression_panel_metabrain_hippocampus", "smr_expression_panel_metabrain_spinalcord",
    "smr_expression_panel_eqtlgen", "smr_protein_panel_rosmap",
    "twas_gsea_drugtargetor", "twas_gsea_drugtargetor_nondirectional",
    "twas_gsea_pathway", "twas_gsea_cmap", "drug_targetor_atc_gls", "gcsc")
}

writeup_uncovered_flags <- function(pipeline_dir, flags = writeup_block_flags()) {
  reg <- writeup_load_registry(pipeline_dir)
  covered <- unique(unlist(lapply(reg$idx, function(e) unlist(e$used_by))))
  setdiff(flags, covered)
}
