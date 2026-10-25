#!/usr/bin/Rscript

# Top-level write-up driver. Assembles a run-level methods section, a per-GWAS
# descriptive results section, a reconciliation footer, and the cited-subset
# reference list (BibTeX + RIS) from the collected `output`, the resolved merged
# config, and the manifest. Pure: returns a list of {relpath -> content} plus
# metadata; the caller (package_results.R) writes the files. This keeps the
# generator unit-testable and re-usable (e.g. by the Shiny app later).

library(glue)

# ---- block specification (fixed publication order) --------------------------
# Each spec: id, methods fn, results fn, and `requested` — a predicate on ctx
# that is TRUE when the run asked for this block (used for reconciliation only).

.writeup_specs <- function() {
  list(
    list(id = "preprocessing", methods = methods_preprocessing, results = NULL,
         requested = function(ctx) TRUE),
    list(id = "ldsc_h2", methods = methods_ldsc_h2, results = results_ldsc_h2,
         requested = function(ctx) cfg_on(ctx, "ldsc")),
    list(id = "sbayesrc", methods = methods_sbayesrc, results = results_sbayesrc,
         requested = function(ctx) cfg_on(ctx, "sbayesrc")),
    list(id = "sbayesrc_r", methods = methods_sbayesrc_r, results = NULL,
         requested = function(ctx) cfg_on(ctx, "sbayesrc_r")),
    list(id = "sumher", methods = methods_sumher, results = results_sumher,
         requested = function(ctx) cfg_on(ctx, "sumher")),
    list(id = "gencor", methods = methods_gencor, results = results_gencor,
         requested = function(ctx) cfg_set(ctx, "gencor_gwas_list") || cfg_on(ctx, "gencor_within_gwas_list")),
    list(id = "clump", methods = methods_clump, results = results_clump,
         requested = function(ctx) cfg_on(ctx, "clump")),
    list(id = "cojo", methods = methods_cojo, results = results_cojo,
         requested = function(ctx) cfg_on(ctx, "cojo") || cfg_on(ctx, "manc_cojo")),
    list(id = "finemap", methods = methods_finemap, results = results_finemap,
         requested = function(ctx) cfg_on(ctx, "finemap")),
    list(id = "magma_gene", methods = methods_magma_gene, results = results_magma_gene,
         requested = function(ctx) cfg_on(ctx, "magma_gene")),
    list(id = "twas", methods = methods_twas, results = results_twas,
         requested = function(ctx) any(vapply(c("twas_panel_psychencode", "twas_panel_fusion", "twas_panel_fetal_brain", "twas_panel_herv", "external_weights"), function(f) cfg_on(ctx, f), logical(1))) || length(cfg_list(ctx, "gtex_weights")) > 0),
    list(id = "pwas", methods = methods_pwas, results = results_pwas,
         requested = function(ctx) cfg_on(ctx, "pwas_panel_rosmap") || cfg_on(ctx, "pwas_panel_banner")),
    list(id = "smr_exp", methods = methods_smr_exp, results = results_smr_exp,
         requested = function(ctx) any(vapply(c("smr_expression_panel_psychencode", "smr_expression_panel_metabrain_basalganglia", "smr_expression_panel_metabrain_cerebellum", "smr_expression_panel_metabrain_cortex", "smr_expression_panel_metabrain_hippocampus", "smr_expression_panel_metabrain_spinalcord", "smr_expression_panel_eqtlgen"), function(f) cfg_on(ctx, f), logical(1)))),
    list(id = "smr_protein", methods = methods_smr_protein, results = results_smr_protein,
         requested = function(ctx) cfg_on(ctx, "smr_protein_panel_rosmap")),
    list(id = "tissue", methods = methods_tissue, results = results_tissue,
         requested = function(ctx) cfg_on(ctx, "tissue_magma")),
    list(id = "enrichment", methods = methods_enrichment, results = results_enrichment,
         requested = function(ctx) any(vapply(c("magma_drugtargetor", "twas_gsea_drugtargetor", "twas_gsea_drugtargetor_nondirectional", "gcsc", "twas_gsea_cmap", "magma_drugtargetor_gls", "drug_targetor_atc_gls"), function(f) cfg_on(ctx, f), logical(1)))),
    list(id = "pathway", methods = methods_pathway, results = results_pathway,
         requested = function(ctx) cfg_on(ctx, "magma_pathway") || cfg_on(ctx, "twas_gsea_pathway"))
  )
}

# ---- headers ----------------------------------------------------------------

writeup_banner <- function(kind, ctx, manifest) {
  ver    <- manifest$pipeline_tag %||% (manifest$pipeline_version %||% "unknown")
  schema <- ctx$schema_version %||% NA
  gstr(
    "> **AUTOMATICALLY GENERATED DRAFT — CHECK AND EDIT BEFORE USE.**\n",
    "> This {kind} text was produced by the GenoDisc write-up generator ",
    "(v{ctx$generator_version}) from the recorded run configuration, software ",
    "provenance and result objects (results schema v{schema}, pipeline {ver}). ",
    "It describes what the pipeline actually ran. It is a starting point, not a ",
    "finished section: verify every number, version and citation, fill the ",
    "`[CITE: ...]` placeholders, and edit for style before use. Citation keys ",
    "(`@key`) refer to `references.bib` in this folder."
  )}

.gwas_meta <- function(manifest, gw) {
  for (m in manifest$gwas) if (!is.null(m$name) && identical(as.character(m$name), gw)) return(m)
  list(name = gw)
}

# ---- reconciliation ---------------------------------------------------------

.reconcile <- function(ctx, specs) {
  missing <- character(0)
  for (s in specs) {
    if (identical(s$id, "preprocessing")) next
    if (isTRUE(s$requested(ctx)) && !(s$id %in% ctx$ran_blocks)) missing <- c(missing, s$id)
  }
  lines <- character(0)
  if (length(missing) > 0) {
    lines <- c(lines,
      "The following analyses were requested in the configuration but produced no results for any submitted GWAS, so they are not described above. This is expected when an analysis is restricted to certain ancestries (e.g. European-only steps run on a non-European GWAS); otherwise check the pipeline logs for failures:",
      paste0("- `", missing, "`"))
  }
  if (length(ctx$warnings) > 0) {
    lines <- c(lines, "", "Generator notes (verify before publication):",
               paste0("- ", ctx$warnings))
  }
  if (length(lines) == 0) return("")
  paste0("\n\n## Reconciliation — check before use\n\n", paste(lines, collapse = "\n"))
}

# ---- assembly ---------------------------------------------------------------

.assemble_methods <- function(ctx, paras, specs, manifest) {
  ordered <- character(0)
  for (s in specs) if (!is.null(paras[[s$id]])) ordered <- c(ordered, paras[[s$id]])
  body <- paste(ordered, collapse = "\n\n")
  closing <- methods_closing(ctx, manifest)
  paste0(writeup_banner("methods", ctx, manifest), "\n\n# Methods\n\n",
         body, "\n\n", closing)
}

.assemble_results <- function(ctx, gw, meta, paras, specs, manifest) {
  label <- meta$label %||% gw
  pop   <- meta$population %||% NA
  pop_clause <- if (!is.na(pop)) gstr(" ({pop}-ancestry)") else ""
  ordered <- character(0)
  for (s in specs) if (!is.null(paras[[s$id]])) ordered <- c(ordered, paras[[s$id]])
  body <- if (length(ordered)) paste(ordered, collapse = "\n\n") else "_No describable results were produced for this GWAS._"
  paste0(writeup_banner("results", ctx, manifest),
         gstr("\n\n# Results — {label}{pop_clause}\n\n"),
         body)
}

# ---- references -------------------------------------------------------------

# Parse references.bib into a named list key -> raw entry text (including @...{).
writeup_parse_bib <- function(pipeline_dir) {
  path <- file.path(pipeline_dir, "writeup", "references.bib")
  lines <- readLines(path, warn = FALSE)
  starts <- grep("^@[A-Za-z]+\\{", lines)
  if (length(starts) == 0) return(list())
  ends <- c(starts[-1] - 1L, length(lines))
  entries <- list()
  for (i in seq_along(starts)) {
    block <- lines[starts[i]:ends[i]]
    # trim trailing comment-only / blank lines between entries
    key <- sub(",$", "", sub("^@[A-Za-z]+\\{", "", regmatches(block[1], regexpr("^@[A-Za-z]+\\{[^,]+,", block[1]))))
    entries[[trimws(key)]] <- paste(block, collapse = "\n")
  }
  entries
}

writeup_subset_bib <- function(ctx, pipeline_dir) {
  entries <- writeup_parse_bib(pipeline_dir)
  keys <- ctx$bibkeys_used
  present <- keys[keys %in% names(entries)]
  missing <- setdiff(keys, names(entries))
  for (k in missing) writeup_warn(ctx, sprintf("Cited key '%s' has no entry in references.bib.", k))
  header <- paste0(
    "% AUTOMATICALLY GENERATED — the subset of the GenoDisc reference list cited by\n",
    "% this run's write-up. Source of truth: pipeline/writeup/references.bib.\n",
    "% Verify every entry (especially any marked VERIFY) before publication.\n\n")
  if (length(present) == 0) return(header)
  paste0(header, paste(vapply(sort(present), function(k) entries[[k]], character(1)), collapse = "\n\n"), "\n")
}

# Minimal BibTeX -> RIS conversion for the cited subset (trivial fields only).
writeup_bib_to_ris <- function(ctx, pipeline_dir) {
  entries <- writeup_parse_bib(pipeline_dir)
  keys <- sort(ctx$bibkeys_used[ctx$bibkeys_used %in% names(entries)])
  field <- function(txt, name) {
    m <- regmatches(txt, regexpr(paste0(name, "\\s*=\\s*\\{"), txt))
    if (length(m) == 0) return(NA_character_)
    # grab balanced-ish content: from after the opening brace to the matching line end brace
    start <- regexpr(paste0(name, "\\s*=\\s*\\{"), txt)
    rest <- substring(txt, start + attr(start, "match.length"))
    sub("\\}.*$", "", sub("\\},?\\s*\\n.*$", "", rest, perl = TRUE))
  }
  out <- character(0)
  for (k in keys) {
    txt <- entries[[k]]
    ty <- if (grepl("^@misc", txt)) "GEN" else "JOUR"
    out <- c(out, paste0("TY  - ", ty))
    auth <- field(txt, "author")
    if (!is.na(auth)) for (a in trimws(strsplit(auth, " and ", fixed = TRUE)[[1]])) out <- c(out, paste0("AU  - ", gsub("[{}]", "", a)))
    ti <- field(txt, "title");   if (!is.na(ti)) out <- c(out, paste0("TI  - ", gsub("[\\{\\}]", "", ti)))
    jo <- field(txt, "journal"); if (!is.na(jo)) out <- c(out, paste0("JO  - ", jo))
    py <- field(txt, "year");    if (!is.na(py)) out <- c(out, paste0("PY  - ", py))
    do <- field(txt, "doi");     if (!is.na(do)) out <- c(out, paste0("DO  - ", do))
    out <- c(out, paste0("ID  - ", k), "ER  - ", "")
  }
  paste0("Automatically generated RIS export of the cited references — verify before use.\n\n",
         paste(out, collapse = "\n"))
}

# ---- entry point ------------------------------------------------------------

generate_writeup <- function(output, merged_config, manifest, pipeline_dir,
                             opt = writeup_fmt_opts()) {
  ctx <- writeup_new_ctx(pipeline_dir, merged_config, opt,
                         schema_version = manifest$schema_version %||% NA)
  specs <- .writeup_specs()
  gwas_names <- setdiff(names(output), "configuration")

  # METHODS: one run-level paragraph per block, built from the first GWAS where it ran.
  methods_paras <- list()
  for (s in specs) {
    para <- NULL
    for (gw in gwas_names) {
      para <- s$methods(ctx, output[[gw]], .gwas_meta(manifest, gw))
      if (!is.null(para)) break
    }
    if (!is.null(para)) methods_paras[[s$id]] <- para
  }

  # RESULTS: per GWAS.
  results_files <- list()
  for (gw in gwas_names) {
    g <- output[[gw]]; meta <- .gwas_meta(manifest, gw)
    paras <- list()
    for (s in specs) {
      if (is.null(s$results)) next
      p <- s$results(ctx, g, meta)
      if (!is.null(p)) paras[[s$id]] <- p
    }
    results_files[[gw]] <- .assemble_results(ctx, gw, meta, paras, specs, manifest)
  }

  recon <- .reconcile(ctx, specs)
  methods_md <- paste0(.assemble_methods(ctx, methods_paras, specs, manifest), recon)

  files <- list()
  files[["methods.md"]] <- methods_md
  for (gw in gwas_names) files[[file.path(gw, "results.md")]] <- results_files[[gw]]
  files[["references.bib"]] <- writeup_subset_bib(ctx, pipeline_dir)
  files[["references.ris"]] <- writeup_bib_to_ris(ctx, pipeline_dir)

  list(
    files = files,
    meta = list(
      generator_version   = ctx$generator_version,
      schema_version_read = ctx$schema_version,
      bibkeys             = sort(ctx$bibkeys_used),
      ran_blocks          = sort(ctx$ran_blocks),
      warnings            = ctx$warnings,
      filenames           = names(files)
    )
  )
}
