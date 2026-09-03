#!/usr/bin/env Rscript
# Inspect structure of results_package.rds

suppressMessages({
  library(data.table)
})

.script_dir <- local({
  args <- commandArgs(trailingOnly = FALSE)
  fa <- grep("^--file=", args, value = TRUE)
  if(length(fa) > 0) dirname(normalizePath(sub("^--file=", "", fa[1]))) else getwd()
})
source(file.path(.script_dir, "_paths.R"))
out_md <- file.path(report_dir, "results_object_structure.md")

x <- readRDS(rds_path)

cat("Top-level names:\n"); print(names(x))

traits <- setdiff(names(x), "configuration")
cat("\nTraits/GWAS:\n"); print(traits)

# Helper: format a value as a short type+dim string
short_desc <- function(v){
  if(is.null(v)) return("NULL")
  if(is.data.frame(v) || is.data.table(v))
    return(paste0("data.frame  ", nrow(v), " x ", ncol(v)))
  if(is.list(v)) return(paste0("list[", length(v), "]"))
  paste0(class(v)[1], "  len=", length(v))
}

# Walk tree and write a markdown structure summary
sink(out_md)
cat("# results_package.rds — structure summary\n\n")
cat("File: `", rds_path, "`\n\n", sep="")
cat("Top-level names: `", paste(names(x), collapse="`, `"), "`\n\n", sep="")
cat("Traits: ", paste(traits, collapse=", "), "\n\n")

cat("## Per-trait structure (showing first trait `", traits[1], "`)\n\n", sep="")
walk_list <- function(node, prefix=""){
  if(is.list(node) && !is.data.frame(node)){
    for(nm in names(node)){
      child <- node[[nm]]
      cat(prefix, "- **", nm, "**: ", short_desc(child), "\n", sep="")
      if(is.list(child) && !is.data.frame(child) && length(child) > 0){
        walk_list(child, paste0(prefix, "  "))
      } else if(is.data.frame(child) || is.data.table(child)){
        cat(prefix, "  - cols: ", paste(names(child), collapse=", "), "\n", sep="")
      }
    }
  }
}
walk_list(x[[traits[1]]])

cat("\n## TWAS-GSEA-related tables — column inventory across traits\n\n")

inspect_table <- function(tab, label){
  cat("### ", label, "\n\n", sep="")
  if(is.null(tab)){ cat("`NULL`\n\n"); return(invisible(NULL)) }
  cat("rows:", nrow(tab), "cols:", ncol(tab), "\n\n")
  cat("Columns: `", paste(names(tab), collapse="`, `"), "`\n\n", sep="")
  cat("First 3 rows:\n\n```\n")
  print(head(as.data.frame(tab), 3))
  cat("```\n\n")
}

for(tr in traits){
  cat("\n#### Trait: ", tr, "\n\n", sep="")
  paths <- list(
    `tx$drug$twas_gsea` = x[[tr]]$tx$drug$twas_gsea,
    `tx$drug$twas_gsea_nondir` = x[[tr]]$tx$drug$twas_gsea_nondir,
    `tx$atc$twas_gsea` = x[[tr]]$tx$atc$twas_gsea,
    `tx$atc$twas_gsea_nondir` = x[[tr]]$tx$atc$twas_gsea_nondir,
    `tx$drug$magma` = x[[tr]]$tx$drug$magma,
    `tx$atc$magma` = x[[tr]]$tx$atc$magma,
    `tx$cmap$drug` = x[[tr]]$tx$cmap$drug,
    `tx$cmap$moa`  = x[[tr]]$tx$cmap$moa
  )
  for(nm in names(paths)){
    tab <- paths[[nm]]
    if(is.null(tab)){
      cat("- **", nm, "**: NULL\n", sep="")
    } else {
      cat("- **", nm, "**: ", nrow(tab), " x ", ncol(tab),
          " — cols: `", paste(names(tab), collapse="`, `"), "`",
          if("Panel" %in% names(tab)) paste0(" — Panels: ", length(unique(tab$Panel))) else "",
          "\n", sep="")
    }
  }
}

cat("\n## Configuration\n\n")
if(!is.null(x$configuration)){
  cat("Repo remote: `", x$configuration$repo$remote, "`\n\n", sep="")
  cat("Repo branch/commit: `", x$configuration$repo$branch, "` / `", x$configuration$repo$commit, "`\n\n", sep="")
  if(!is.null(x$configuration$gwas_list)){
    cat("\nGWAS list:\n\n```\n")
    print(x$configuration$gwas_list)
    cat("```\n")
  }
}

sink()
cat("Wrote: ", out_md, "\n")
