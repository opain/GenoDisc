#!/usr/bin/Rscript

# Number / P-value / string formatting for the automated write-up generator.
# Pure, dependency-free (base R). Every quantitative value that reaches the
# generated text passes through one of these so sig-figs and P-value style are
# consistent across blocks and configurable in one place.

# glue -> character in the caller's environment. Used instead of `glue(...) |>
# as.character()` because the pipeline runs R 4.0.2, which predates the `|>` pipe.
gstr <- function(...) as.character(glue::glue(..., .envir = parent.frame()))

# Default formatting options. Override by passing a modified list as `opt` to the
# formatters (the write-up driver threads a single opt list through).
writeup_fmt_opts <- function(sig = 3L, p_style = c("sci", "plain"),
                             p_sci_below = 1e-4, p_digits = 2L) {
  list(
    sig        = as.integer(sig),      # significant figures for generic numbers
    p_style    = match.arg(p_style),   # "sci" -> scientific for small P; "plain" -> fixed decimals
    p_sci_below = p_sci_below,         # P below this is rendered in scientific notation
    p_digits   = as.integer(p_digits)  # sig figs inside the P-value mantissa / decimals
  )
}

# Generic number: `sig` significant figures, NA-safe.
fmt_num <- function(x, opt = writeup_fmt_opts(), sig = NULL) {
  if (length(x) == 0 || is.null(x)) return("NA")
  s <- if (is.null(sig)) opt$sig else as.integer(sig)
  vapply(x, function(v) {
    if (is.na(v)) return("NA")
    if (!is.finite(v)) return(as.character(v))
    trimws(formatC(signif(v, s), format = "fg", digits = s, big.mark = ""))
  }, character(1))
}

# Integer with thousands separators, NA-safe.
fmt_int <- function(x) {
  if (length(x) == 0 || is.null(x) || is.na(x)) return("NA")
  formatC(round(as.numeric(x)), format = "d", big.mark = ",")
}

# P-value, NA-safe. Scientific for small P (default), otherwise fixed decimals.
fmt_p <- function(p, opt = writeup_fmt_opts()) {
  if (length(p) == 0 || is.null(p) || is.na(p)) return("NA")
  p <- as.numeric(p)
  if (!is.finite(p)) return(as.character(p))
  if (p == 0) return("< 1e-300")
  if (opt$p_style == "sci" && p < opt$p_sci_below) {
    return(formatC(signif(p, opt$p_digits), format = "e", digits = opt$p_digits - 1L))
  }
  # plain: round to a sensible number of decimals (never show 0.000 for a tiny P)
  dec <- max(opt$p_digits, 3L)
  out <- formatC(signif(p, opt$p_digits), format = "f", digits = dec)
  if (as.numeric(out) == 0) {
    return(formatC(signif(p, opt$p_digits), format = "e", digits = opt$p_digits - 1L))
  }
  out
}

# Estimate (SE) pair, e.g. "0.066 (0.009)".
fmt_est_se <- function(est, se, opt = writeup_fmt_opts()) {
  paste0(fmt_num(est, opt), " (", fmt_num(se, opt), ")")
}

# Interval "lo-hi".
fmt_ci <- function(lo, hi, opt = writeup_fmt_opts()) {
  paste0(fmt_num(lo, opt), "–", fmt_num(hi, opt))
}

# Percentage from a fraction (0-1) or an already-percent value if `frac=FALSE`.
fmt_pct <- function(x, opt = writeup_fmt_opts(), frac = TRUE) {
  if (length(x) == 0 || is.null(x) || is.na(x)) return("NA")
  v <- if (frac) as.numeric(x) * 100 else as.numeric(x)
  paste0(fmt_num(v, opt), "%")
}

# English pluralisation helper: pluralise("gene", 1) -> "gene"; (n=2) -> "genes".
pluralise <- function(word, n, plural = NULL) {
  if (!is.na(n) && n == 1) return(word)
  if (!is.null(plural)) return(plural)
  paste0(word, "s")
}

# "was"/"were" agreement.
was_were <- function(n) if (!is.na(n) && n == 1) "was" else "were"

# Join a character vector into an English list: "a", "a and b", "a, b and c".
oxford <- function(x, conj = "and") {
  x <- x[!is.na(x) & nzchar(x)]
  n <- length(x)
  if (n == 0) return("")
  if (n == 1) return(x[1])
  if (n == 2) return(paste(x[1], conj, x[2]))
  paste0(paste(x[-n], collapse = ", "), " ", conj, " ", x[n])
}

# Safe nested access: pluck(output, "mol_assoc", "exp", "fusion"). Returns NULL if
# any level is missing, so block builders never error on an absent branch.
pluck <- function(x, ...) {
  keys <- c(...)
  for (k in keys) {
    if (is.null(x)) return(NULL)
    if (!is.list(x) || is.null(x[[k]])) return(NULL)
    x <- x[[k]]
  }
  x
}

# Row count that mirrors manifest's .gd_row_count: nrow of a data.frame, or
# nrow of $res / $results / $table inside a list. NA when uncountable.
wu_nrow <- function(x) {
  if (is.null(x)) return(NA_integer_)
  if (is.data.frame(x)) return(nrow(x))
  if (is.list(x)) {
    for (k in c("res", "results", "table")) {
      if (!is.null(x[[k]]) && is.data.frame(x[[k]])) return(nrow(x[[k]]))
    }
  }
  NA_integer_
}

# Count rows of `df` whose column `col` is non-NA and < thr (the pipeline/viewer
# FDR rule). Returns 0L for an empty/absent table.
count_below <- function(df, col, thr = 0.05) {
  if (is.null(df) || !is.data.frame(df) || !col %in% names(df) || nrow(df) == 0) return(0L)
  v <- suppressWarnings(as.numeric(df[[col]]))
  sum(!is.na(v) & v < thr)
}
