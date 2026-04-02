#' Generate JS callback for scientific notation formatting
#'
#' @param col_indices Integer vector of 0-based column indices to format
#' @return Character vector containing the JS rowCallback function
sci_notation_js <- function(col_indices) {
  body_lines <- character(0)
  for (idx in col_indices) {
    var_name <- paste0("v", idx)
    body_lines <- c(body_lines,
                    paste0("  var ", var_name, " = data[", idx, "];"),
                    paste0("  $('td:eq(", idx, ")', row).html(", var_name, ".toExponential(2));"))
  }
  c(
    "function(row, data, displayNum, index){",
    body_lines,
    "}"
  )
}

#' Render a datatable with standard options used throughout the app
#'
#' @param data Data frame or data.table to display
#' @param sci_cols Integer vector of 0-based column indices for scientific notation, or NULL
#' @param center_targets Column targets for centering (default "_all")
#' @param width_specs Named list of width specifications, e.g. list("60px" = 4:5). NULL for none.
#' @param escape Passed to datatable(). TRUE by default.
#' @param ... Additional arguments passed to datatable()
#' @return A DT datatable object
render_standard_dt <- function(data, sci_cols = NULL, center_targets = "_all",
                               width_specs = NULL, escape = TRUE, ...) {

  opts <- list()

  if (!is.null(sci_cols)) {
    js <- sci_notation_js(sci_cols)
    opts$rowCallback <- DT::JS(js)
  }

  col_defs <- list(list(className = 'dt-center', targets = center_targets))

  if (!is.null(width_specs)) {
    for (w in names(width_specs)) {
      col_defs <- c(col_defs, list(list(width = w, targets = width_specs[[w]])))
    }
  }

  opts$columnDefs <- col_defs

  DT::datatable(data, rownames = FALSE, options = opts, escape = escape, ...)
}
