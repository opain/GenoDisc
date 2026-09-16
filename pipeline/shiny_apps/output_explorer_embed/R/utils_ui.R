#' Interpretive legend block shown under a table or figure.
#'
#' Renders a small grey "What these mean" heading followed by a bulleted list of
#' term -> explanation pairs, styled via the `.gd-legend` CSS class in app.R.
#'
#' @param items Named list / character vector mapping term -> explanation.
#' @param heading Heading shown above the list.
#' @return A shiny tag (div).
gd_legend <- function(items, heading = "What these mean") {
  lis <- lapply(seq_along(items), function(i) {
    tags$li(tags$strong(names(items)[i]), " — ", items[[i]])
  })
  tags$div(
    class = "gd-legend",
    tags$strong(heading),
    tags$ul(lis)
  )
}

# ============================================================================
# Download-size unit handling
#
# Shared helpers for the "Units: Inches / Centimetres / Pixels" selector that
# every plot download-options panel exposes. All three helpers are
# module-agnostic: each caller supplies its own input IDs so the existing
# per-module naming conventions (prefix-based, suffix-based, plain) survive.
# ============================================================================

#' Choice list for the Units selectInput. Displayed labels in full, values
#' compatible with grDevices' `units` argument.
.dl_units_choices <- c("Inches" = "in",
                        "Centimetres" = "cm",
                        "Pixels" = "px")

#' Convert a plot dimension between "in", "cm", and "px".
#'
#' `dpi` is the resolution used for px conversions. Callers that back a
#' user-picked DPI slider should pass that value so switching between inches
#' and px keeps the output file the same size at the current DPI (e.g. 12 in
#' at 300 dpi -> 3600 px, not 1152 px at a fixed 96 dpi).
#' Returns `NA_real_` for non-finite or missing input; passes `x` through
#' unchanged when `from == to`.
.convert_dim <- function(x, from, to, dpi = 96) {
  x <- suppressWarnings(as.numeric(x))
  if (length(x) != 1L || !is.finite(x)) return(NA_real_)
  if (identical(from, to)) return(x)
  in_val <- switch(from, "in" = x, "cm" = x / 2.54, "px" = x / dpi, NA_real_)
  switch(to,   "in" = in_val, "cm" = in_val * 2.54, "px" = in_val * dpi, NA_real_)
}

#' Register the auto-convert observer for a Units / Width / Height triple.
#'
#' Call once per download-options block, after building the widgets. The
#' observer watches `input[[units_id]]` and, on each change, converts the
#' width and height values so the OUTPUT FILE size the user set keeps its
#' meaning across the switch (8 in <-> 20.32 cm; 12 in at 300 dpi <-> 3600 px).
#'
#' @param input,session Shiny session objects from `moduleServer`.
#' @param units_id,width_id,height_id Input ids (module-scoped, no `ns()`).
#' @param prev_val A `reactiveVal(initial_units)` the caller creates and
#'   passes so the observer can tell the previous unit apart from the new
#'   one. Reset it if you programmatically change units.
#' @param dpi_id Optional input id of a numericInput carrying the current
#'   DPI. When supplied, px conversions use that DPI; otherwise falls back
#'   to `default_dpi`.
#' @param default_dpi Fallback DPI used when `dpi_id` is NULL or the input
#'   value is empty. 300 matches the module defaults.
dl_units_auto_convert <- function(input, session,
                                   units_id, width_id, height_id,
                                   prev_val,
                                   dpi_id = NULL, default_dpi = 300) {
  observeEvent(input[[units_id]], {
    new_u <- input[[units_id]]
    old_u <- prev_val()
    if (identical(new_u, old_u)) return()
    dpi <- if (!is.null(dpi_id)) (input[[dpi_id]] %||% default_dpi) else default_dpi
    w <- .convert_dim(input[[width_id]],  from = old_u, to = new_u, dpi = dpi)
    h <- .convert_dim(input[[height_id]], from = old_u, to = new_u, dpi = dpi)
    rd <- if (identical(new_u, "px")) 0 else 2
    if (is.finite(w)) updateNumericInput(session, width_id,  value = round(w, rd))
    if (is.finite(h)) updateNumericInput(session, height_id, value = round(h, rd))
    prev_val(new_u)
  }, ignoreInit = TRUE)
}

#' Open the appropriate `grDevices` device for a given format and units.
#'
#' PNG honours the user's units natively (so "800 px" gives an 800-px file
#' rather than being round-tripped through inches). PDF / SVG only take
#' inches, so px / cm values are converted at the current DPI so a px value
#' that was auto-converted from inches at that DPI round-trips back to the
#' same inch value.
#'
#' Caller is responsible for `print(p)` / `grid.draw(gt)` and
#' `grDevices::dev.off()`.
open_plot_device <- function(file, fmt, w, h, units = "in", dpi = 300) {
  if (identical(fmt, "png")) {
    grDevices::png(file, width = w, height = h, units = units, res = dpi)
  } else {
    w_in <- .convert_dim(w, from = units, to = "in", dpi = dpi)
    h_in <- .convert_dim(h, from = units, to = "in", dpi = dpi)
    switch(fmt,
      pdf = grDevices::pdf(file, width = w_in, height = h_in),
      svg = grDevices::svg(file, width = w_in, height = h_in)
    )
  }
}
