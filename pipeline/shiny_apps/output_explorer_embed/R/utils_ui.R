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
