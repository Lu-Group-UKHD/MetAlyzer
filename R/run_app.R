#' Launch the Shiny App
#'
#' This function launches the Shiny application.
#'
#' @export
run_app <- function() {
  app_dir <- system.file("shinyapp", package = "MetAlyzer")
  if (app_dir == "") {
    stop("Could not find app directory. Try re-installing the package.", call. = FALSE)
  }
  shiny::runApp(app_dir, launch.browser = TRUE)
}
