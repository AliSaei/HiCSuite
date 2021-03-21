#' Run scaffoldR app
#'
#' Launch the app for scaffolding contigs
#' @export
run_scaffoldR <- function() {
  appDir <- system.file("shiny-apps", "scaffoldR", package = "HiCSuite")
  if (appDir == "") {
    stop("Could not find app directory. Try re-installing `HiCSuite`.", call. = FALSE)
  }
  
  shiny::runApp(appDir, display.mode = "normal")
}