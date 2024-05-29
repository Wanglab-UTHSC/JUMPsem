#' @export quiet
#' @param x The sentence needs to be quiet.

quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}
