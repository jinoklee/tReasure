#' window
#'
#' @param ()
#' @return tReasure function
#' @export
#' stopFuture
stopFuture <- function(x){
  tools::pskill(x,signal = tools::SIGTERM)
  tools::pskill(x,signal = tools::SIGKILL)
}