#' window
#'
#' @param ()
#' @return loadcl_name
#' @export
loadcl_name <- function(){
  cl_name <-  read.table(system.file("extdata", "class_name.txt", package = "tReasure",mustWork = TRUE), sep = "\t", fill = T,header = T, as.is = T)
  return(cl_name)
}

