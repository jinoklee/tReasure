#' window
#'
#' @param ()
#' @return p3Nms
#' @export
p3Nms <- function(d, refv,envir=.GlobalEnv){
  cl_name <-  read.table(system.file("extdata", "class_name.txt", package = "tReasure",mustWork = TRUE), sep = "\t", fill = T,header = T, as.is = T)
  unique(cl_name$P3[grep(refv, cl_name$P2)], envir=envir)
}  