#' window
#'
#' @param ()
#' @return p2Nms
#' @export
p2Nms <- function(d, refv, envir=.GlobalEnv){
  cl_name <-  read.table(system.file("extdata", "class_name.txt", package = "tReasure",mustWork = TRUE), sep = "\t", fill = T,header = T, as.is = T)
  unique(cl_name$P2[grep(refv, cl_name$P1)], envir=envir) 
}  
