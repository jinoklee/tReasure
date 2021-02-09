f <- function(h,...){
  setwd(svalue(fq_dir))
  dir <- getwd()
  if(!dir.exists(paste(dir, "/trim", sep = ""))){
    dir.create(paste(dir, "/trim", sep = ""), recursive = TRUE)}
  if(!dir.exists(paste(dir, "/stat/plot", sep = ""))){
    dir.create(paste(dir, "/stat/plot", sep = ""), recursive = TRUE)}
}
