#' window
#'
#' @param ()
#' @return sel_sample
#' @export

sel_sample <-  function(h,...){
  sample2 <- read.table(svalue(sel_button), header = T, sep = "\t")
  sampath <- gsub("sample.txt", "", svalue(sel_button))
  setwd(sampath)
  dir <- getwd()
  if(!dir.exists(paste(dir, "/pre", sep = ""))){
    dir.create(paste(dir, "/pre", sep = ""), recursive = TRUE)}
  if(!dir.exists(paste(dir, "/stat/plot", sep = ""))){
    dir.create(paste(dir, "/stat/plot", sep = ""), recursive = TRUE)}
  
  sample2$Batch <- as.character(sample2$Batch)
  tbl[] <- sample2
  insert(st, " ", do.newline = TRUE)
  insert(st, "Complete uploading the sample list.", do.newline = TRUE)
  if(length(grep("control", sample2$Group)) < 2 | length(grep("test", sample2$Group)) < 2){
    insert(st, "Warnning : There should be at least 2 samples per group for differentially expression analysis.")}
  # save update sample list
  addHandlerChanged(tbl, handler = function(h, ...){
    new_sample <- tbl[]
    SampleFile <- file.path(svalue(fq_dir), "sample.txt")
    write.table(new_sample, SampleFile, sep = "\t", quote = FALSE, row.names = FALSE)
    insert(st, " ", do.newline = TRUE)
    insert(st, ".", do.newline = TRUE)
  })
}
