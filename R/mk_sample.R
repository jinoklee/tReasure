#' window
#'
#' @param ()
#' @return mk_sample
#' @export
mk_sample <- function(h,...){
  dir <- svalue(fq_dir)
  if(identical(dir,character(0))){
    gmessage("Warning : Select a directory of FASTQ files")
  }else{
    fq_list <- list.files(svalue(fq_dir), pattern = "fastq", full=TRUE)
    if(length(grep("trim",fq_list)) == 0){
      sfile <- list.files(svalue(fq_dir), pattern = "fastq")
    }else{
      fq_list <- fq_list[-grep("trim", fq_list)]
      sfile <- list.files(svalue(fq_dir), pattern = "fastq")
      sfile <- sfile[-grep("trim", sfile)]
    }
    
    sname <- gsub(paste(".","fastq", sep = ""), "", sfile)
    sample <- data.frame(FileName = as.character(fq_list),
                         SampleName = as.character(sname),
                         Group = c("control", "test",rep("NA", (length(fq_list)-2))),
                         Batch = c(rep("NA", length(fq_list))))
    sample$Group <- as.factor(sample$Group)
    sample$Batch <- as.character(sample$Batch)
    tbl[] <- sample
    
    SampleFile <- file.path(dir, "sample.txt")
    write.table(sample, SampleFile, sep = "\t", quote = FALSE, row.names = FALSE)
    if(length(grep("control", sample$Group)) < 2 | length(grep("test", sample$Group)) < 2){
      insert(st, "Warnning : There should be at least 2 samples per group for differentially expression analysis.")}
    
    addHandlerChanged(tbl, handler = function(h, ...){
      new_sample <- tbl[]
      SampleFile <- file.path(svalue(fq_dir), "sample.txt")
      write.table(new_sample, SampleFile, sep = "\t", quote = FALSE, row.names = FALSE)
      insert(st, " ", do.newline = TRUE)
    })
    insert(st, ".", do.newline = TRUE)
  }
}