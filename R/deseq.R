#' window
#'
#' @param ()
#' @return deseq
#' @export

deseq<- function(count){
  rownames(count) <- count$Names
  sFile<- read.delim("sample.txt")
  colData  <- sFile[,c("Group","Batch")]
  colnames(colData) <- c("condition","type")
  colData$condition <- factor(colData$condition)
  colData$type <- factor(colData$type)
  
  dds <-DESeqDataSetFromMatrix(countData = count[,-1], colData = colData , design = ~ condition)
  
  dds$condition <- factor(dds$condition, levels=c("control","test"))
  dds <- DESeq(dds)
  
  res_bh <- results(dds, contrast = c("condition","test","control"),pAdjustMethod= "BH")
  res_bh <- as.data.frame(res_bh)
  colnames(res_bh)[ncol(res_bh)] <- "Benjamini"
  
  res_bonf <- results(dds, contrast = c("condition","test","control"),pAdjustMethod= "bonferroni")
  res_bonf<- as.data.frame(res_bonf)
  colnames(res_bonf)[ncol(res_bonf)] <- "Bonferroni"
  
  res_fdr <- results(dds, contrast = c("condition","test","control"),pAdjustMethod= "fdr")
  res_fdr <- as.data.frame(res_fdr)
  colnames(res_fdr)[ncol(res_fdr)] <- "FDR"
  
  res <- cbind(res_fdr[,], Bonferroni = res_bonf$Bonferroni, Benjamini=res_bh$Benjamini)
  colnames(res)[2] <- "logFC"
  res <- data.frame(Names = rownames(res), res[,])
  return(res)
}