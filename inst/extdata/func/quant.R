quant<- function(count){
  rownames(count) <- count$Names
  sFile<- read.delim("sample.txt")
  case <- factor(sFile$Group)
  case <- relevel(case, ref="control")
  if(nlevels(sFile$Batch) != 0){batch <- factor(sFile$Batch)
  design <- model.matrix(~batch+case)
  coef = nlevels(batch)+1
  }else{design <- model.matrix(~case)
  coef = 2}
  rownames(design) <- colnames(count)[-1]
  
  v <- voom(count[,-1], design, normalize = "quantile" )
  fit <- lmFit(v, design)
  fit <- eBayes(fit)
  
  out_bh <- topTable(fit, coef=ncol(design), n= Inf, adjust.method = "BH")
  colnames(out_bh)[grep("adj",colnames(out_bh))] <- "Benjamini"
  
  out_bonf <- topTable(fit, coef=ncol(design), n= Inf, adjust.method = "bonferroni")
  colnames(out_bonf)[grep("adj",colnames(out_bonf))] <- "Bonferroni"
  
  out_fdr <- topTable(fit, coef=ncol(design), n= Inf, adjust.method = "fdr")
  colnames(out_fdr)[grep("adj",colnames(out_fdr))] <- "FDR"
  
  out <- cbind(out_fdr[,c(1,2,3,4,6,5)], Bonferroni = out_bonf$Bonferroni, Benjamini=out_bh$Benjamini)
  out <- data.frame(Names = rownames(out), out[,])
  return(out)
}