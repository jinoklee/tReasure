#' window
#'
#' @param ()
#' @return edger
#' @export

edger<- function(count, edgeRmethod){
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
  
  y <- DGEList(count[,-1], group = case)
  y <- calcNormFactors(y)
  y <- estimateDisp(y, design, robust = TRUE)
  
  if( edgeRmethod== "Exact test"){
    test <- exactTest(y)
  }else if(edgeRmethod == "Quasi-likelihood F-test"){
    fit <- glmQLFit(y, design)
    test <- glmQLFTest(fit, coef=coef)
  }else{
    fit <- glmFit(y, design)
    test <- glmLRT(fit, coef=coef)
  }
  
  out_bh <- topTags(test, n = Inf, adjust.method="BH")
  out_bh <- out_bh$table
  colnames(out_bh)[ncol(out_bh)] <- "Benjamini"
  
  out_bonf <- topTags(test, n= Inf, adjust.method = "bonferroni")
  out_bonf <- out_bonf$table
  colnames(out_bonf)[ncol(out_bonf)] <- "Bonferroni"
  
  out_fdr <- topTags(test, n= Inf, adjust.method = "fdr")
  out_fdr <- out_fdr$table
  colnames(out_fdr)[ncol(out_fdr)] <- "FDR"
  
  out <- cbind(out_fdr[,], Bonferroni = out_bonf$Bonferroni, Benjamini=out_bh$Benjamini)
  out <- data.frame(Names = rownames(out), out[,])
  return(out)
}