pyramid <- function(width, height, res){
  trna <- read.delim("./rc/readcount_isodecoders.txt")
  out <- read.delim("./stat/stat_isoacceptor_list.txt")
  pval <- svalue(widget_list$pval)
  fc <- svalue(widget_list$FC)

  if(svalue(widget_list$fdr_s) == "Benjamini-Hochberg"){
    dw <- filter(out, logFC < -fc, out$Benjamini < pval)
    up <- filter(out, logFC > fc, out$Benjamini <  pval)
  }else if(svalue(widget_list$fdr_s) =="Bonferroni"){
    dw <- filter(out, logFC < -fc, out$Bonferroni<  pval)
    up <- filter(out, logFC > fc, out$Bonferroni <  pval)
  }else{
    dw <- filter(out, logFC < -fc, out$FDR <  pval)
    up <- filter(out, logFC > fc, out$FDR <  pval)
  }


  dw <- as.character(dw$Names)
  up <- as.character(up$Names)


  aa_codon <- function(name){
    out <- data.frame(do.call('rbind', strsplit(as.character(name), "-")))
    colnames(out) <- c("trna","aa")
    out <- data.frame(out <- table(out$aa))
    return(out)
  }

  u_gene <- aa_codon(up)
  colnames(u_gene) <- c("Var1","Up_DEtRNA")
  d_gene <- aa_codon(dw)
  colnames(d_gene) <- c("Var1", "Down_DEtRNA")
  geneplot <- merge(d_gene, u_gene, by="Var1", all= T)
  geneplot <- geneplot[c(order(geneplot$Up_DEtRNA,decreasing = T)),]
  png("./stat/plot/pyramid_isoaccepter%02d.png", width=width, height=height,  res=res)

  p <- function(){pyramid.plot(geneplot$Down_DEtRNA, geneplot$Up_DEtRNA,labels= geneplot$Var1,lxcol="#67A9CF", rxcol="#EF8A62",unit = "Freqency",gap=0.3, space=0.15, top.labels = c("Down_DEtRNAs", "tRNA-AA","Up_DEtRNAs"),laxlab=c(0,1,2,3), raxlab=c(0,1,2,3))}
  print(p())
  dev.off()
  save(p, geneplot, file="./stat/plot/Pyramid_Plot.RData")
}
