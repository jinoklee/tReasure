#' window
#'
#' @param ()
#' @return volcano
#' @export
volcano <- function(width, height, res, pval, fc, statmethod){
  out <- read.delim("./stat/stat_trna_list.txt")
  
  if(statmethod == "Bejamini-Hochberg"){
    outt <- select(out, Names, logFC, Benjamini)
    colnames(outt)[ncol(outt)]<-"FDR"
  }else if(statmethod == "Bonferroni"){
    outt <- select(out, Names, logFC, Bonferroni)
    colnames(outt)[ncol(outt)]<-"FDR"
  }else{
    outt <- select(out, Names, logFC, FDR)
  }
  

  detRNA <- mutate(outt, Sig=ifelse(outt$FDR<= pval & outt$logFC <= -fc, "Down_DEtRNA", ifelse(outt$FDR <= pval & outt$logFC>= fc,"Up_DEtRNA", "Non_DEtRNA")))
  
  png("./stat/plot/volcanoplot_trna%02d.png", height=height, width=width, res=res)
  p <- function(){
    ggplot(detRNA, aes(x = logFC, y = -log10(FDR)))+
      geom_point(size=1.5, aes(col=Sig)) +
      xlab(" log2 Fold Change") +
      ylab("-log10 Adjusted P value ") +
      geom_vline(xintercept = c(-fc,fc),col = "red",linetype = "dotted",size = 0.5) +
      geom_hline(yintercept = c(-log10(pval)),col = "red", linetype = "dotted",size = 0.5) +
      theme_classic()+
      theme(legend.position = "top")+
      scale_colour_manual(values = c("Non_DEtRNA"="grey89", "Up_DEtRNA"="tomato","Down_DEtRNA"="#67A9CF"))}
  print(p())
  dev.off()
  save(p,pval,fc,detRNA,file="./stat/plot/Volcano_Plot.RData")
  
}