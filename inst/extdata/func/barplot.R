barplot <- function(width, height, res){
  trna <- read.delim("./rc/readcount_isodecoders.txt")
  out <- read.delim("./stat/stat_isodecoder_list.txt")

  pval <- svalue(widget_list$pval)
  fc <- svalue(widget_list$FC)

  if(svalue(widget_list$fdr_s) == "Benjamini-Hochberg"){
    dw <- filter(out, logFC < -pval, out$Benjamini < pval)
    up <- filter(out, logFC > fc, out$Benjamini < pval)
  }else if(svalue(widget_list$fdr_s) =="Bonferroni"){
    dw <- filter(out, logFC < -fc, out$Bonferroni< pval)
    up <- filter(out, logFC > fc, out$Bonferroni < pval)
  }else{
    dw <- filter(out, logFC < -fc, out$FDR < pval)
    up <- filter(out, logFC > fc, out$FDR < pval)
  }

  # dw <- filter(out, logFC < -1.5, out$Benjamini < 0.5)
  # up <- filter(out, logFC > 1.5, out$Benjamini < 0.5)

  dw <- as.character(dw$Names)
  up <- as.character(up$Names)
  no <- as.character(out$Names[!(out$Names %in% c(dw,up))])
  su <- c(up, dw,no)
  #rm <- as.character(trna$Names[!(trna$Names %in% su)])

  tRNA_aa_codon <- function(data, t){
    out <- data.frame(aa_codon=substr(data,6,13))
    out$aa_codon <- gsub("-$","", out$aa_codon)
    out <- data.frame(table(out$aa_codon))
    outt <- data.frame(do.call('rbind', strsplit(as.character(out$Var1), split = "-")))
    colnames(outt) <- c("aa","codon")
    out <- cbind(out,outt)
    out$Sig <- rep(t,nrow(out))
    return(out)
  }


  up <- tRNA_aa_codon(up, "Up_DEtRNA")
  dw <- tRNA_aa_codon(dw, "Down_DEtRNA")
  no <- tRNA_aa_codon(no, "Non_DEtRNA")
  #rm <- tRNA_aa_codon(rm, "filter")
  #rm$Freq <- rep(0, nrow(rm))
  c <- rbind(up,dw, no)
  c <- c %>% mutate_if(is.factor, as.character)
  c<- c[order(c$aa),]
  n <- length(grep("iMet", c$aa))
  c$codon[grep("iMet", c$aa)] <- rep("iCAT", n)
  c$Var1 <- gsub("iMet-CAT","iMet-iCAT", c$Var1)
  empty_bar <- 1

  to_add <- data.frame(matrix(NA, empty_bar*nlevels(as.factor(c$aa)), ncol(c)))
  colnames(to_add) <- colnames(c)
  to_add$aa <- rep(levels(as.factor(c$aa)), each=empty_bar)
  to_add$codon <- paste("a", seq(1, nrow(to_add)))
  to_add$Var1 <- paste("a", seq(1, nrow(to_add)))
  c <- rbind(c, to_add)


  c$aa <- as.character(c$aa)
  c$codon <- as.character(c$codon)

  c <- c %>% arrange(aa,codon)

  c$Freq[is.na(c$Freq)] <- 0

  base <- data.frame(Var1=c$Var1)
  base$Var1 <- as.character(base$Var1)
  base_data <- data.frame(Var1= base[!duplicated(base$Var1),])
  out <- data.frame(do.call('rbind',strsplit(as.character(base_data$Var1), split="-")))

  base_data$aa <- out$X1
  base_data$codon <- out$X2
  base_data$id <- seq(1, nrow(base_data))

  b1_data<- base_data%>%
    group_by(aa) %>%
    summarize(start=min(id), end=max(id)) %>%
    rowwise() %>%
    mutate(title=mean(c(start, end)))

  b1_data<- b1_data[!grepl("^a", b1_data$aa),]

  b2_data <- base_data%>%
    group_by(codon) %>%
    summarize(start=min(id), end=max(id)) %>%
    rowwise() %>%
    mutate(title=mean(c(start, end)))

  b2_data<- b2_data[!grepl("^a", b2_data$codon),]

  png("./stat/plot/barplot_isodecoder%02d.png", height=height, width=width, res=res)
  p <- function(){
    ggplot(c, aes(x=codon, y=Freq, fill=Sig))+
      geom_bar( stat="identity")+
      scale_fill_manual(values=c("Up_DEtRNA"="#EF8A62", "Down_DEtRNA"="#67A9CF", "Non_DEtRNA"="grey89", filter="white"), breaks = c("Up_DEtRNA","Non_DEtRNA","Down_DEtRNA"))+
      scale_x_discrete(limits=unique(c$codon))+
      theme_minimal()+
      scale_y_continuous(breaks = seq(0,max(c$Freq),2))+
      theme(axis.text = element_blank(),
            panel.grid.major.x = element_blank(),
            legend.position = "top")+
      ylab("Frequency")+
      xlab("Aminoacid : Anticodon")+
      #annotate("text",x=rep(1, max(c$Freq)/2), y=seq(2,max(c$Freq),2), label=c("2","4","6","8","10","12","14"), color="black", size = 3)+
      geom_text(data=b1_data, aes(x = title, y = -1.5, label=aa), colour ="black", size=3, fontface="bold", inherit.aes = FALSE)+
      geom_segment(data=b1_data, aes(x=start, y = -1.1, xend=end, yend=-1.1), colour="black", alpha=1, size=1,inherit.aes = FALSE)+
      geom_text(data=b2_data, aes(x = title, y= -0.45, label=codon), size=3, colour = "black", angle=90, inherit.aes = FALSE)}
  print(p())
  dev.off()
  #ggsave("./stat/plot/barplot_isodecoder%02d.png")
  save(p, pval,fc, c, file="./stat/plot/Bar_Plot.RData")
}
