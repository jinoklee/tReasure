anl_filter <- function(h,...){
  if(identical(svalue(selrc_button),character(0))){
    dir <- getwd()
  }else{
    setwd(svalue(selrc_button))
    dir <- getwd()
    if(!dir.exists(paste(dir, "/pre", sep = ""))){
      dir.create(paste(dir, "/pre", sep = ""), recursive = TRUE)}
    if(!dir.exists(paste(dir, "/stat/plot", sep = ""))){
      dir.create(paste(dir, "/stat/plot", sep = ""), recursive = TRUE)}
  }

  sFile <- read.delim("sample.txt")
  tcount <- read.table(file.path(dir, "rc","readcount_trnas.txt"), header = T)
  ccount <- read.table(file.path(dir, "rc","readcount_isodecoders.txt"), header = T)
  acount <- read.table(file.path(dir, "rc","readcount_isoacceptors.txt"), header = T)

  filter_gene <- function(count){
    rownames(count) <- count$Names
    count <- count[,colnames(count) %in% sFile$SampleName]
    x <- DGEList(count, group = sFile$Group)
    isexpr <- rowSums(cpm(x) > svalue(cfilv)) >= svalue(sfil)*nrow(sFile)/100
    x <- x[isexpr,]
    return(x)
  }


  filter_t <- filter_gene(tcount)
  filter_c <- filter_gene(ccount)
  filter_a <- filter_gene(acount)

  filter_tt <- data.frame(Names = rownames(filter_t), filter_t$counts)
  filter_cc <- data.frame(Names = rownames(filter_c), filter_c$counts)
  filter_aa <- data.frame(Names = rownames(filter_a), filter_a$counts)


  f_kep[] <- data.frame(No.reads=c("total", "keep", "remove"),
                        tRNAs=c(nrow(tcount), nrow(filter_tt), nrow(tcount)-nrow(filter_tt)),
                        isodecoders=c(nrow(ccount), nrow(filter_cc), nrow(ccount)-nrow(filter_cc)),
                        isoacceptors=c(nrow(acount), nrow(filter_aa), nrow(acount)-nrow(filter_aa)))



  write.table(filter_tt, file.path(dir, "rc","filtered_readcount_trnas.txt"),sep = "\t", quote = F, row.names = F)
  write.table(filter_cc, file.path(dir, "rc","filtered_readcount_isodecoders.txt"),sep = "\t", quote = F, row.names = F)
  write.table(filter_aa, file.path(dir, "rc","filtered_readcount_isoacceptors.txt"),sep = "\t", quote = F, row.names = F)


  y <- calcNormFactors(filter_t)
  cols <- as.character(sFile$Group)
  cols[cols=="control"] = "blue"
  cols[cols=="test"] ="red"
  png( './stat/plot/plotMDS.png', width=500, height=500,  res=65)
  p <- function(){
    plotMDS(y,col=cols)
  }
  print(p())
  dev.off()
  save(p,y,cols, file = "./stat/plot/MDS_plot.RData")
  insert(st, " ", do.newline = TRUE)
  insert(st,"Done : Filtering. Click! Next tab of DEtRNA Detection.", do.newline = TRUE )
  insert(st, ".", do.newline = TRUE)
}
