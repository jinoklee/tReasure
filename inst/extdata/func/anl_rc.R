anl_rc <- function(ref_name){
  dir <- getwd()
  workDir <- file.path(dir, "post")
  bamlist <- list.files(workDir, pattern = "\\.bam$", recursive = T)
  bamlist <- basename(bamlist)

  rc_df <- function(dfbn){
    chr <- data.frame(stringr::str_split_fixed(dfbn$seqnames, "[.]", 2))
    out <- data.frame(stringr::str_split_fixed(chr$X2, "_", 4))
    out$X1 <- paste(chr$X1, out$X1,sep = ".")
    colnames(out) <- c("tRNAscan.SE.id","GtRNAdb.id","Iso","Confidence")
    dfbn <- cbind(out, dfbn[,-1])
    return(dfbn)
  }

  dfbn1 <- idxstatsBam(file.path(workDir,bamlist[1]))
  dfbn1 <- dfbn1[,c("seqnames"), drop = F]

  for ( bn in bamlist){
    df <- idxstatsBam(file.path(workDir,bn))
    df <- df[,c("seqnames", "mapped")]
    colnames(df)[2] <- gsub("\\_.*","",bn)
    dfbn1 <- full_join(dfbn1, df)
  }

  dfbn <- rc_df(dfbn1)
  write.table(dfbn, file.path(dir,"rc", "readcount_all_tRNAs.txt"), sep = "\t", row.names = F, col.names = T, quote = F)

  hi <- filter(dfbn, Confidence =="high")
  tcount <- hi[,-c(1,3,4)]
  colnames(tcount)[1] <- "Name"

  write.table(tcount,file.path(dir,"rc","readcount_trnas.txt"), sep = "\t", row.names = F, quote = F)

  gtable(tcount, container=RC_trna)
  insert(st, " ", do.newline = TRUE)
  insert(st,"Done : Counting of individual tRNAs", do.newline = TRUE )


  ## isodecoder
  clu<- c()
  for (i in unique(hi$Iso)){
    df <- hi[grep(i, hi$Iso), -c(1:4)]
    if(nrow(df) > 1){
      df <- colSums(df)
    }
    df$Name <- i
    df <- data.frame(df)
    clu <- rbind(clu, df)
  }
  ccount <- clu[,c(ncol(clu), seq(1,ncol(clu)-1))]
  write.table(ccount ,file.path(dir,"rc", "readcount_isodecoders.txt"), sep = "\t", row.names = F, quote = F)
  gtable(ccount, container=RC_codon)
  insert(st, " ", do.newline = TRUE)
  insert(st,"Done : Counting of isodecoders", do.newline = TRUE )


  ## isoacceptor
  acount <- data.frame(Name=substr(ccount$Name, 1,12), ccount[,-1])
  clu <- c()
  for( k in unique(acount$Name)){
    df <- acount[grep(k, acount$Name),-1]
    if(nrow(df)>1){
      df <- colSums(df)
    }
    df$Name <- k
    df <- data.frame(df)
    clu <- rbind(clu,df)
  }

  acount <- clu[,c(ncol(clu), seq(1,ncol(clu)-1))]
  write.table(acount,file.path(dir,"rc", "readcount_isoacceptors.txt"), sep = "\t", row.names = F, quote = F)
  gtable(acount, container=RC_aa)
  insert(st,"Done : Counting of isoacceptors", do.newline = TRUE )
  insert(st, " ", do.newline = TRUE)
}
