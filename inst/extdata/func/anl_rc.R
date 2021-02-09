anl_rc <- function(ref_name){
  dir <- getwd()
  workDir <- file.path(dir, "post")
  bamlist <- list.files(workDir, pattern = "\\.bam$", recursive = T)
  bamlist <- basename(bamlist)

  # predic tRNADB ID
  cluinfo <-  read.table(system.file("extdata", file.path("refer",ref_name, "GtRNAdb_clusterInfo.txt"),
                                     package = "tReasure",mustWork = TRUE), sep = "\t", fill = T,header = T)
  cluinfo <- cluinfo[, c(1,2,4)]


  rc_df <- function(dfbn){
    out <- data.frame(do.call('rbind', strsplit(as.character(dfbn$seqnames), split = ":")))
    out$X1 <- gsub("^chr\\d{1,3}|^chr\\D{1}|^chr\\d{1,2}_\\D.*\\.", "", out$X1)
    out$X1 <- gsub(".tRNA", "trna", out$X1)
    out$X1 <- gsub("tRNA", "trna", out$X1)
    out$X1 <- paste(out$X3, out$X1, sep = ".")
    out <- data.frame(tRNAscan.SE.ID= out$X1)
    dfbn <- cbind(out, dfbn[,c(1:3)])
    dfbn <- dfbn[,-5]
    colnames(dfbn)[4] <- "Counts"
    return(dfbn)
  }

  dfbn1 <- idxstatsBam(file.path(workDir,bamlist[1]))
  dfbn <- rc_df(dfbn1)
  dfbn <- dfbn[,c(1,2)]

  for ( bn in bamlist){
    df <- idxstatsBam(file.path(workDir,bn))
    df <- rc_df(df)
    colnames(df)[4] <- gsub("\\_.*","",bn)
    dfbn <- left_join(dfbn, df)
  }

  write.table(dfbn, file.path(dir,"rc", "readcount_all_tRNAs.txt"), sep = "\t", row.names = F, quote = F)

  id <- data.frame(do.call('rbind', strsplit(as.character(dfbn$tRNAscan.SE.ID), split= "-")))
  dfbn$tRNAscan.SE.ID <- id$X1

  dfbn <- dfbn[dfbn$tRNAscan.SE.ID%in%cluinfo$tRNAscan.SE.ID,]
  dfbn <- left_join(cluinfo, dfbn)

  write.table(dfbn[,-c(1,3,4,5)], file.path(dir,"rc", "readcount_highconfidence_tRNAs.txt"), sep = "\t", row.names = F, quote = F)


  tcount <- dfbn[,-c(1,3,4,5)]
  colnames(tcount)[1] <- "Name"

  write.table(tcount,file.path(dir,"rc","readcount_trnas.txt"), sep = "\t", row.names = F, quote = F)

  gtable(tcount, container=RC_trna)
  insert(st, " ", do.newline = TRUE)
  insert(st,"Done : Counting of individual tRNAs", do.newline = TRUE )


  ## isodecoder
  sdfbn <- dfbn[,-c(1,2,4,5)]
  clu <- c()
  uni_clu <- unique(sdfbn$iso)
  for (i in uni_clu){
    dd <- sdfbn[grep(i, sdfbn$iso), ]
    if(nrow(dd) > 1){
      d <- colSums(dd[,-c(1)])
    }else{
      d <- data.frame(dd[,-c(1)])
    }
    clu <- rbind(clu, d)
  }
  ccount <- data.frame(iso = uni_clu, clu)
  colnames(ccount )[1] <- "Name"
  write.table(ccount ,file.path(dir,"rc", "readcount_isodecoders.txt"), sep = "\t", row.names = F, quote = F)
  gtable(ccount, container=RC_codon)
  insert(st, " ", do.newline = TRUE)
  insert(st,"Done : Counting of isodecoders", do.newline = TRUE )


  ## isoacceptor
  aa <- c()
  sdfbn$iso <- gsub("\\-\\d{1,2}.*", "", sdfbn$iso)
  uni_aa <- unique(sdfbn$iso)
  for (i in uni_aa){
    dd <- sdfbn[grep(i, sdfbn$iso), ]
    if(nrow(dd) > 1){
      d <- colSums(dd[,-1])
    }else{
      d <- data.frame(dd[,-1])
    }
    aa <- rbind(aa, d)
  }

  aa <- data.frame(Name= uni_aa, aa)
  acount <- aa
  write.table(acount,file.path(dir,"rc", "readcount_isoacceptors.txt"), sep = "\t", row.names = F, quote = F)
  gtable(acount, container=RC_aa)
  insert(st,"Done : Counting of isoacceptors", do.newline = TRUE )

  insert(st, " ", do.newline = TRUE)
}
