#' window
#'
#' @param ()
#' @return tReasure01
#' @export
tReasure01 <- function(){
  # install
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

  bpkg <- c("pillar","QuasR","DESeq2","edgeR", "Rsamtools","seqinr","ShortRead")

  ibpak <- function(bpkg){
    new.pkg <- bpkg[!(bpkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg))
      BiocManager::install(new.pkg)
    sapply(bpkg, require, character.only = TRUE)
  }

  ibpak(bpkg)
  pkg <- c("gWidgets2","cairoDevice","plotrix","tidyverse",
           "gridExtra","ggplot2","grid","dplyr","statmod",
           "future", "stringr")
  
  
  ipak <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg))
      install.packages(new.pkg, dependencies = TRUE, repos = "http://cran.r-project.org")
    sapply(pkg, require, character.only = TRUE)
  }
  
  ipak(pkg)
  
  # R_path <-.libPaths()[1]
  # if(!dir.exists(file.path(R_path,"gWidgets2RGtk2"))){
  #   download.file("http://treasure.pmrc.re.kr/data/gWidgets2RGtk2.zip", destfile = "gWidgets2RGtk2.zip")
  #   unzip(zipfile = "gWidgets2RGtk2.zip", exdir= R_path)
  # }
  # 
  # if(!dir.exists(file.path(R_path,"RGtk2"))){
  #   download.file("http://treasure.pmrc.re.kr/data/RGtk2.zip", destfile = "RGtk2.zip")
  #   unzip(zipfile = "RGtk2.zip", exdir= R_path)
  # }
  
  # libraray
  pkg <- c("gWidgets2","cairoDevice","plotrix","tidyverse", "RGtk2","gWidgets2RGtk2",
           "gridExtra","ggplot2","grid","dplyr","statmod","future", "stringr",
           "QuasR","DESeq2","edgeR", "Rsamtools","seqinr","ShortRead", "tReasure")
  sapply(pkg, require, character.only = TRUE)

  # Sys.setlocale('LC_ALL','C')
  #-------------------------------------------------------------------------------------
  # setting the PATH for TEST : system.file('', package="tReasure")
  #......................................................................................#
  intro <- system.file("extdata", "intro.png", package = "tReasure", mustWork = TRUE)
  cl_name <-  read.table(system.file("extdata", "class_name.txt", package = "tReasure",mustWork = TRUE), sep = "\t", fill = T,header = T, as.is = T)
  #-------------------------------------------------------------------------------------
  # function
  #......................................................................................#

  stopFuture <- function(x){
    tools::pskill(x,signal = tools::SIGTERM)
    tools::pskill(x,signal = tools::SIGKILL)
  }

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

  sel_sample <-  function(h,...){
    sample2 <- read.table(svalue(sel_button), header = T, sep = "\t")
    sampath <- gsub("sample.txt", "", svalue(sel_button))
    setwd(sampath)
    dir <- getwd()
    if(!dir.exists(paste(dir, "/pre", sep = ""))){
      dir.create(paste(dir, "/pre", sep = ""), recursive = TRUE)}
    if(!dir.exists(paste(dir, "/stat/plot", sep = ""))){
      dir.create(paste(dir, "/stat/plot", sep = ""), recursive = TRUE)}

    sample2$Batch <- as.character(sample2$Batch)
    tbl[] <- sample2
    insert(st, " ", do.newline = TRUE)
    insert(st, "Complete uploading the sample list.", do.newline = TRUE)
    if(length(grep("control", sample2$Group)) < 2 | length(grep("test", sample2$Group)) < 2){
      insert(st, "Warnning : There should be at least 2 samples per group for differentially expression analysis.")}
    # save update sample list
    addHandlerChanged(tbl, handler = function(h, ...){
      new_sample <- tbl[]
      SampleFile <- file.path(svalue(fq_dir), "sample.txt")
      write.table(new_sample, SampleFile, sep = "\t", quote = FALSE, row.names = FALSE)
      insert(st, " ", do.newline = TRUE)
      insert(st, ".", do.newline = TRUE)
    })
  }

  anl_trim <-  function(h,...){
    plan(multiprocess)
    if(file.exists("sample.txt") == FALSE & file.exists("*.fastq") == FALSE){
      gmessage("Warning : No file found matching sample list or fastq files")
    }else{
      insert(st, "", do.newline = TRUE)
      insert(st,"Start : Quality Control", do.newline = TRUE )
      dir <- getwd()

      sam_trim <- file.path(dir,"pre","sample.txt")
      sam_trim <- sub(".txt", "_trim.txt", sam_trim)

      sFile1 <- read.delim("sample.txt", header = TRUE, as.is = TRUE)

      sFileq <- sFile1[,c("FileName","SampleName")]
      sFileq$FileName <- file.path(dir, 'pre', sFileq$SampleName)
      sFileq$FileName <- gsub("$","_qc.fastq", sFileq$FileName)

      sFile2 <- sFile1[,c("FileName","SampleName")]
      sFile2$FileName <- file.path(dir, 'pre', sFile2$SampleName)
      sFile2$FileName <- gsub("$","_trim.fastq", sFile2$FileName)

      write.table(sFile2, sam_trim, sep = "\t", quote = FALSE, row.names = FALSE)

      # preprocessing future-------------------
      qnumber <- as.numeric(svalue(q))
      qf <- function(){
        for( i in sFile1$FileName){
          f <- FastqStreamer(i,readerBlockSize=1000)
          while(length(fq <- yield(f))){
            qPerBase = as(quality(fq), "matrix")
            qcount = rowSums( qPerBase <= qnumber )
            qcount[is.na(qcount)] = 0
            writeFastq(fq[qcount == 0],
                       file.path(dir, "pre", paste0(gsub(".fastq","_qc.fastq", basename(i)))), mode="a")}}
      }
      qff <- future(
        # qpid <- Sys.getpid()
        # save(qpid, file = file.path(dir,"qpid"))
        qf()
      )

      save(qff, file = "q.RData")

      qc_status <- function(){
        for(i in sFileq$SampleName){
          a <- sFile1$FileName[grep(i, sFile2$SampleName)]
          t <- sFileq$FileName[grep(i, sFile2$SampleName)]
          insert(st, paste("Screen : ", a), do.newline = TRUE)
          repeat{
            Sys.sleep(1)
            insert(st, ".", do.newline = FALSE)
            if(file.exists(t) == TRUE) break
          }
          insert(st, " ", do.newline = TRUE)
        }
      }

      qc_status()


      pre <- function(){
        apid <- Sys.getpid()
        save(apid, file = file.path(dir,"apid"))
        resa <- preprocessReads(filename = sFileq$FileName,outputFilename = sFile2$FileName,
                                minLength = 10, Rpattern = aseq)
        return(resa)}

      preno <- function(){resno <- preprocessReads(filename = sFileq$FileName, outputFilename = sFile2$FileName,
                                                   minLength = 10)
      return(resno)}

      if(svalue(adapt) == "Illumina smallRNA 3' adapter"){
        aseq <- "TGGAATTCTCGGGTGCCAAGG"
        #minLength = as.numeric(svalue(min))
        a <- future(pre())
      }else if(svalue(adapt) =="Illumina universal adapter"){
        aseq <- "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
        a <- future(pre())
      }else if(svalue(adapt) =="SOLiD Adapter"){
        aseq <- "CGCCTTGGCCGTACAGCAG"
        a <-future(pre())
      }else{
        a <-future(preno())
      }


      # preprocessing status-------------------

      preprocess_status <- function(){
        for(i in sFile2$SampleName){
          a <- sFileq$FileName[grep(i, sFile2$SampleName)]
          t <- sFile2$FileName[grep(i, sFile2$SampleName)]
          insert(st, paste("Filtering : ", a), do.newline = TRUE)
          repeat{
            Sys.sleep(1)
            insert(st, ".", do.newline = FALSE)
            if(file.exists(t) == TRUE) break
          }
          insert(st, " ", do.newline = TRUE)
        }
      }

      #status <- future(preprocess_status())


      preprocess_status()
      res <- value(a)

      # display
      rest <- data.frame(Names = row.names(res),res)
      colnames(rest) <- gsub(".fastq","", colnames(rest))
      write.table(rest, "./pre/trim_res.txt", sep="\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
      rest <- rest[-c(2,5,6),]
      res_trim <- gtable(data.frame(rest), container = ggr21)
      #
      file.remove(sFileq$FileName)
      #
      insert(st, "", do.newline = TRUE)
      insert(st,"Done : Quality Control. Click! Next tab of Alignment & Counting.", do.newline = TRUE )
      insert(st, ".", do.newline = TRUE)
    }
  }

  anl_align <- function(h,...){
    insert(st,"Start : Alignment & Read Counting", do.newline = TRUE )
    dir <- svalue(fq_dir)
    insert(st, paste("Check genome.... "), do.newline = TRUE)
    # download_refer
    dw_refer <- function(url, sub, sub.zip){
      R_path <- system.file( "extdata", package = "tReasure", mustWork = TRUE)
      if(!dir.exists(file.path(R_path,"refer",sub))){
        dir.create(file.path(R_path,"refer",sub), recursive = TRUE)
      }
      Rsub_path <- file.path(R_path,"refer",sub)
      if(!file.exists(file.path(Rsub_path, paste0(sub,"_artificial.fa")))){
        insert(st, "Download and unzip genome indexes. It takes several minutes. Please wait....", do.newline = TRUE)
        download.file(url, destfile = file.path(Rsub_path, sub.zip))
        unzip(zipfile = file.path(Rsub_path, sub.zip), exdir= Rsub_path)
      }
    }

    ref_name <- cl_name$fa[grep(svalue(gsel_button),cl_name$P4)]

    url <- file.path("http://treasure.pmrc.re.kr/data/genome", svalue(ref_P1), ref_name, paste0(ref_name, ".zip"))
    dw_refer(url, ref_name, paste0(ref_name, ".zip"))
    genome = system.file( "extdata", file.path("refer",ref_name, paste0(ref_name, "_artificial.fa")), package = "tReasure", mustWork = TRUE)

    sFile <- file.path(dir,"pre", "sample_trim.txt")
    sFile2 <- read.delim(sFile)

    # alignment future----------------------

    bowtie <- function(sFile, genome){
      plan(multiprocess)
      bpid <- Sys.getpid()
      save(bpid, file = file.path(dir,"bpid"))
      proj <- qAlign(sFile, genome = genome, aligner = "Rbowtie",alignmentParameter = "-v 3 --best", clObj = NULL)
      return(proj)
    }

    bow <- future({bowtie(sFile, genome)})

    # alignment status----------------------

    algin_status_pre <- function(){
      for(i in sFile2$SampleName){
        a <- sFile2$FileName[grep(i, sFile2$SampleName)]
        Sys.sleep(10)
        insert(st, paste("pre-mapping   : ", a), do.newline = TRUE)
        Sys.sleep(10)
        repeat{
          Sys.sleep(1)
          insert(st, ".", do.newline = FALSE)
          dir <- getwd()
          list <- list.files(file.path(dir,"pre"), pattern=".bam$", full.names = TRUE)
          if(length(grep(i, list))>0)break
        }
        insert(st, " ", do.newline = TRUE)
      }
    }

    algin_status_pre()

    # pre-mapping
    projt_pre <- value(bow)
    alinstat <- data.frame(Names=row.names(alignmentStats(projt_pre)),alignmentStats(projt_pre))
    write.table(alinstat, "./pre/preproccessing_align_stat.txt",sep = "\t")

    # QC
    save(projt_pre, file = "./pre/premappingQC.RData")
    #qQCReport(projt_pre, pdfFilename = "./pre/qc_report.pdf", useSampleNames = TRUE)

    # remove reads(premature-tRNAs, non_tRNAs)
    insert(st, ".", do.newline = TRUE)
    insert(st,"Removing reads aligned premauture tRNA and no tRNAs", do.newline = TRUE )
    anl_rm(ref_name)


    # post-mapping
    insert(st, ".", do.newline = TRUE)
    clu_genome = system.file( "extdata", file.path("refer", ref_name, paste0(ref_name, ".tRNAscan_mature.fa")), package = "tReasure", mustWork = TRUE)

    matfq_list <- list.files(file.path(dir, "post"), pattern = "_trim_mature.fastq$", full=TRUE)
    matsFile2 <- data.frame(FileName = matfq_list, SampleName = gsub("_trim_mature.fastq","",basename(matfq_list)))
    write.table(matsFile2, file.path(dir, "post", "sample_mature.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
    matsFile <- file.path(dir,"post", "sample_mature.txt")

    bow.post <- future({bowtie(matsFile,clu_genome)})

    # alignment status----------------------

    algin_status_post <- function(){
      for(i in matsFile2$SampleName){
        a <- matsFile2$FileName[grep(i, matsFile2$SampleName)]
        Sys.sleep(1)
        insert(st, paste("post-mapping   : ", a), do.newline = TRUE)
        repeat{
          Sys.sleep(1)
          insert(st, ".", do.newline = FALSE)
          dir <- getwd()
          list <- list.files(file.path(dir,"post"), pattern=".bam$", full.names = TRUE)
          if(length(grep(i, list))>0)break
        }
        insert(st, " ", do.newline = TRUE)
      }
    }

    algin_status_post()

    # premmapping
    projt_post <- value(bow.post)


    alinstat <- data.frame(Names=row.names(alignmentStats(projt_post)),alignmentStats(projt_post))
    write.table(alinstat, "./post/postprocessing_align_stat.txt",sep = "\t")

    # Readcounting
    insert(st,"counting", do.newline = TRUE )
    anl_rc(ref_name)
    insert(st,"Done : Alignment & Read Counting. Click! Next tab of Filtering.", do.newline = TRUE )
    insert(st, ".", do.newline = TRUE)

  }

  anl_rc <- function(ref_name){
    dir <- getwd()
    workDir <- file.path(dir, "post")
    bamlist <- list.files(workDir, pattern = "\\.bam$", recursive = T)
    bamlist <- basename(bamlist)

    rc_df <- function(dfbn){
      dfbn$seqnames <-gsub("_t",":t", dfbn$seqnames) 
      dfbn$seqnames <-gsub("_N",":N", dfbn$seqnames) 
      dfbn$seqnames <-gsub("_none",":none", dfbn$seqnames) 
      dfbn$seqnames <-gsub("_high",":high", dfbn$seqnames)
      chr <- data.frame(stringr::str_split_fixed(dfbn$seqnames, "[.]", 2))
      out <- data.frame(stringr::str_split_fixed(chr$X2, ":", 4))
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


  anl_rm <- function(ref_name){
    dir <- getwd()
    genomeDir <- system.file("extdata", "refer", ref_name , package = "tReasure", mustWork = TRUE)
    workDir <- file.path(dir, "pre")

    ## bed file
    bed <- read.table(file.path(genomeDir, paste0(ref_name, ".tRNAscan_pre-tRNAs.bed12")), sep = "\t")
    len <- data.frame(refname = bed$V4, reflen= bed$V11)
    len %>% mutate_if(is.factor, as.character) -> len

    for ( i in 1:nrow(len)){
      if(length(grep("[,]",len$reflen[i])) != 0){
        d <- data.frame(do.call('rbind',strsplit(as.character(len$reflen[i]), split=",")), stringsAsFactors = F)
        d <- lapply(d,as.numeric)
        len$reflen[i]<- as.character(d[[1]] + d[[2]])
      }
    }
    len$reflen <- as.numeric(len$reflen)

    # bam file

    bamlist <- list.files(workDir, pattern = "\\.bam$", recursive = T)
    bamlist <- basename(bamlist)

    # function ---------
    mfun <- function(x){
      tid <- gsub("\\s.*","" ,data.frame(ShortRead::id(x))[,1])
      x[tid%in%tb1$qname]
    }
    matcher <- function(pattern, x) {
      ind = gregexpr(pattern, x)[[1]]
      start = as.numeric(ind)
      end = start + attr(ind, "match.length")- 2
      apply(cbind(start,end), 1, function(y) substr(x, start=y[1], stop=y[2]));
    }
    doone <- function(c, cigar) {
      pat <- paste("\\d+", c , sep="")
      sum(as.numeric(matcher(pat, cigar)), na.rm=T)
    }
    ## takes a cigar string and parses it, not very fast but...
    cigarsums <- function(cigar, chars=c("D","I")) {
      sapply (chars, doone, cigar)
    }
    fun <- function(x){
      tid <- gsub("\\s.*","" ,data.frame(ShortRead::id(x))[,1])
      x[tid%in%df$qname]
    }

    # screen bamlist ---------
    for( bn in bamlist){
      b1 <- scanBam(file.path(workDir,bn))
      b1 <- data.frame(b1[[1]], stringsAsFactors = F)
      b1$rname <- as.character(b1$rname)
      tb1 <- b1[grep(".tRNA", b1$rname),]
      gb1 <- b1[!grepl(".tRNA", b1$rname),]

      tb1 <- filter(tb1, tb1$pos > 50 )
      tb1$refname <- gsub("\\::.*", "", tb1$rname)
      sname <- sub("\\_.*", "", bn)
      fq <- file.path(workDir, paste0(sname, "_trim", ".fastq"))

      # cigar ------------
      con <- unique(c(grep("D",tb1$cigar),grep("I", tb1$cigar)))
      if(length(con) == 0 ){
        tb1$end <- tb1$pos + tb1$qwidth -1
      }else{
        tb1t <- tb1[con,]
        tb1o <- tb1[-con,]
        tb1tt <- sapply(tb1t$cigar, cigarsums)
        tb1tt <- data.frame(t(tb1tt))
        tb1t$end <- tb1t$pos + tb1t$qwidth -1 + tb1tt$D - tb1tt$I
        tb1o$end <- tb1o$pos + tb1o$qwidth -1
        tb1 <- rbind(tb1t, tb1o)
      }
      df <- left_join(tb1, len)
      df <- filter(df, end <= (df$reflen-44))

      filterFastq(fq, destinations = file.path(dir, "post", paste0(sname,  "_trim_mature.fastq")), filter =fun , compress= FALSE)
      insert(st, paste("Completed : ", paste0(sub("\\_.*", "", bn),  "_trim_mature.fastq")), do.newline = TRUE)
    }
  }

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
      rownames(count) <- count$Name
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

  edger<- function(count){
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

    if(svalue(method_edgeR) == "Exact test"){
      test <- exactTest(y)
    }else if(svalue(method_edgeR) == "Quasi-likelihood F-test"){
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

  anl_deseq <- function(h,...){
    insert(st,"Statistical analysis  : DESeq2", do.newline = TRUE )
    dir <- getwd()
    if(!dir.exists(paste(dir, "/stat/plot", sep = ""))){
      dir.create(paste(dir, "/stat/plot", sep = ""), recursive = TRUE)}

    tcount <- read.delim(file.path(dir, "rc","filtered_readcount_trnas.txt"))
    ccount <- read.delim(file.path(dir, "rc","filtered_readcount_isodecoders.txt"))
    acount <- read.delim(file.path(dir, "rc","filtered_readcount_isoacceptors.txt"))

    t_out <- deseq(tcount)
    c_out <- deseq(ccount)
    a_out <- deseq(acount)

    stat_trna[] <- t_out
    stat_codon[] <- c_out
    stat_aa[] <- a_out

    write.table(t_out, "./stat/stat_trna_list.txt", sep="\t", quote = FALSE, row.names = F)
    write.table(c_out, "./stat/stat_isodecoder_list.txt", sep="\t", quote = FALSE, row.names = F)
    write.table(a_out, "./stat/stat_isoacceptor_list.txt", sep="\t", quote = FALSE, row.names = F)
    insert(st, " ", do.newline = TRUE)
    insert(st,"Done : Statistical analysis. Set the threshold value.", do.newline = TRUE )
    insert(st, ".", do.newline = TRUE)
  }

  anl_edger <- function(h,...){
    insert(st,"Statistical analysis  : EdgeR", do.newline = TRUE )
    insert(st,"It takes a few minutes. Please wait....", do.newline = TRUE )
    dir <- getwd()

    tcount <- read.delim(file.path(dir, "rc", "filtered_readcount_trnas.txt"))
    ccount <- read.delim(file.path(dir, "rc", "filtered_readcount_isodecoders.txt"))
    acount <- read.delim(file.path(dir, "rc", "filtered_readcount_isoacceptors.txt"))

    t_out <- edger(tcount)
    c_out <- edger(ccount)
    a_out <- edger(acount)

    stat_trna[] <- t_out
    stat_codon[] <- c_out
    stat_aa[] <- a_out

    write.table(t_out, "./stat/stat_trna_list.txt", sep="\t", quote = FALSE, row.names = F)
    write.table(c_out, "./stat/stat_isodecoder_list.txt", sep="\t", quote = FALSE, row.names = F)
    write.table(a_out, "./stat/stat_isoacceptor_list.txt", sep="\t", quote = FALSE, row.names = F)
    insert(st, " ", do.newline = TRUE)
    insert(st,"Done : Statistical analysis. Set the threshold value.", do.newline = TRUE )
    insert(st, ".", do.newline = TRUE)
  }

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

  volcano <- function(width, height, res){
    out <- read.delim("./stat/stat_trna_list.txt")

    if(svalue(widget_list$fdr_s) == "Bejamini-Hochberg"){
      outt <- select(out, Names, logFC, Benjamini)
      colnames(outt)[ncol(outt)]<-"FDR"
    }else if(svalue(widget_list$fdr_s) == "Bonferroni"){
      outt <- select(out, Names, logFC, Bonferroni)
      colnames(outt)[ncol(outt)]<-"FDR"
    }else{
      outt <- select(out, Names, logFC, FDR)
    }

    pval <- svalue(widget_list$pval)
    fc <- svalue(widget_list$FC)


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

  barplot <- function(width, height, res){
    trna <- read.delim("./rc/readcount_isodecoders.txt")
    out <- read.delim("./stat/stat_isodecoder_list.txt")

    pval <- svalue(widget_list$pval)
    fc <- svalue(widget_list$FC)

    if(svalue(widget_list$fdr_s) == "Benjamini-Hochberg"){
      dw <- filter(out, logFC < -fc, out$Benjamini < pval)
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


  #-------------------------------------------------------------------------------------
  #  tReasure
  #......................................................................................#
  window <- gwindow("tReasure")
  addHandlerUnrealize(window, handler = function(h,...) {
    val <- gconfirm("Are you sure you want to exit tReausre?", parent=h$obj)
    if(as.logical(val)){
      if(dir.exists("apid")){
        load("apid")
        stopFuture(apid)}
      if(dir.exists("bpid")){
        load("bpid")
        stopFuture(bpid)}
      dispose(window)
      gtkMainQuit()
    }else{
      return(TRUE)
    }
  })



  mother <- ggroup(container = window, horizontal = FALSE)
  size(mother) <- c(1000,740)

  # Sub window (Notebook)------------
  notebook <- gnotebook(container = mother)
  gr0 <- ggroup(container = notebook, label ="  Introduction  ")
  gr1 <- ggroup(container = notebook, label ="  Uploading Samples  ", horizontal = FALSE)
  gr2 <- ggroup(container = notebook, label ="  Quality Control  ",horizontal = FALSE)
  gr3 <- ggroup(container = notebook, label ="  Alignment & Counting  ",horizontal = FALSE)
  gr4 <- ggroup(container = notebook, label ="  Filtering   ",horizontal = FALSE)
  #gr5 <- ggroup(container = notebook, label ="  Exploration   ",horizontal = FALSE)
  gr6 <- ggroup(container = notebook, label ="  DEtRNA Detection ",horizontal = FALSE)
  gr7 <- ggroup(container = notebook, label ="  Visualization ",horizontal = FALSE)
  child <- gframe("  Progress  ", container = mother, horizontal = FALSE)
  size(child) <- c(1000,100)
  st <- gtext("  ",container = child)
  size(st) <- c(1000,100)

  #-------------------------------------------------------------------------------------
  #  gr0. Introduction
  #......................................................................................
  ggr1 <- ggroup(container = gr0, horizontal = TRUE, fill = "both", expand = TRUE)
  gimage(intro, container = ggr1)
  #-------------------------------------------------------------------------------------
  #  gr1. Uploading Samples
  #......................................................................................
  addSpace(gr1, 20)
  ggr1 <- ggroup(container = gr1,  spacing = 10) ; size(ggr1) <- c(896,532)
  addSpace(ggr1, 10)
  paned.1 <- gpanedgroup(container = ggr1, horizontal = FALSE, spacing = 10)
  tmp.1 <- gframe(" Create the sample list ", container = paned.1, horizontal = FALSE, spacing = 10 ) ; size(tmp.1) <- c(250,284)
  addSpace(tmp.1, 10)
  glabel(" Select a directory of FASTQ files : ", container = tmp.1, anchor = c(-1,0))
  addSpace(tmp.1, 10)
  fq_dir <- gfilebrowse(text = " ", quote = FALSE, type = "selectdir", container = tmp.1)
  addSpace(tmp.1, 20)
  make_button <- gbutton("RUN", container = tmp.1, handler = mk_sample)

  tmp.11 <- gframe(" [*OPTION] Upload the sample list ", container = paned.1, horizontal = FALSE, spacing = 10) ; size(tmp.11) <- c(250,200)
  addSpace(tmp.11, 10)
  glabel(" Select a file of \n the pre-made sample list : ", container = tmp.11, anchor = c(-1,0))
  addSpace(tmp.11, 20)

  sel_button <- gfilebrowse(text=" ", quote = FALSE, type="open", container = tmp.11,
                            filter=list("*.txt" = list(patterns = c("*.txt")), "*.*" = list(patterns = c("*"))),
                            handler = sel_sample)

  tbl <- gdf(data.frame(FileName=I(list("Path/samplename.fastq", NA, NA)),
                        SampleName=I(list("samplename", NA, NA)),
                        Group=I(list("control", NA, NA)),
                        Batch=I(list("batch I", NA, NA))), container = ggr1, multiple = FALSE)


  addHandlerChanged(fq_dir, handler = function(h,...){
    setwd(svalue(fq_dir))
    dir <- getwd()
    if(!dir.exists(paste(dir, "/pre", sep = ""))){
      dir.create(paste(dir, "/pre", sep = ""), recursive = TRUE)}

    if(!dir.exists(paste(dir, "/post", sep = ""))){
      dir.create(paste(dir, "/post", sep = ""), recursive = TRUE)}

    if(!dir.exists(paste(dir, "/stat/plot", sep = ""))){
      dir.create(paste(dir, "/stat/plot", sep = ""), recursive = TRUE)}
    if(!dir.exists(paste(dir, "/rc", sep = ""))){
      dir.create(paste(dir, "/rc", sep = ""), recursive = TRUE)}
  })



  #-------------------------------------------------------------------------------------
  #  gr2. Quality Control
  #......................................................................................#
  addSpace(gr2, 20)
  ggr21 <- ggroup(container = gr2,  spacing = 10); size(ggr21) <- c(896,532)
  addSpace(ggr21, 10)
  paned.2 <- gpanedgroup(container = ggr21, horizontal = FALSE, spacing = 10)
  tmp.2 <- gframe("  Settings  ", container = paned.2, horizontal = FALSE, spacing = 10); size(tmp.2) <- c(250,532)
  addSpace(tmp.2, 10)

  glabel(" Information of adapter : ", container = tmp.2, anchor = c(-1,0))
  addSpace(tmp.2, 10)

  adapt <- gradio(c("Illumina smallRNA 3' adapter", "Illumina universal adapter","SOLiD adapter","No adapter"), container = tmp.2)
  addSpace(tmp.2, 10)

  glabel(" Minimum quality : ", container = tmp.2, anchor = c(-1,0))
  addSpace(tmp.2, 10)
  q <- gradio(c("25", "30"), container = tmp.2)
  addSpace(tmp.2, 10)

  glabel(" Minimum lenght : ", container = tmp.2, anchor = c(-1,0))
  addSpace(tmp.2, 10)
  min <- gradio(c("10"), container = tmp.2)
  addSpace(tmp.2, 20)

  trim_button <- gbutton("RUN", container = tmp.2, expand = FALSE, handler = anl_trim)


  #-------------------------------------------------------------------------------------
  #  gr3. Alignment & Counting
  #......................................................................................
  addSpace(gr3, 20)
  ggr31 <- ggroup(container = gr3,  spacing = 10); size(ggr31) <- c(896,532)
  addSpace(ggr31, 10)
  paned.3 <- gpanedgroup(container = ggr31, horizontal = FALSE, spacing = 10)
  tmp.3 <- gframe("  Settings  ", container = paned.3 , horizontal = FALSE, spacing = 10); size(tmp.3) <- c(250,260)
  pangr <- ggroup(container = paned.3, horizontal =FALSE)
  addSpace(tmp.3, 10)

  RC_view <- gnotebook(container = ggr31 ); size(RC_view) <- c(700, 500)
  RC_trna <- ggroup(container = RC_view, horizontal = FALSE, label = " tRNA level ")
  RC_codon <- ggroup(container = RC_view, horizontal = FALSE, label = " Isodecoder level")
  RC_aa <- ggroup(container = RC_view, horizontal = FALSE, label = "Isoacceptor level")

  glabel(" Genome assembly : ", container = tmp.3, anchor = c(-1,0))
  addSpace(tmp.3, 10)


  ### child_window--------------------
  child_window <- gwindow("Search genome assembly", parent = window, width = 600, height = 300, visible = FALSE )
  ggroup(container = child_window, spacing = 10)
  gr_child <- ggroup(container = child_window, spacing = 10)

  addSpace(gr_child, 20)
  gr_frame <- gframe(" ", container = gr_child, horizontal = FALSE, spacing = 10)
  size(gr_frame) <- c(550, 250)
  addSpace(gr_child, 20)

  addSpace(gr_frame, 10)
  c1 <- glabel("  Domain : ", container = gr_frame, anchor = c(-1,0))
  ref_P1 <- gcombobox(c(" ", "Eukaryota" ,"Archaea", "Bacteria"),container = gr_frame) # container = tmp.3
  addSpace(gr_frame, 10)

  c2 <- glabel("  Next 1 : ", container = gr_frame, anchor = c(-1,0))
  ref_P2 <- (v1names <- gcombobox("  ",container = gr_frame)) # container = tmp.3
  addSpace(gr_frame, 10)

  c3 <- glabel("  Next 2 : ", container = gr_frame, anchor = c(-1,0))
  ref_P3 <- (v2names <- gcombobox("  ",container = gr_frame)) # container = tmp.3
  addSpace(gr_frame, 10)

  c4 <- glabel("  Species : ", container = gr_frame, anchor = c(-1,0))
  ref_P4 <- (v3names <- gcombobox(" ",container = gr_frame)) # container = tmp.3

  p2Nms <- function(d, envir=.GlobalEnv)  unique(cl_name$P2[grep(svalue(ref_P1), cl_name$P1)], envir=envir)
  p3Nms <- function(d, envir=.GlobalEnv)  unique(cl_name$P3[grep(svalue(ref_P2), cl_name$P2)], envir=envir)
  p4Nms <- function(d, envir=.GlobalEnv)  unique(cl_name$P4[grep(svalue(ref_P3), cl_name$P3)], envir=envir)

  addHandlerChanged(ref_P1, handler = function(h,...){
    ref_P2[] <- p2Nms(svalue(ref_P2))
  })

  addHandlerChanged(ref_P2, handler = function(h,...){
    ref_P3[] <- p3Nms(svalue(ref_P3))
  })

  addHandlerChanged(ref_P3, handler = function(h,...){
    ref_P4[] <- p4Nms(svalue(ref_P4))
  })

  addSpace(gr_frame, 20)
  g_run <- gbutton("OK", container = gr_frame, handler = function(h, ...){
    svalue(gsel_button) <<- svalue(ref_P4)
    visible(child_window) <- FALSE
    insert(st, "Click the RUN", do.newline = TRUE)
  })
  addSpace(gr_child, 20)

  gsel_button <- gbutton("Search genome assembly", container = tmp.3, handler = function(h, ...){
    visible(child_window) <- TRUE
  })
  addSpace(tmp.3, 20)

  anl_button <- gbutton("RUN", container=tmp.3, handler = anl_align)

  addSpace(pangr, 10)
  # tmp.31 <- gframe(" [*OPTION] Working directory ", container = pangr , horizontal = FALSE, spacing = 10); size(tmp.31) <- c(250,70)
  # addSpace(tmp.31, 10)
  #
  # glabel(" Select the directory of the FASTQ files : ", container = tmp.31, anchor = c(-1,0))
  # addSpace(tmp.31, 10)

  tmp.32 <- gframe(" [*OPTION] Load the readcounts files ", container = pangr , horizontal = FALSE, spacing = 10); size(tmp.32) <- c(250,120)
  addSpace(tmp.32, 10)

  glabel(" Select the file of readcounts : ", container = tmp.32, anchor = c(-1,0))
  addSpace(tmp.32, 10)

  # selfq_button <- gfilebrowse(text = "", quote = FALSE, type = "selectdir", container = tmp.31,
  #                             filter=list("*.fa" = list(patterns = c("*.fa")), "*.*" = list(patterns = c("*"))))

  seltrn_button <- gfilebrowse(text = "tRNA level", quote = FALSE, type = "open", container = tmp.32,
                               filter=list("*.txt" = list(patterns = c("*.txt")), "*.*" = list(patterns = c("*"))))

  selco_button <- gfilebrowse(text = "Isodecoder level", quote = FALSE, type = "open", container = tmp.32,
                              filter=list("*.txt" = list(patterns = c("*.txt")), "*.*" = list(patterns = c("*"))))

  selaa_button <- gfilebrowse(text = "Isoacceptor level", quote = FALSE, type = "open", container = tmp.32,
                              filter=list("*.txt" = list(patterns = c("*.txt")), "*.*" = list(patterns = c("*"))))


  addHandlerChanged(seltrn_button, handler = function(h,...){
    tcount <- read.delim(svalue(seltrn_button))
    gtable(tcount, container=RC_trna)
  })

  addHandlerChanged(selco_button, handler = function(h, ...){
    ccount = read.delim(svalue(selco_button))
    gtable(ccount, container=RC_codon)
  })

  addHandlerChanged(selaa_button,handler = function(h, ...){
    acount = read.delim(svalue(selaa_button))
    gtable(acount, container=RC_aa)
  })


  #-------------------------------------------------------------------------------------
  #  gr4. Filtering
  #......................................................................................
  addSpace(gr4, 20)
  ggr41 <- ggroup(container = gr4,  spacing = 10); size(ggr41) <- c(896,532)
  addSpace(ggr41, 10)
  paned.4 <- gpanedgroup(container = ggr41, horizontal = FALSE, spacing = 10)
  tmp.4 <- gframe("  Settings   ", container = paned.4, horizontal = FALSE, spacing = 10)
  size(tmp.4) <- c(250,284)
  addSpace(tmp.4, 20)

  prelyt <- gformlayout(container = tmp.4, spacing = 1.5)
  cfilv <- gcombobox(seq(0,10,by=1),label = "cpm value >=  ", container = prelyt)
  sfil <- gcombobox(seq(0,100, by=10), label = "sample (%) >=  ", container = prelyt)

  addSpace(tmp.4, 20)
  Pre_button <- gbutton("RUN", container = tmp.4, handler = anl_filter)

  tmp.42 <- gframe(" [*OPTION] Working directory ", container = paned.4, horizontal = FALSE, spacing = 10)
  size(tmp.42) <- c(250,200)
  addSpace(tmp.42, 10)
  glabel(" Select the directory of the readcounts files :", container = tmp.42, anchor = c(-1,0))
  addSpace(tmp.42, 10)

  selrc_button <- gfilebrowse(text = "", quote = FALSE, type = "selectdir", container = tmp.42,
                              filter=list("*.txt" = list(patterns = c("*.txt")), "*.*" = list(patterns = c("*"))))

  f_kep <- gtable(data.frame(No.reads=c("total", "keep", "remove")),container = ggr41, multiple = FALSE)

  #-------------------------------------------------------------------------------------
  #  gr5. Exploration
  #......................................................................................
  # addSpace(gr5, 20)
  # ggr51 <- ggroup(container = gr5, spacing = 10, fill=TRUE); #size(ggr51) <- c(996,550)
  # addSpace(ggr51, 10)
  # mds_view <- gnotebook(container = ggr51, spacing = 10, fill=TRUE, expand=TRUE); #size(mds_view) <- c(950, 550)
  # mds_plot <- ggroup(container = mds_view, horizontal = FALSE, label = " MDS plot ")
  # addSpace(ggr51, 10)
  # addSpace(gr5, 20)


  #-------------------------------------------------------------------------------------
  #  gr6. DEtRNA Detection
  #......................................................................................
  addSpace(gr6, 20)
  ggr61 <- ggroup(container = gr6,  spacing = 10); size(ggr61) <- c(896,532)
  addSpace(ggr61, 10)
  paned.6 <- gpanedgroup(container = ggr61, horizontal = FALSE, spacing = 10)
  addSpace(ggr61, 10)

  ggr62 <- gnotebook(container=paned.6, tab.pos  = 3) ; size(ggr62) <- c(250, 200)
  stat2 <- ggroup(container = ggr62, label = " DEseq2 ")
  stat1 <- ggroup(container = ggr62, label = " EdgeR ")

  tmp.61 <- gframe("  Statics  ", container = stat1, horizontal = FALSE, spacing = 10); size(tmp.61) <- c(235,150)
  addSpace(tmp.61, 10)
  glabel(" Stat.Method : ", container = tmp.61, anchor = c(-1,0))
  addSpace(tmp.61, 10)
  method_edgeR <- gradio(c("Exact test", "Likelihood F-test", "Quasi-likelihood F-test"), container = tmp.61)
  addSpace(tmp.61, 10)

  tmp.62 <- gframe("  Statics  ", container = stat2, horizontal = FALSE, spacing = 10); size(tmp.62) <- c(235,150)
  addSpace(tmp.62, 10)
  glabel(" Stat.Method : ", container = tmp.62, anchor = c(-1,0))
  addSpace(tmp.62, 10)
  method_deseq <- gradio("Walt test", container = tmp.62)
  addSpace(tmp.62, 10)


  statedgeR_button <- gbutton ("RUN EdgeR", container = tmp.61, handler = anl_edger)
  statdeseq_button <- gbutton("RUN DESeq2", container = tmp.62, handler = anl_deseq)


  widget_list <- list()
  ggr.63 <- gframe("  Threshold  ", container = paned.6, horizontal = FALSE, spacing = 10)
  addSpace(ggr.63, 10)
  fdr_label<- glabel(" Adjust.Method : ", container = ggr.63, anchor = c(-1,0))
  addSpace(ggr.63, 10)
  widget_list$fdr_s <- gcombobox(c("FDR", "Bonferroni","Benjamini-Hochberg"  ), container = ggr.63)
  addSpace(ggr.63, 10)

  statlyt <- gformlayout(container = ggr.63, spacing = 1.5)
  widget_list$pval <- gcombobox(c(0,0.001, 0.05, 0.01),label = "Adj P.value  < : ", container = statlyt)
  widget_list$FC <- gcombobox(c(0,1,1.5,2),label = "Fold change  > : ", container = statlyt)

  addSpace(ggr.63, 20)

  DE_view <- gnotebook(container = ggr61, expand=TRUE ); #size(DE_view) <- c(700, 532)
  DE_trna <- ggroup(container = DE_view, horizontal = FALSE, label=" tRNA level ")
  DE_codon <- ggroup(container = DE_view, horizontal = FALSE, label=" Isodecoder level")
  DE_aa <- ggroup(container=DE_view, horizontal = FALSE, label= "Isoacceptor level")

  stat_trna <- gtable(data.frame(), container = DE_trna)
  stat_codon <- gtable(data.frame(), container = DE_codon)
  stat_aa <- gtable(data.frame(), container = DE_aa)


  # Threshold handler  --------------------
  update_trna_frame <- function(...){
    data <- read.delim("./stat/stat_trna_list.txt")
    if(svalue(widget_list$fdr_s) == "Benjamini-Hochberg"){
      stat_trna[] <- filter(data, data$Benjamini < svalue(widget_list$pval) & abs(data$logFC) > svalue(widget_list$FC))
    }else if(svalue(widget_list$fdr_s) =="Bonferroni"){
      stat_trna[] <- filter(data, data$Bonferroni < svalue(widget_list$pval) & abs(data$logFC) > svalue(widget_list$FC))
    }else{stat_trna[] <- filter(data, data$FDR < svalue(widget_list$pval) & abs(data$logFC) > svalue(widget_list$FC))}
  }

  update_codon_frame <- function(...){
    data <- read.delim("./stat/stat_isodecoder_list.txt")
    if(svalue(widget_list$fdr_s) == "Benjamini-Hochberg"){
      stat_codon[] <- filter(data, data$Benjamini < svalue(widget_list$pval) & abs(data$logFC) > svalue(widget_list$FC))
    }else if(svalue(widget_list$fdr_s) =="Bonferroni"){
      stat_codon[] <- filter(data, data$Bonferroni < svalue(widget_list$pval) & abs(data$logFC) > svalue(widget_list$FC))
    }else{stat_codon[] <- filter(data, data$FDR < svalue(widget_list$pval) & abs(data$logFC) > svalue(widget_list$FC))}
  }

  update_aa_frame <- function(...){
    data <- read.delim("./stat/stat_isoacceptor_list.txt")
    if(svalue(widget_list$fdr_s) == "Benjamini-Hochberg"){
      stat_aa[] <- filter(data, data$Benjamini < svalue(widget_list$pval) & abs(data$logFC) > svalue(widget_list$FC))
    }else if(svalue(widget_list$fdr_s) =="Bonferroni"){
      stat_aa[] <- filter(data, data$Bonferroni < svalue(widget_list$pval) & abs(data$logFC) > svalue(widget_list$FC))
    }else{stat_aa[] <- filter(data, data$FDR < svalue(widget_list$pval) & abs(data$logFC) > svalue(widget_list$FC))}
  }

  callback <- function(h, ...){
    update_trna_frame()
    update_codon_frame()
    update_aa_frame()
  }

  sapply(widget_list, addHandlerChanged, handler = callback)
  DE_button <- gbutton(" Plot ", container = ggr.63)

  # Plot save  handler --------------------
  addHandlerClicked(DE_button, handler = function(h, ...){
    DEtRNA <- stat_trna[]
    DEcodon <- stat_codon[]
    DEaa <- stat_aa[]

    write.table(DEtRNA, "./stat/DEtrna_list.txt", sep = "\t", quote = FALSE, row.names = F)
    write.table(DEcodon, "./stat/DEisodecoder_list.txt", sep = "\t", quote = FALSE, row.names = F)
    write.table(DEaa, "./stat/DEisoacceptor_list.txt", sep = "\t", quote = FALSE, row.names = F)

    volcano(650, 460, 80)
    barplot(850, 460,80)
    pyramid(650, 460,80)

    insert(st, " ", do.newline = TRUE)
    insert(st,"Done : The plots of DEtRNA Detection. Click! Next tab of Visualization.", do.newline = TRUE )
    insert(st, ".", do.newline = TRUE)
  })


  #-------------------------------------------------------------------------------------
  #  gr7. Visualiztion
  #......................................................................................
  addSpace(gr7, 20)
  ggr71 <- ggroup(container = gr7,  spacing = 10); size(ggr71) <- c(996,550)
  addSpace(ggr71, 10)
  Plot_view <- gnotebook(container = ggr71,spacing = 10, expand=TRUE ); #size(Plot_view) <- c(950, 550)
  devs <- lapply(c("MDS plot","Plot_trnas", "Plot_isodecoders", "Plot_isoacceptors"), function(i) ggraphics(container = Plot_view, visible = TRUE, label = as.character(i)))

  # Plot viewe handler--------------------
  addSpace(ggr71, 10)
  addHandlerChanged(Plot_view, handler = function(h,...){
    gg <- h$obj[h$page.no]
    visible(gg) <- TRUE
    #print(p())
    if(h$page.no == "1"){
      load("./stat/plot/MDS_plot.RData")
      print(p())
    }else if(h$page.no == "2" ){
      load("./stat/plot/Volcano_Plot.RData")
      print(p())
    }else if(h$page.no == "3"){
      load("./stat/plot/Bar_Plot.RData")
      print(p())
    }else{
      load("./stat/plot/Pyramid_Plot.RData")
      print(p())
    }
  })

  gtkMain()

}



