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
