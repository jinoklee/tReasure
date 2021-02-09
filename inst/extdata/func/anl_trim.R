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
    sFile2 <- sFile1[,c("FileName","SampleName")]

    sFile2$FileName <- file.path(dir, 'pre', sFile2$SampleName)
    sFile2$FileName <- gsub("$",".fastq", sFile2$FileName)
    sFile2$FileName <- sub(paste0(".","fastq$"), paste0("_trim.","fastq"),sFile2$FileName)

    write.table(sFile2, sam_trim, sep = "\t", quote = FALSE, row.names = FALSE)

    # preprocessing future-------------------

    pre <- function(){
      apid <- Sys.getpid()
      save(apid, file = file.path(dir,"apid"))
      resa <- preprocessReads(filename = sFile1$FileName,outputFilename = sFile2$FileName,
                              minLength = minLength, Rpattern = aseq)
      return(resa)}

    preno <- function(){resno <- preprocessReads(filename = sFile1$FileName, outputFilename = sFile2$FileName,
                                                 minLength = minLength)
    return(resno)}

    if(svalue(adapt) == "Illumina smallRNA 3' adapter"){
      aseq <- "TGGAATTCTCGGGTGCCAAGG"
      minLength = as.numeric(svalue(min))
      a <- future(pre())
    }else if(svalue(adapt) =="Illumina universal adapter"){
      aseq <- "AGATCGGAAGAG"
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
        a <- sFile1$FileName[grep(i, sFile2$SampleName)]
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
    insert(st, "", do.newline = TRUE)
    insert(st,"Done : Quality Control. Click! Next tab of Alignment & Counting.", do.newline = TRUE )
    insert(st, ".", do.newline = TRUE)
  }
}
