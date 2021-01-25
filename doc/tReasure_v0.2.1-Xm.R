#install.packages("statmod")
options(guiToolkit="RGtk2")
library(gWidgets2)
library(gWidgets2RGtk2)
library(cairoDevice)
library(QuasR)
library(Rsubread)
library(edgeR)
library(ggplot2)
library(grid)
library(gridExtra)
library(plotrix)
library(tidyverse)
library(stringr)
library(Rsamtools)
library(seqinr)
library(ShortRead)
library(tidyverse)
library(dplyr)
library(DESeq2)
library(future)
plan(multiprocess)
library(statmod)
#plan(multicore)


#......................................................................................#
# setting the PATH for TEST : system.file('', package="tReasure")
#......................................................................................#

intro <- system.file("extdata", "intro.png", package = "tReasure", mustWork = TRUE)

#......................................................................................#
# function
#......................................................................................#
source(system.file("extdata", "func/stopFuture.R", package = "tReasure", mustWork = TRUE))
source(system.file("extdata", "func/deseq_function.R", package = "tReasure", mustWork = TRUE))
source(system.file("extdata", "func/edge_function.R", package = "tReasure", mustWork = TRUE))
source(system.file("extdata", "func/pyramid.R", package = "tReasure", mustWork = TRUE))
source(system.file("extdata", "func/volcano.R", package = "tReasure", mustWork = TRUE))
source(system.file("extdata", "func/remove_reads.R", package = "tReasure", mustWork = TRUE))
cl_name <- read.table("/data6/jinoklee/tReasure/inst/extdata/class_name.txt", sep = "\t", fill = T,header = T)

#......................................................................................#
#  tReasure
#......................................................................................#

# Main window----------------------
window <- gwindow("tReasure")
addHandlerUnrealize(window, handler = function(h,...) {
  val <- gconfirm("Really close window", parent=h$obj)
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
gr5 <- ggroup(container = notebook, label ="  DEtRNA Detection ",horizontal = FALSE)
gr6 <- ggroup(container = notebook, label ="  Visualiztion ",horizontal = FALSE)
child <- gframe("  Progress  ", container = mother, horizontal = FALSE)
size(child) <- c(1000,100)
st <- gtext("  ",container = child)
size(st) <- c(1000,100)

# gr0 : Introduction---------------
ggr1 <- ggroup(container = gr0, horizontal = TRUE, fill = "both", expand = TRUE)
gimage(intro, container = ggr1)

# gr1 : Sample List----------------

# design-------------------
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


make_button <- gbutton("RUN", container = tmp.1)

# option-------------------

tmp.11 <- gframe(" [*OPTION] Upload the sample list ", container = paned.1, horizontal = FALSE, spacing = 10) ; size(tmp.11) <- c(250,200)
addSpace(tmp.11, 10)

glabel(" Select a file of \n the pre-made sample list : ", container = tmp.11, anchor = c(-1,0))
addSpace(tmp.11, 20)

sel_button <- gfilebrowse(text=" ", quote = FALSE, type="open", container = tmp.11,
                          filter=list("*.txt" = list(patterns = c("*.txt")), "*.*" = list(patterns = c("*"))))

tbl <- gdf(data.frame(FileName=c("Path/samplename.fastq"),SampleName=c("samplename"),
                      Group=c("control"), Batch=c("batch I")), container = ggr1, multiple = FALSE)

# handler-------------------
addHandlerChanged(fq_dir, handler = function(h,...){
  setwd(svalue(fq_dir))
  dir <- getwd()
  if(!dir.exists(paste(dir, "/trim", sep = ""))){
    dir.create(paste(dir, "/trim", sep = ""), recursive = TRUE)}
  if(!dir.exists(paste(dir, "/stat/plot", sep = ""))){
    dir.create(paste(dir, "/stat/plot", sep = ""), recursive = TRUE)}
})


addHandlerChanged(make_button, handler = function(h,...){
  b <- svalue(fq_dir)
  if(identical(b,character(0))){
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
    if(length(fq_list)%%2 == 0){len = length(fq_list)/2}else{len=round(length(fq_list)/2)}
    sname <- gsub(paste(".","fastq", sep = ""), "", sfile)

    sample <- data.frame(FileName = fq_list,
                         SampleName = sname, Group = c(rep("control", length(fq_list)-len), rep("test", len)),
                         Batch = c(rep("NA", length(fq_list))))

    sample$Batch <- as.character(sample$Batch)
    tbl[] <- sample

    SampleFile <- file.path(svalue(fq_dir), "sample.txt")
    write.table(sample, SampleFile, sep = "\t", quote = FALSE, row.names = FALSE)


    new_sample <- tbl[]
    write.table(new_sample, SampleFile, sep = "\t", quote = FALSE, row.names = FALSE)
    insert(st, " ", do.newline = TRUE)
    insert(st, "Complete making the sample list.", do.newline = TRUE)
    if(length(grep("control", new_sample$Group)) < 2 | length(grep("test", new_sample$Group)) < 2){
      insert(st, "Warnning : There should be at least 2 samples per group for differentially expression analysis.")
    }
    insert(st, ".", do.newline = TRUE)
  }
})


addHandlerChanged(sel_button, handler = function(h,...){
  sample2 <- read.table(svalue(sel_button), header = T, sep = "\t")
  sampath <- gsub("sample.txt", "", svalue(sel_button))
  setwd(sampath)
  dir <- getwd()
  if(!dir.exists(paste(dir, "/trim", sep = ""))){
    dir.create(paste(dir, "/trim", sep = ""), recursive = TRUE)}
  if(!dir.exists(paste(dir, "/stat/plot", sep = ""))){
    dir.create(paste(dir, "/stat/plot", sep = ""), recursive = TRUE)}

  sample2$Batch <- as.character(sample2$Batch)
  tbl[] <- sample2
  insert(st, " ", do.newline = TRUE)
  insert(st, "Complete uploading the sample list.", do.newline = TRUE)
  if(length(grep("control", sample2$Group)) < 2 | length(grep("test", sample2$Group)) < 2){
    insert(st, "Warnning : There should be at least 2 samples per group for differentially expression analysis.")}
  insert(st, "Click! Next tab of Quality Control.", do.newline = TRUE)
  # save update sample list
  addHandlerChanged(tbl, handler = function(h, ...){
    new_sample <- tbl[]
    SampleFile <- file.path(svalue(fq_dir), "sample.txt")
    write.table(new_sample, SampleFile, sep = "\t", quote = FALSE, row.names = FALSE)
    insert(st, " ", do.newline = TRUE)
    insert(st, "Complete making the sample list.", do.newline = TRUE)
    if(length(grep("control", new_sample$Group)) < 2 | length(grep("test", new_sample$Group)) < 2){
      insert(st, "Warnning : There should be at least 2 samples per group for differentially expression analysis.")
    }
    insert(st, ".", do.newline = TRUE)
  })
})

# gr2 : Preprocessing----------------

# design-------------------
addSpace(gr2, 20)
ggr21 <- ggroup(container = gr2,  spacing = 10); size(ggr21) <- c(896,532)
addSpace(ggr21, 10)
paned.2 <- gpanedgroup(container = ggr21, horizontal = FALSE, spacing = 10)
tmp.2 <- gframe("  Settings  ", container = paned.2, horizontal = FALSE, spacing = 10); size(tmp.2) <- c(250,532)
addSpace(tmp.2, 10)

# option-------------------
glabel(" Information of adapter : ", container = tmp.2, anchor = c(-1,0))
addSpace(tmp.2, 10)

adapt <- gradio(c("Illumina smallRNA 3' adapter", "Illumina universal adapter","SOLiD Adapter","No Adapter"), container = tmp.2)
addSpace(tmp.2, 10)

glabel(" Minimum lenght : ", container = tmp.2, anchor = c(-1,0))
addSpace(tmp.2, 10)

# min <- gradio(c("15", "16","17"), container = tmp.2)
min <- gradio(c("10"), container = tmp.2)
addSpace(tmp.2, 10)


glabel(" Minimum quality : ", container = tmp.2, anchor = c(-1,0))
addSpace(tmp.2, 10)

q <- gradio(c("25", "30"), container = tmp.2)
addSpace(tmp.2, 20)

trim_button <- gbutton("RUN", container = tmp.2, expand = FALSE)

# handler-------------------
addHandlerClicked(trim_button, handler = function(h,...){
  if(file.exists("sample.txt") == FALSE & file.exists("*.fastq") == FALSE){
    gmessage("Warning : No file found matching sample list or fastq files")
  }else{
    insert(st, "", do.newline = TRUE)
    insert(st,"Start : Quality Control", do.newline = TRUE )
    dir <- getwd()

    sam_trim <- file.path(dir,"trim","sample.txt")
    sam_trim <- sub(".txt", "_trim.txt", sam_trim)

    sFile1 <- read.delim("sample.txt", header = TRUE, as.is = TRUE)
    sFile2 <- sFile1[,c("FileName","SampleName")]

    sFile2$FileName <- file.path(dir, 'trim', sFile2$SampleName)
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
    write.table(rest, "./trim/trim_res.txt", sep="\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    rest <- rest[-c(2,5,6),]
    res_trim <- gtable(data.frame(rest), container = ggr21)
    insert(st, "", do.newline = TRUE)
    insert(st,"Done : Quality Control. Click! Next tab of Alignment & Counting.", do.newline = TRUE )
    insert(st, ".", do.newline = TRUE)
  }
})

# gr3 : Analysis-------------------

# design-------------------
addSpace(gr3, 20)
ggr31 <- ggroup(container = gr3,  spacing = 10); size(ggr31) <- c(896,532)
addSpace(ggr31, 10)
paned.3 <- gpanedgroup(container = ggr31, horizontal = FALSE, spacing = 10)
pangr <- ggroup(container = paned.3, horizontal =FALSE)
tmp.3 <- gframe("  Settings  ", container = pangr , horizontal = FALSE, spacing = 10); size(tmp.3) <- c(250,260)
addSpace(tmp.3, 10)

RC_view <- gnotebook(container = ggr31 ); size(RC_view) <- c(700, 500)
RC_trna <- ggroup(container = RC_view, horizontal = FALSE, label = " tRNA level ")
RC_codon <- ggroup(container = RC_view, horizontal = FALSE, label = " Isodecoder level")
RC_aa <- ggroup(container = RC_view, horizontal = FALSE, label = "Isoacceptor level")

# option-------------------
glabel(" Genome assembly : ", container = tmp.3, anchor = c(-1,0))
addSpace(tmp.3, 10)
glabel("  select phylogeny group : ", container = tmp.3, anchor = c(-1,0))
ref_P1 <- gcombobox(c("groups", "eukaryota" ,"archaea", "bacteria"),container = tmp.3)
lsNms <- function(d, envir=.GlobalEnv)  unique(cl_name$P2[grep(svalue(ref_P1), cl_name$P1)], envir=envir)
ref_P2 <- (vnames <- gcombobox("2nd groups",container = tmp.3))
ref_P3 <- (vnames <- gcombobox("3rd groups",container = tmp.3))

addHandlerChanged(ref_P1, handler = function(h,...){
  ref_P2[] <- lsNms(svalue(ref_P2))
})

addSpace(tmp.3, 10)
#glabel("  select subgroup : ", container = tmp.3, anchor = c(-1,0))
#ref_P2 <- gcombobox(" ",container = tmp.3)
#addSpace(tmp.3, 10)
#glabel("  select final genome : ", container = tmp.3, anchor = c(-1,0))
#ref <- gcombobox(c("Homo sapiens (GRCh37/hg19)", "Homo sapiens (GRCh38/hg38)"),container = tmp.3)
#addSpace(tmp.3, 10)
addSpring(tmp.3)
anl_button <- gbutton("RUN", container=tmp.3)

addSpace(pangr, 10)
tmp.31 <- gframe(" [*OPTION] Working directory ", container = pangr , horizontal = FALSE, spacing = 10); size(tmp.31) <- c(250,70)
addSpace(tmp.31, 10)

glabel(" Select the directory of the FASTQ files : ", container = tmp.31, anchor = c(-1,0))
addSpace(tmp.31, 10)

tmp.32 <- gframe(" [*OPTION] Load the readcounts files ", container = pangr , horizontal = FALSE, spacing = 10); size(tmp.32) <- c(250,120)
addSpace(tmp.32, 10)

glabel(" Select the file of readcounts : ", container = tmp.32, anchor = c(-1,0))
addSpace(tmp.32, 10)

selfq_button <- gfilebrowse(text = "", quote = FALSE, type = "selectdir", container = tmp.31,
                            filter=list("*.fa" = list(patterns = c("*.fa")), "*.*" = list(patterns = c("*"))))

seltrn_button <- gfilebrowse(text = "tRNA level", quote = FALSE, type = "open", container = tmp.32,
                             filter=list("*.txt" = list(patterns = c("*.txt")), "*.*" = list(patterns = c("*"))))

selco_button <- gfilebrowse(text = "Isodecoder level", quote = FALSE, type = "open", container = tmp.32,
                            filter=list("*.txt" = list(patterns = c("*.txt")), "*.*" = list(patterns = c("*"))))

selaa_button <- gfilebrowse(text = "Isoacceptor level", quote = FALSE, type = "open", container = tmp.32,
                            filter=list("*.txt" = list(patterns = c("*.txt")), "*.*" = list(patterns = c("*"))))


# handler --------------
addHandlerChanged(selfq_button, handler = function(h,...){
  setwd(svalue(selfq_button))
  dir <- getwd()
  if(!dir.exists(paste(dir, "/trim", sep = ""))){
    dir.create(paste(dir, "/trim", sep = ""), recursive = TRUE)}
  if(!dir.exists(paste(dir, "/stat/plot", sep = ""))){
    dir.create(paste(dir, "/stat/plot", sep = ""), recursive = TRUE)}
})

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

addHandlerClicked(anl_button,handler = function(h, ...){
  insert(st,"Start : Alignment & Counting of reads", do.newline = TRUE )
  b <- svalue(fq_dir)
  c <- svalue(selfq_button)

  # download_refer
  dw_refer <- function(url, sub, sub.zip){
    R_path <- system.file( "extdata", package = "tReasure", mustWork = TRUE)
    if(!dir.exists(file.path(R_path,"refer",sub))){
      dir.create(file.path(R_path,"refer",sub), recursive = TRUE)
    }
    Rsub_path <- file.path(R_path,"refer",sub)
    if(!file.exists(file.path(Rsub_path, paste(sub,"fa$",sep = ".")))){
      insert(st, "Download and unzip genome indexes. It takes several minutes. Please wait....", do.newline = TRUE)
      download.file(url, destfile = file.path(Rsub_path, sub.zip))
      unzip(zipfile = file.path(Rsub_path, sub.zip), exdir= Rsub_path)
      #file.remove(file.path(Rsub_path))
      insert(st, "Complete downloaded.", do.newline = TRUE)
    }
  }

  if(
    identical(b,character(0)) & identical(c,character(0))){
    gmessage("Warning : Select a directory of FASTQ files")}
  else{
  if(identical(b,character(0))){
    dir <- svalue(selfq_button)}else{dir <- getwd()}


    if(svalue(ref) == "Homo sapiens (GRCh38/hg38)"){
      ref_name <- "Hsapi38"
    }

    url <- file.path("http://treasure.pmrc.re.kr/data/genome/ucsc", svalue(ref_P1), ref_name, paste0(ref_name, ".zip"))
    dw_refer(url, ref_name, paste0(ref_name, ".zip"))
    genome = system.file( "extdata", file.path("refer",ref_name, paste0(ref_name, "_artificial.fa")), package = "tReasure", mustWork = TRUE)

    sFile <- file.path(dir,"trim", "sample_trim.txt")
    sFile2 <- read.delim(sFile)

    # alignment future----------------------

    # dir <- getwd()
    # bowtie <- function(){
    #   bpid <- Sys.getpid()
    #   save(bpid, file = file.path(dir,"bpid"))
    #   proj <- qAlign(sFile, genome = genome, aligner = "Rbowtie",alignmentParameter = "-v 2 --best", clObj = NULL)
    #   return(proj)
    # }
    #bow <- future(bowtie())


    # alignment status----------------------

    # algin_status <- function(){
    #   for(i in sFile2$SampleName){
    #     a <- sFile2$FileName[grep(i, sFile2$SampleName)]
    #     Sys.sleep(10)
    #     insert(st, paste("Aligning   : ", a), do.newline = TRUE)
    #     repeat{
    #       Sys.sleep(1)
    #       insert(st, ".", do.newline = FALSE)
    #       dir <- getwd()
    #       list <- list.files(file.path(dir,"trim"), pattern=".bam$", full.names = TRUE)
    #       if(length(grep(i, list))>0)break
    #     }
    #     insert(st, " ", do.newline = TRUE)
    #     a <- gsub(".fastq$","_*.bam",a)
    #   }
    #   insert(st, paste0("complete"), do.newline=TRUE)
    # }

    #algin_status()
    # projt <- bowtie()
    #
    # alinstat <- data.frame(Names=row.names(alignmentStats(projt)),alignmentStats(projt))
    # write.table(alinstat, "./trim/alignment_stat.txt",sep = "\t")


    # read counts

    bFile <- list.files(file.path(dir,"trim"), pattern = ".bam$", full.names = T)

    # tRNA level
    readcount <- function(trna, gene, output){
      trna <- featureCounts(bFile, annot.ext = annot_ext, isGTFAnnotationFile=TRUE, GTF.featureType = "gene", GTF.attrType = gene)
      tname <- data.frame(do.call('rbind',strsplit(as.character(colnames(trna$counts)), split="[.]")))
      colnames(trna$counts) <- tname$X1
      tcount <- data.frame(Names= row.names(trna$counts),trna$counts)
      write.table(tcount, file.path(dir,paste0("readcount_",output,".txt")), quote = F, sep = "\t", row.names = F)
      write.table(trna$annotation,file.path(paste0("annotation_",output,".txt")), quote = F, sep = "\t")
      return(tcount)
    }

    #tcount <- readcount(trna, "gene_id", "trna")
    #gtable(tcount, container=RC_trna)
    insert(st, " ", do.newline = TRUE)
    insert(st,"Done : Counting of individual tRNAs", do.newline = TRUE )

    #ccount <- readcount(decoder, "isodecoder", "isodecoder")
    #gtable(ccount, container=RC_codon)
    insert(st, " ", do.newline = TRUE)
    insert(st,"Done : Counting of isodecoders", do.newline = TRUE )

    #acount <- readcount(acceptor, "isoacceptor", "isoacceptor")
    #gtable(acount, container=RC_aa)
    insert(st,"Done : Counting of isoacceptors", do.newline = TRUE )

    insert(st, " ", do.newline = TRUE)
    insert(st,"Done : Alignment & Counting of reads. Click! Next tab of Filtering.", do.newline = TRUE )
    insert(st, ".", do.newline = TRUE)
  }

})

# gr4 : Preview--------------------

# design--------------------
# left
addSpace(gr4, 20)
ggr41 <- ggroup(container = gr4,  spacing = 10); size(ggr41) <- c(896,532)
addSpace(ggr41, 10)
paned.4 <- gpanedgroup(container = ggr41, horizontal = FALSE, spacing = 10)
tmp.4 <- gframe("  Preview  ", container = paned.4, horizontal = FALSE, spacing = 10)
size(tmp.4) <- c(250,284)
addSpace(tmp.4, 10)

tmp.41 <- gframe( " Settings : ", container = tmp.4, anchor =c(-1,0), horizontal = FALSE)
addSpace(tmp.41, 10)
prelyt <- gformlayout(container = tmp.41, spacing = 1.5)
cfilv <- gcombobox(seq(0,10,by=1),label = "cpm value >=  ", container = prelyt)
sfil <- gcombobox(seq(0,100, by=10), label = "sample (%) >=  ", container = prelyt)
addSpace(tmp.41, 10)

addSpace(tmp.4, 20)
Pre_button <- gbutton("RUN", container = tmp.4)

# right
Pre_view <- gnotebook(container = ggr41 ); size(Pre_view) <- c(700, 532)
Pre_MDS <- ggroup(container = Pre_view, horizontal = FALSE, label=" MDS Plot ")

tmp.42 <- gframe(" [*OPTION] Working directory ", container = paned.4, horizontal = FALSE, spacing = 10)
size(tmp.42) <- c(250,200)
addSpace(tmp.42, 10)
glabel(" Select the directory of the readcounts files :", container = tmp.42, anchor = c(-1,0))
addSpace(tmp.42, 10)

# option--------------------
selrc_button <- gfilebrowse(text = "", quote = FALSE, type = "selectdir", container = tmp.42,
                            filter=list("*.txt" = list(patterns = c("*.txt")), "*.*" = list(patterns = c("*"))))

# handler--------------------
addHandlerChanged(Pre_button, handler = function(h, ...){
  b <- svalue(selrc_button)
  if(identical(b,character(0))){invisible()
  }else{
    setwd(svalue(selrc_button))
    dir <- getwd()
    if(!dir.exists(paste(dir, "/trim", sep = ""))){
      dir.create(paste(dir, "/trim", sep = ""), recursive = TRUE)}
    if(!dir.exists(paste(dir, "/stat/plot", sep = ""))){
      dir.create(paste(dir, "/stat/plot", sep = ""), recursive = TRUE)}
  }

  sFile <- read.delim("sample.txt")
  tcount<- read.delim("readcount_trna.txt")
  ccount <- read.delim("readcount_isodecoder.txt")
  acount <- read.delim("readcount_isoacceptor.txt")

  # tRNA levels
  filter_gene <- function(count){
    rownames(count) <- count$Names
    count<- count[,colnames(count) %in% sFile$SampleName]
    x <- DGEList(count, group = sFile$Group)
    isexpr <- rowSums(cpm(x) > 1) >= 90*nrow(sFile)/100
    x <- x[isexpr,]
    return(x)
  }

  filter_t <- filter_gene(tcount)
  filter_c <- filter_gene(ccount)
  filter_a <- filter_gene(acount)

  filter_tt <- data.frame(Names = rownames(filter_t), filter_t$counts)
  filter_cc <- data.frame(Names = rownames(filter_c), filter_c$counts)
  filter_aa <- data.frame(Names = rownames(filter_a), filter_a$counts)

  write.table(filter_tt, "filtered_readcount_trna.txt",sep = "\t", quote = F, row.names = F)
  write.table(filter_cc, "filtered_readcount_isodecoder.txt",sep = "\t", quote = F, row.names = F)
  write.table(filter_aa, "filtered_readcount_isoacceptor.txt",sep = "\t", quote = F, row.names = F)

  y <- calcNormFactors(filter_t)
  cols <- as.character(sFile$Group)
  cols[cols=="control"] = "blue"
  cols[cols=="test"] ="red"
  png( './stat/plot/plotMDS.png', width=650, height=500,  res=80)
  plotMDS(y,col=cols)
  graphics.off()
  gimage("./stat/plot/plotMDS.png", container = Pre_MDS)

  insert(st, " ", do.newline = TRUE)
  insert(st,"Done : Filtering. Click! Next tab of DEtRNA Detection.", do.newline = TRUE )
  insert(st, ".", do.newline = TRUE)
})



# gr5 : DEtRNA list----------------

# Statistical design--------------------
addSpace(gr5, 20)
ggr51 <- ggroup(container = gr5,  spacing = 10); size(ggr51) <- c(896,532)
addSpace(ggr51, 10)
paned.5 <- gpanedgroup(container = ggr51, horizontal = FALSE, spacing = 10)

ggr52 <- gnotebook(container=paned.5, tab.pos  = 3) ; size(ggr52) <- c(250, 200)
stat2 <- ggroup(container = ggr52, label = " DEseq2 ")
stat1 <- ggroup(container = ggr52, label = " EdgeR ")

tmp.51 <- gframe("  Statics  ", container = stat1, horizontal = FALSE, spacing = 10); size(tmp.51) <- c(235,150)
addSpace(tmp.51, 10)
glabel(" Stat.Method : ", container = tmp.51, anchor = c(-1,0))
addSpace(tmp.51, 10)
method_edgeR <- gradio(c("Exact test", "Likelihood F-test", "Quasi-likelihood F-test"), container = tmp.51)
addSpace(tmp.51, 10)

tmp.52 <- gframe("  Statics  ", container = stat2, horizontal = FALSE, spacing = 10); size(tmp.52) <- c(235,150)
addSpace(tmp.52, 10)
glabel(" Stat.Method : ", container = tmp.52, anchor = c(-1,0))
addSpace(tmp.52, 10)
method_deseq <- gradio("Walt test", container = tmp.52)
addSpace(tmp.52, 10)

# Plot design--------------------
addSpace(gr6, 20)
ggr61 <- ggroup(container = gr6,  spacing = 10); size(ggr61) <- c(996,550)
addSpace(ggr61, 10)
Plot_view <- gnotebook(container = ggr61 ); size(Plot_view) <- c(950, 550)
Plot_trna <- ggroup(container = Plot_view, horizontal = FALSE, label = " tRNA level ")
Plot_codon <- ggroup(container = Plot_view, horizontal = FALSE, label = " Isodecoder level ")
Plot_aa <- ggroup(container = Plot_view, horizontal = FALSE, label = " Isoacceptor level ")

tmp.61 <- gframe("",container = Plot_trna, anchor =c(-1,0), horizontal = TRUE)
tmp.611 <- ggroup(container = tmp.61, horizontal = TRUE)
wid_label <-glabel(" Width : ", container = tmp.611, anchor=c(-1,0)); size(wid_label) <- c(110,20)
wid_v <- gtext("650",container = tmp.611, width = 120, height = 20 )
tmp.612 <- ggroup(container = tmp.61, horizontal = TRUE)
hei_label<- glabel(" Height : ", container = tmp.612, anchor = c(-1,0)); size(hei_label) <- c(110,20)
hei_v <- gtext("460",container = tmp.612, width = 120, height = 20)
tmp.613 <- ggroup(container = tmp.61, horizontal = TRUE)
res_label<- glabel(" Resolution : ", container = tmp.613, anchor = c(-1,0)); size(res_label) <- c(110,20)
res_v <- gtext("100",container = tmp.613, width = 120, height = 20)


tmp.62 <- gframe(container = Plot_codon, anchor =c(-1,0), horizontal = TRUE)
tmp.621 <- ggroup(container = tmp.62, horizontal = TRUE)
wid_label <-glabel(" Width : ", container = tmp.621, anchor=c(-1,0)); size(wid_label) <- c(110,20)
wid_b <- gtext("850",container = tmp.621, width = 120, height = 20 )
tmp.622 <- ggroup(container = tmp.62, horizontal = TRUE)
hei_label<- glabel(" Height : ", container = tmp.622, anchor = c(-1,0)); size(hei_label) <- c(110,20)
hei_b <- gtext("460",container = tmp.622, width = 120, height = 20)
tmp.623 <- ggroup(container = tmp.62, horizontal = TRUE)
res_label<- glabel(" Resolution : ", container = tmp.623, anchor = c(-1,0)); size(res_label) <- c(110,20)
res_b <- gtext("100",container = tmp.623, width = 120, height = 20)


tmp.63 <- gframe(container = Plot_aa, anchor =c(-1,0), horizontal = TRUE)
tmp.631 <- ggroup(container = tmp.63, horizontal = TRUE)
wid_label <-glabel(" Width : ", container = tmp.631, anchor=c(-1,0)); size(wid_label) <- c(110,20)
wid_p <- gtext("650",container = tmp.631, width = 120, height = 20 )
tmp.632 <- ggroup(container = tmp.63, horizontal = TRUE)
hei_label<- glabel(" Height : ", container = tmp.632, anchor = c(-1,0)); size(hei_label) <- c(110,20)
hei_p <- gtext("460",container = tmp.632, width = 120, height = 20)
tmp.633 <- ggroup(container = tmp.63, horizontal = TRUE)
res_label<- glabel(" Resolution : ", container = tmp.633, anchor = c(-1,0)); size(res_label) <- c(110,20)
res_p <- gtext("100",container = tmp.633, width = 120, height = 20)


v_plot <- ggraphics(container = Plot_trna)
b_plot <- ggraphics(container = Plot_codon)
p_plot <- ggraphics(container = Plot_aa)

# Statistical Option --------------------
statedgeR_button <- gbutton ("RUN EdgeR", container = tmp.51)

statdeseq_button <- gbutton("RUN DESeq2", container = tmp.52)

# Statistical handler --------------------
addHandlerClicked(statedgeR_button, handler = function(h, ...){
  insert(st,"Statistical analysis  : EdgeR", do.newline = TRUE )
  insert(st,"It takes a few minutes. Please wait....", do.newline = TRUE )

  tcount <- read.delim("filtered_readcount_trna.txt")
  ccount <- read.delim("filtered_readcount_isodecoder.txt")
  acount <- read.delim("filtered_readcount_isoacceptor.txt")

  t_out <- edge_function(tcount)
  c_out <- edge_function(ccount)
  a_out <- edge_function(acount)

  stat_trna[] <- t_out
  stat_codon[] <- c_out
  stat_aa[] <- a_out

  write.table(t_out, "./stat/stat_trna_list.txt", sep="\t", quote = FALSE, row.names = F)
  write.table(c_out, "./stat/stat_isodecoder_list.txt", sep="\t", quote = FALSE, row.names = F)
  write.table(a_out, "./stat/stat_isoacceptor_list.txt", sep="\t", quote = FALSE, row.names = F)
  insert(st, " ", do.newline = TRUE)
  insert(st,"Done : Statistical analysis. Set the threshold value.", do.newline = TRUE )
  insert(st, ".", do.newline = TRUE)
})

addHandlerClicked(statdeseq_button, handler = function(h,...){
  insert(st,"Statistical analysis  : DESeq2", do.newline = TRUE )
  dir <- getwd()
  if(!dir.exists(paste(dir, "/stat/plot", sep = ""))){
    dir.create(paste(dir, "/stat/plot", sep = ""), recursive = TRUE)}

  tcount <- read.delim("filtered_readcount_trna.txt")
  ccount <- read.delim("filtered_readcount_isodecoder.txt")
  acount <- read.delim("filtered_readcount_isoacceptor.txt")

  t_out <- deseq_function(tcount)
  c_out <- deseq_function(ccount)
  a_out <- deseq_function(acount)

  stat_trna[] <- t_out
  stat_codon[] <- c_out
  stat_aa[] <- a_out

  write.table(t_out, "./stat/stat_trna_list.txt", sep="\t", quote = FALSE, row.names = F)
  write.table(c_out, "./stat/stat_isodecoder_list.txt", sep="\t", quote = FALSE, row.names = F)
  write.table(a_out, "./stat/stat_isoacceptor_list.txt", sep="\t", quote = FALSE, row.names = F)
  insert(st, " ", do.newline = TRUE)
  insert(st,"Done : Statistical analysis. Set the threshold value.", do.newline = TRUE )
  insert(st, ".", do.newline = TRUE)
})

# Threshold Design --------------------
widget_list <- list()
ggr.53 <- gframe("  Threshold  ", container = paned.5, horizontal = FALSE, spacing = 10)
addSpace(ggr.53, 10)
fdr_label<- glabel(" Adjust.Method : ", container = ggr.53, anchor = c(-1,0))
addSpace(ggr.53, 10)
widget_list$fdr_s <- gcombobox(c("FDR", "Bonferroni","Benjamini-Hochberg"  ), container = ggr.53)
addSpace(ggr.53, 10)

statlyt <- gformlayout(container = ggr.53, spacing = 1.5)
widget_list$pval <- gcombobox(c(0,0.001, 0.05, 0.01),label = "Adj P.value  < : ", container = statlyt)
widget_list$FC <- gcombobox(c(0,1,1.5,2),label = "Fold change  > : ", container = statlyt)

addSpace(ggr.53, 20)

DE_view <- gnotebook(container = ggr51 ); size(DE_view) <- c(700, 532)
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


# Plot option --------------------
DE_button <- gbutton(" Plot ", container = ggr.53)

# Plot handler --------------------
addHandlerClicked(DE_button, handler = function(h, ...){
  DEtRNA <- stat_trna[]
  DEcodon <- stat_codon[]
  DEaa <- stat_aa[]

  write.table(DEtRNA, "./stat/DEtrna_list.txt", sep = "\t", quote = FALSE)
  write.table(DEcodon, "./stat/DEisodecoder_list.txt", sep = "\t", quote = FALSE)
  write.table(DEaa, "./stat/DEisoacceptor_list.txt", sep = "\t", quote = FALSE)

  volcano(650, 460, 80)
  v_plot <- gimage("./stat/plot/volcanoplot_trna01.png", container = Plot_trna)

  barplot(850, 460,80)
  b_plot <- gimage("./stat/plot/barplot_isodecoder01.png", container = Plot_codon)

  pyramid(650, 460,80)
  p_plot <- gimage("./stat/plot/pyramid_isoaccepter01.png", container = Plot_aa)

  insert(st, " ", do.newline = TRUE)
  insert(st,"Done : The plots of DEtRNA Detection. Click! Next tab of Visualization.", do.newline = TRUE )
  insert(st, ".", do.newline = TRUE)
})


# gr6 : Plots----------------------

# re-size saved plot handler----------------------
tsave <- gbutton(" Save ", container = tmp.61,expand=TRUE, fill="x", handler= function(h, ...){
  volcano(as.numeric(svalue(wid_v)), as.numeric(svalue(hei_v)), as.numeric(svalue(res_v)))
}); size(tsave) <- c(200,20)

csave <- gbutton(" Save ", container = tmp.62,expand=TRUE, fill="x",handler= function(h, ...){
  barplot(as.numeric(svalue(wid_b)), as.numeric(svalue(hei_b)), as.numeric(svalue(res_b)))
}); size(csave) <- c(200,20)

asave <- gbutton(" Save ", container = tmp.63,expand=TRUE, fill="x", handler= function(h, ...){
  pyramid(as.numeric(svalue(wid_p)), as.numeric(svalue(hei_p)), as.numeric(svalue(res_p)))
}); size(asave) <- c(200,20)


# Main window----------------------

gtkMain()



