anl_align_nomul <- function(h,...){
  insert(st,"Start : Alignment & Read Counting", do.newline = TRUE )
  dir <- svalue(fq_dir)
  #c <- svalue(selfq_button)

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

  # if(
  #   identical(dir,character(0)) & identical(c,character(0))){
  #   gmessage("Warning : Select a directory of FASTQ files")}
  # else{
    # if(identical(dir,character(0))){
    #   dir <- svalue(selfq_button)}else{dir <- getwd()}


    # if(svalue(gsel_button) == "Homo sapiens (hg38 - GRCh38 Dec 2013)"){
    #   ref_name <- "Hsapi38"
    # }

    ref_name <- cl_name$fa[grep(svalue(gsel_button),cl_name$P4)]

    url <- file.path("http://treasure.pmrc.re.kr/data/genome", svalue(ref_P1), ref_name, paste0(ref_name, ".zip"))
    dw_refer(url, ref_name, paste0(ref_name, ".zip"))
    genome = system.file( "extdata", file.path("refer",ref_name, paste0(ref_name, "_artificial.fa")), package = "tReasure", mustWork = TRUE)

    sFile <- file.path(dir,"pre", "sample_trim.txt")
    sFile2 <- read.delim(sFile)

    # alignment future----------------------
    bowtie <- function(sFile, genome){
      bpid <- Sys.getpid()
      save(bpid, file = file.path(dir,"bpid"))
      proj <- qAlign(sFile, genome = genome, aligner = "Rbowtie",alignmentParameter = "-v 3 --best", clObj = NULL)
      return(proj)
    }

    #bow <- future(bowtie())
    # alignment status----------------------

    # algin_status <- function(){
    #   for(i in sFile2$SampleName){
    #     a <- sFile2$FileName[grep(i, sFile2$SampleName)]
    #     Sys.sleep(10)
    #     insert(st, paste("pre   : ", a), do.newline = TRUE)
    #     repeat{
    #       Sys.sleep(1)
    #       insert(st, ".", do.newline = FALSE)
    #       dir <- getwd()
    #       list <- list.files(file.path(dir,"pre"), pattern=".bam$", full.names = TRUE)
    #       if(length(grep(i, list))>0)break
    #     }
    #     insert(st, " ", do.newline = TRUE)
    #     a <- gsub(".fastq$","_*.bam",a)
    #   }
    #   insert(st, paste0("complete"), do.newline=TRUE)
    # }
    #
    #algin_status()======================

    # premmapping
    projt_pre <- bowtie(sFile, genome)
    alinstat <- data.frame(Names=row.names(alignmentStats(projt_pre)),alignmentStats(projt_pre))
    write.table(alinstat, "./pre/preproccessing_align_stat.txt",sep = "\t")

    # QC
    #insert(st,"Start : QC", do.newline = TRUE )
    save(projt_pre, file = "./pre/premappingQC.RData")
    # qQCReport(projt_pre, pdfFilename = "./pre/qc_report.pdf", useSampleNames = TRUE)

    # remove reads(premature-tRNAs, non_tRNAs)
    insert(st, ".", do.newline = TRUE)
    insert(st,"Removing reads aligned premauture tRNA and no tRNAs", do.newline = TRUE )
    anl_rm(ref_name)


    # postmapping
    insert(st, ".", do.newline = TRUE)
    insert(st,"postprocessing", do.newline = TRUE )
    clu_genome = system.file( "extdata", file.path("refer", ref_name, paste0(ref_name, ".tRNAscan_mature.fa")), package = "tReasure", mustWork = TRUE)

    matfq_list <- list.files(file.path(dir, "post"), pattern = "_trim_mature.fastq$", full=TRUE)
    matsFile <- data.frame(FileName = matfq_list, SampleName = gsub("_trim_mature.fastq","",basename(matfq_list)))
    write.table(matsFile, file.path(dir, "post", "sample_mature.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
    matsFile <- file.path(dir,"post", "sample_mature.txt")

    projt_post <- bowtie(matsFile,clu_genome)
    alinstat <- data.frame(Names=row.names(alignmentStats(projt_post)),alignmentStats(projt_post))
    write.table(alinstat, "./post/postprocessing_align_stat.txt",sep = "\t")

    # Readcounting
    insert(st,"counting", do.newline = TRUE )
    anl_rc(ref_name)
    insert(st,"Done : Alignment & Read Counting. Click! Next tab of Filtering.", do.newline = TRUE )
    insert(st, ".", do.newline = TRUE)

  }
#}
