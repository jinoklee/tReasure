#' window
#'
#' @param ()
#' @return anl_align
#' @export
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