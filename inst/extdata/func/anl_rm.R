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
