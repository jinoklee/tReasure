anl_quantile <- function(h,...){
  insert(st,"Statistical analysis  : Quqntile Normalization using limma voom ", do.newline = TRUE )
  insert(st,"It takes a few minutes. Please wait....", do.newline = TRUE )
  dir <- getwd()
  
  tcount <- read.delim(file.path(dir, "rc", "filtered_readcount_trnas.txt"))
  ccount <- read.delim(file.path(dir, "rc", "filtered_readcount_isodecoders.txt"))
  acount <- read.delim(file.path(dir, "rc", "filtered_readcount_isoacceptors.txt"))
  
  t_out <- quant(tcount)
  c_out <- quant(ccount)
  a_out <- quant(acount)
  
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