
pkg <- c("gWidgets2","gWidgets2RGtk2","cairoDevice","plotrix","tidyverse",
         "gridExtra","ggplot2","grid","dplyr","statmod",
         "future", "stringr")


ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE, repos="http://cran.us.r-project.org")
  sapply(pkg, require, character.only = TRUE)
}

ipak(pkg)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

bpkg <- c("QuasR","DESeq2","edgeR", "Rsamtools","seqinr","ShortRead")

ibpak <- function(bpkg){
  new.pkg <- bpkg[!(bpkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    BiocManager::install(new.pkg)
  sapply(bpkg, require, character.only = TRUE)
}

ibpak(bpkg)
