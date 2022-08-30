#' window
#'
#' @param ()
#' @return tReasure install and load
#' @export

install.tReasure <- function(){
  # install cran pakcage
  pkg <- c("gWidgets2","plotrix","tidyverse",
           "gridExtra","ggplot2","grid","dplyr","statmod",
           "future", "stringr")
  ipak <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg))
      install.packages(new.pkg, dependencies = TRUE, repos = "http://cran.r-project.org")
    #sapply(pkg, require, character.only = TRUE)
  }
  ipak(pkg)
  
  # install RGTK2 GUI
  if(Sys.info()[names(Sys.info())== "sysname"] == "Linux"){
    if(!"RGtk2"%in%installed.packages()[,"Package"]){
      install.packages(
        "https://cran.microsoft.com/snapshot/2021-11-08/src/contrib/RGtk2_2.20.36.2.tar.gz",repos=NULL) 
    }
    if(!"gWidgets2RGtk2"%in%installed.packages()[,"Package"]){
      install.packages(
        "https://cran.microsoft.com/snapshot/2021-11-08/src/contrib/gWidgets2RGtk2_1.0-7.tar.gz",repos=NULL)
    }
    if(!"cairoDevice"%in%installed.packages()[,"Package"]){
      install.packages(
        "https://cran.microsoft.com/snapshot/2021-11-08/src/contrib/cairoDevice_2.28.2.1.tar.gz",repos=NULL)
    }
    if(!"gWidgets2"%in%installed.packages()[,"Package"]){
      install.packages(
        "https://cran.microsoft.com/snapshot/2021-11-08/src/contrib/gWidgets2_1.0-8.tar.gz", repos=NULL)
    }
  }else if(Sys.info()[names(Sys.info())== "sysname"] == "Windows"){
    if(grepl("3.6",R.version.string)){
      if(!"RGtk2"%in%installed.packages()[,"Package"]){
        install.packages(
          "https://cran.microsoft.com/snapshot/2021-11-08/bin/windows/contrib/3.6/RGtk2_2.20.36.zip",repos=NULL)
      }
      if(!"gWidgets2RGtk2"%in%installed.packages()[,"Package"]){
        install.packages(
          "https://cran.microsoft.com/snapshot/2021-11-08/bin/windows/contrib/3.6/gWidgets2RGtk2_1.0-7.zip",repos=NULL)
      }
      if(!"cairoDevice"%in%installed.packages()[,"Package"]){
        install.packages(
          "https://cran.microsoft.com/snapshot/2021-11-08/bin/windows/contrib/3.6/cairoDevice_2.28.2.zip",repos=NULL)}
      if(!"gWidgets2"%in%installed.packages()[,"Package"]){
        install.packages(
          "https://cran.microsoft.com/snapshot/2021-11-08/bin/windows/contrib/3.6/gWidgets2_1.0-8.zip",repos=NULL)}
      }else{
        ver <- regmatches(R.version.string, regexpr("\\d{1}\\.\\d{1}", R.version.string))
        if(!"RGtk2"%in%installed.packages()[,"Package"]){
          install.packages(paste0("https://cran.microsoft.com/snapshot/2021-11-08/bin/windows/contrib/", ver,"/RGtk2_2.20.36.2.zip") , repos= NULL)
        }
        if(!"gWidgets2RGtk2"%in%installed.packages()[,"Package"]){
          install.packages(paste0("https://cran.microsoft.com/snapshot/2021-11-08/bin/windows/contrib/", ver,"/gWidgets2RGtk2_1.0-7.zip") , repos= NULL)
        }
        if(!"cairoDevice"%in%installed.packages()[,"Package"]){
          install.packages(paste0("https://cran.microsoft.com/snapshot/2021-11-08/bin/windows/contrib/", ver,"/cairoDevice_2.28.2.1.zip") , repos= NULL)
        }
        
        if(!"gWidgets2"%in%installed.packages()[,"Package"]){
          install.packages(paste0("https://cran.microsoft.com/snapshot/2021-11-08/bin/windows/contrib/", ver,"/gWidgets2_1.0-8.zip ") , repos= NULL)}
      }
    }else{
      ver <- regmatches(R.version.string, regexpr("\\d{1}\\.\\d{1}", R.version.string))
      if(ver != "4.2"){
        if(!"RGtk2"%in%installed.packages()[,"Package"]){
          install.packages("https://cran.r-project.org/src/contrib/Archive/RGtk2/RGtk2_2.20.36.3.tar.gz",repos=NULL, type = "source")
        }
        if(!"gWidgets2RGtk2"%in%installed.packages()[,"Package"]){
          install.packages("https://cran.r-project.org/src/contrib/Archive/gWidgets2RGtk2/gWidgets2RGtk2_1.0-7.tar.gz",repos=NULL, type = "source")
        }
        if(!"cairoDevice"%in%installed.packages()[,"Package"]){
          install.packages("https://cran.r-project.org/src/contrib/Archive/cairoDevice/cairoDevice_2.28.2.2.tar.gz",repos=NULL, type = "source")
        }
        }
        if(!"gWidgets2"%in%installed.packages()[,"Package"]){
          install.packages("https://cran.r-project.org/src/contrib/Archive/gWidgets2/gWidgets2_1.0-8.tar.gz",repos=NULL, type = "source")
        }
        }else{
          
          if(!"RGtk2"%in%installed.packages()[,"Package"]){
            install.packages("https://cran.microsoft.com/snapshot/2021-11-08/bin/macosx/contrib/r-release/RGtk2_2.20.36.2.tgz",
                             repos=NULL)
          }
          if(!"gWidgets2RGtk2"%in%installed.packages()[,"Package"]){
            install.packages("https://cran.microsoft.com/snapshot/2021-11-08/bin/macosx/contrib/r-release/gWidgets2RGtk2_1.0-7.tgz",
                             repos=NULL)
          }
          if(!"cairoDevice"%in%installed.packages()[,"Package"]){
            install.packages("https://cran.microsoft.com/snapshot/2021-11-08/bin/macosx/contrib/r-release/cairoDevice_2.28.2.1.tgz",
                             repos=NULL)
          }
          if(!"gWidgets2"%in%installed.packages()[,"Package"]){
            install.packages("https://cran.microsoft.com/snapshot/2021-11-08/bin/macosx/contrib/r-release/gWidgets2_1.0-8.tgz",
                             repos=NULL)
          }
          }
  }
  
  # install Biomanger
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  bpkg <- c("pillar","QuasR","DESeq2","edgeR", "Rsamtools","seqinr","ShortRead")
  
  ibpak <- function(bpkg){
    new.pkg <- bpkg[!(bpkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg))
      BiocManager::install(new.pkg)
    # sapply(bpkg, require, character.only = TRUE)
  }
  
  ibpak(bpkg)
  library(rstudioapi)
  restartSession(command = "tReasure::tReasure()")
  
}