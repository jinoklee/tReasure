
# tReasure
tReasure (tRna Expression Analysis Software Utilizing R for Easy use) is user-friendly tool for the tRNA expression analysis of deep-sequencing experiments for small RNAs using R packages. 

    tReasure currently implements the following tasks:
    1) Making sample list for analysis
    2) Pre-processing of trimming adapter and filtering reads
    3) Alignment and quantification
    4) Statistical analysis for differentially expressed tRNAs gene
    5) Visualization 

<br/><br/>
  ### Anlaysis flow
   ![Flow](https://github.com/jinoklee/tReasure/blob/master/inst/extdata/flow.png?raw=true)<br/><br/>
   



Installation guide
---------------------------
**tReasure** is a package for the R computing environment and it is assumed that you have already installed R. See the R project at (http://www.r-project.org). tReasure requires the gwidget2 graphical library to run and a few additional packages for the analysis of RNA-seq. 


  
1.  Install R from http://cran.r-project.org/ according to your operating system.
2.  Download tReasure packages from https://github.com/jinoklee/tReasure and additional files.  
**STEP1. Install the devtools package**
To install tReasure package, start by installing the devtools package.
Open R or Rstudio and type on console
'''
install.packages(“devtools”)
library(devtools)
install_github(“jinoklee/tReasure”)
library(tReasure)
'''
You can also install tReasure packages from local source.
    
Download tReasure_1.0.0.tar.gz
Open R or Rstudio and type on console
    
'''
install.packages(“~/Download/tReasure_1.0.0.tar.gz”, repo=NULL, type = “source”)
library(tReasure)
'''
**STEP2. Download files contained genome index and sample raw data. **
Before the start tReasure packages, you should download the genome index files what you want ( hg38, hg19, mm10 or all). 
* Download a bundle of genome index files and move the files inside tReasure package folder (~Documents/R/win-library/tReasure/WholeGenomeFasta). 
* Caution! Before moving, you should make a folder named “WholeGenomeFasta”.
You can also download sample raw files for a test. 

**STEP3. Download tReasure Rscript for analysis using Rscript execution.**
For a standalone tReasure for window user, download script’s file (tReasure.R). 
Download a bundle of files contained Rscript (tReasure.R) and a batch file named shortcut_install.bat.
Click on the shortcut_install.bat, it creates a shortcut icon of tReasure on the desktop.

3.  Preliminaries
3.1 install all the packages
'''
install.packages("gWidgets2")
install.packages("gWidgets2RGtk2")
install.packages("cairoDevice")
install.packages("plotrix")
install.packages("tidyverse")
install.packages("gridExtra")
install.packages("ggplot2")
install.packages("grid")
install.packages("dplyr")
install.packages("statmod")
install.packages(“future”)

if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("QuasR")
BiocManager::install("Rsubread")
BiocManager::install("DESeq2")
BiocManager::install("edgeR")
    
check that all the packages listed above have been installed correctly.
'''
options(guiToolkit="RGtk2")
library(gWidgets2)
library(gWidgets2RGtk2)
library(cairoDevice)
library(plotrix)
library(tidyverse)
library(gridExtra)
library(ggplot2)
library(grid)
library(dplyr)
library(QuasR)
library(Rsubread)
library(DESeq2)
library(edgeR)
library(statmod)
library(future)
plan(multiprocess) for Window or plan(multicore) for Linux
    
If you have successfully gone through the installation you are ready to use tReasure.

Start
---------------------------
1.  Using R or Rstudio
Open R or Rstudio and type on console
'''
library(tReasure)
tReasure()
'''
2. using Rscript 
2.1 For Window users 
Double click the icon of tReasure on the desktop.
2.2 For Linux or MacOS users
Click on the tReasure and type on command
'''
Rscript tReasure.R
'''

