
<br/>
<br/>
<br/>
<br/>

# tReasure

**tReasure (tRna Expression Analysis Software Utilizing R for Easy use)** is user-friendly tool for the tRNA expression analysis of deep-sequencing experiments for small RNAs using R packages. 

    tReasure currently implements the following tasks:
    1) Making sample list for analysis
    2) Pre-processing of trimming adapter and filtering reads
    3) Alignment and quantification
    4) Statistical analysis for differentially expressed tRNAs gene
    5) Visualization 

<br/><br/>
  ## Analysis flow
   ![Flow](https://github.com/jinoklee/tReasure/blob/master/inst/extdata/flow.png?raw=true)<br/><br/>
   
<br/>
<br/>



## Installation guide

tReasure is a package for the R computing environment and it is assumed that you have already installed R according to your operating system. See the R project at (http://www.r-project.org). tReasure requires the gwidget2 graphical library to run and a few additional packages for the analysis of RNA-seq. 


<br/>

### **STEP1. Install the devtools package**

To install tReasure package, start by installing the devtools package. Open R or Rstudio and type on console as below.
<pre>
<code>
install.packages(“devtools”)
library(devtools)

install_github(“jinoklee/tReasure”) 
library(tReasure) 
</code>
</pre>

### **[OR]**
You can also install tReasure packages from local source. Download  [tReasure_1.0.0.tar.gz](https://www.dropbox.com/s/3vix75t5yh5okbp/tReasure_1.0.0.tar.gz?dl=0) and saved it on Download folders. Open R or Rstudio and type on console as below.
<pre>
<code>
install.packages(“~/Download/tReasure_1.0.0.tar.gz”, repo=NULL, type = “source”)
library(tReasure)
</code>
</pre>

<br/>

### **STEP2. Download files contained genome index and sample raw data**

Before the start tReasure packages, you should download the genome index files what you want ( [hg38](https://www.dropbox.com/s/f9qefc0sjxfzjkl/hg38.zip?dl=0), [hg19](https://www.dropbox.com/s/2wlfv8oh1qq7vzh/hg19.zip?dl=0), [mm10](https://www.dropbox.com/s/6z86iaszuvuijy7/mm10.zip?dl=0) or all). 
+ Download [a bundle of genome index files](https://www.dropbox.com/sh/1aikvdszjlvncic/AADzL8G55ayI3lRfzZ6LYjvPa?dl=0) and move the files inside tReasure package folder (~Documents/R/win-library/tReasure/WholeGenomeFasta)._Caution! Before moving, you should make a folder named “WholeGenomeFasta”._
+ You can also download [sample raw files](https://www.dropbox.com/s/iskwxurw06fjzkh/tReasure.R?dl=0) for a test. 

<br/>

### **STEP3. Download tReasure Rscript for analysis using Rscript execution**

For a standalone tReasure for window user, download script’s file. Download a bundle of files contained Rscript ([tReasure.R](https://www.dropbox.com/s/iskwxurw06fjzkh/tReasure.R?dl=0)) and a batch file named shortcut_install.bat.Click on the shortcut_install.bat, it creates a shortcut icon of tReasure on the desktop.

<br/>
<br/>
<br/>


## Preliminaries

Install the required packages for tReasure as below.
<pre>
<code>
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
</code>
</pre>

Check that all the packages listed above have been installed correctly as below.
<pre>
<code>
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
</code>
</pre>

If you have successfully gone through the installation you are ready to use tReasure.

<br/>
<br/>
<br/>


## Start

#### Using R or Rstudio

Open R or Rstudio and type on console as below.
<pre>
<code>
library(tReasure)
tReasure()
</code>
</pre>

#### Using Rscript 
+ For Window users  : Double click the icon of tReasure on the desktop.

+ For Linux or MacOS users : Click on the tReasure and type on command as below.
<pre>
<code>
Rscript tReasure.R
</code>
</pre>

