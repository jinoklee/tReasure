


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
   
<br/>
<br/>



## Installation guide

**tReasure** is a package for the R computing environment and it is assumed that you have already installed R according to your operating system. See the R project at (http://www.r-project.org). tReasure requires the gwidget2 graphical library to run and a few additional packages for the analysis of RNA-seq. 


<br/>

### _**STEP1. Install the devtools package**_

To install tReasure package, start by installing the devtools package. Open R or Rstudio and type on console as below.
<pre>
<code>
install.packages(“devtools”)
library(devtools)

install_github(“jinoklee/tReasure”) 
library(tReasure) 
</code>
</pre>

### _**[Other methods]**_
You can also install tReasure packages from local source. Download <p style= "color:red">tReasure_1.0.0.tar.gz</p> and saved it on Download folders.Open R or Rstudio and type on console as below.
<pre>
<code>
install.packages(“~/Download/tReasure_1.0.0.tar.gz”, repo=NULL, type = “source”)
library(tReasure)
</code>
</pre>

<br/>

### _**STEP2. Download files contained genome index and sample raw data**_

Before the start tReasure packages, you should download the genome index files what you want ( hg38, hg19, mm10 or all). 
+ Download <span style= "color:red">a bundle of genome index files</span> and move the files inside tReasure package folder (~Documents/R/win-library/tReasure/WholeGenomeFasta)._Caution! Before moving, you should make a folder named “WholeGenomeFasta”._
+ You can also download sample raw files for a test. 

<br/>

### _**STEP3. Download tReasure Rscript for analysis using Rscript execution**_

For a standalone tReasure for window user, download script’s file (tReasure.R). Download a bundle of files contained Rscript (tReasure.R) and a batch file named shortcut_install.bat.Click on the shortcut_install.bat, it creates a shortcut icon of tReasure on the desktop.

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

