

<img src = "https://github.com/jinoklee/tReasure/blob/master/inst/extdata/tresure.png" width="100" height="100" />

# tReasure
***
**tReasure (tRNA Expression Analysis Software Utilizing R for Easy use)** is user-friendly tool for the tRNA expression analysis of deep-sequencing experiments for small RNAs using R packages. 

    tReasure currently implements the following tasks:
    1) Making sample list for analysis
    2) Pre-processing of trimming adapter and filtering reads
    3) Alignment and quantification
    4) Statistical analysis for differentially expressed tRNAs gene
    5) Visualization 
<br/>

## 🌱 Installation       
tReasure is a package for the R computing environment and it is assumed that you have already installed R according to your operating system. See the R project at (http://www.r-project.org). tReasure requires the gwidget2 graphical library to run and a few additional packages for the analysis of RNA-seq. 


### 1. Preliminaries  

  1) For window users :  
     User need to install [Rtools](https://cran.r-project.org/bin/windows/Rtools/rtools40.html) before installing tReasure.
  
  5) For linux users :  
     Users need to install **libcurl4-openssl-dev** and **libssl-dev** for devtools.  
     - Example for user Ubuntu
     
     ```
     sudo apt-get install libcurl4-openssl-dev 
     sudo apt-get install libssl-dev
     ```
     
     Also, user need to install **libgtk2.0-dev** and **libxml2-dev** for tReasure.  
     - Example for user Ubuntu
     
     ```
     sudo apt-get update -y
     sudo apt-get install -y libxml2-dev
     sudo apt-get install -y libgtk2.0-dev
     ```
     
  2) For Mac OS users :  
     Users need to install **opnssl** for devtools.   
     - Example using brew
     ```
     brew update
     brew install openssl
     ```
     Users need to install [XQuartz](https://www.xquartz.org), **gtk2** and **cairo** before installing tReasure.  
     - Example using brew 
     ```
     brew install cairo  
     brew install gtk+
     ```
 
     
### 2. Installation of tReasure : Open Rstudio or R and type as below.
    It takes a few minutes to install for the first time.  
 
   
     install.packages("devtools")
     library("devtools")
    
     devtools::install_github("jinoklee/tReasure",force = TRUE)
     library("tReasure")
     install.tReasure()
   
    
   ***Note*** If the descriptioin shows as below during the installation, choose **"Install GTK+"** 
    
    
    Need GTK+? (Restart Required)
        Install GTK+
        Do not install GTK+
        
    
    
   ***Note*** If the installation was successful, the tReasure window appears. If the window dose not appear, restart R and type as below.  
   
    tReasure::tReasure()
   

## 🌱 Running tReasure  
### Open Rstuido or R and type as below

     library("tReasure")
     tReasure::tReasure()


## 🔔 User Manual
* [Download](doc/tReasure-Manual-220109.pdf)

## 🔔 Detailed Mapping Methods
* [Click](doc/Detailed-Mapping-Methods.pdf)
## 🔔 Sample Data  
* Small RNA-seq dataset from human breast tissues (part of GSE68085)
    - [Download](https://www.dropbox.com/sh/phkerfxxq3jmgo9/AAC3sR1rWWo5DsTZAD3_VUANa?dl=0)
    - Control (Normal): SRR1982473, SRR1982474, SRR1982475
    - Test (Cancer): SRR1982580, SRR1982581, SRR1982582

## 🔔 Tutorial Videos
 * Installation and setup: [Window](https://www.dropbox.com/s/ssxux5ad7jwvxk7/win_install.mkv?dl=0) 
 * Full version of analysis: [Window](https://www.dropbox.com/s/vinwwdl1umw74l7/tReasure.Analysis.mp4?dl=0) [Mac](https://www.dropbox.com/s/ylt89pxcw9bf4eq/tReasure.Analysis.Mac.mov?dl=0)
 * Results folders: [Click](https://www.dropbox.com/s/3kwx3i45sllnnhu/tReasure.Outputs.mp4?dl=0)

