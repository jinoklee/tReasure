

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

## Installation       
tReasure is a package for the R computing environment and it is assumed that you have already installed R according to your operating system. See the R project at (http://www.r-project.org). tReasure requires the gwidget2 graphical library to run and a few additional packages for the analysis of RNA-seq. 

<br/>   

### **Method 1. Install tReasure from GitHub**  : Open Rstudio or R and type as below   
   
    install.packages("devtools")
    library("devtools")
    devtools::install_github("jinoklee/tReasure")
    
   ***Note*** If the descriptioin shows as below during the installation, choose 1 
    
    
    Select a GUI toolkit
        1: gWidgets2RGK2
        2: gWidtets2ticltk
    

### **Method 2. Install tReasure from source**   
    In case, addtional R tools (http://www.r-project.org) must be installed.
    The difference from installing with GitHub is that it is installed as a standalone tools

   **STEP 1.** Download and unzip a file
   + Windows: [tReasrue_win.zip]
   + Linux/Mac OS: [tReasure_src.zip]      

   **STEP 2.** Double-click or type on command window an install-file  
   + Windows: install_win.bat
   + Mac OS: install.sh.command
   + Linux: open commnad window and type as below       
        ~~~   
        sh install_src_v1.sh
        ~~~   

Installation may take several minutes. 
+ It automatically installs tReasure package.
+ It creates a folder named tReasure_v1 on Documents. 
+ It creates **shortcut for tReasure** on Desktop (only Window).

<br/>

## Start   
### **Method 1. In case downloading tReasure from GitHub**  : Open Rstuido or R and type as below  
  
    
    library(tReasure)
    tReasure()
       

### **Method 2.  In case downloading tReasrue from source**  
   + Windows: Double-click the icon of tReasure on Desktop.  

   + MacOS: Ddouble-Click **run_tReasure.sh.command**   

   + Linux: type as belows         
       ~~~
        chmod 777 run_tReasure.sh.command
        sh run_tReasure.sh.command
       ~~~
    
<br/>   


### Sample data  
* Human breast cancer [Download](https://www.dropbox.com/sh/phkerfxxq3jmgo9/AAC3sR1rWWo5DsTZAD3_VUANa?dl=0)  




