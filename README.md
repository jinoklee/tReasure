
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


  #### Workflow
   ![Flow](https://github.com/jinoklee/tReasure/blob/master/inst/extdata/flow.png?raw=true)
   
<br/>
<br/>
<br/>


## Installation guide

tReasure is a package for the R computing environment and it is assumed that you have already installed R according to your operating system. See the R project at (http://www.r-project.org). tReasure requires the gwidget2 graphical library to run and a few additional packages for the analysis of RNA-seq. 

<br/>


***

### **STEP 1.** Download and unzip a file
   > ##### Win : [easy_int_tReasrue_win_v1.0.zip](https://www.dropbox.com/s/gnq144mputz4fdm/easy_int_win_v1.0.zip?dl=0)
   > ##### Linux/Mac OS : 
### **STEP 2.** Double-click a install file 
   > ##### Win : easy_int_win_v1.0.bat
   > ##### Linux/Mac OS : 
Installation may take several munutes. 
+ It automatically installs **tReasure package**.
+ It creates **shortcut for tReasure** on Desktop. 
+ It creates a folder named **tReasure_v1** on Document. 
### **STEP 3.** Download the genome index as needed
   > ##### Donwload : [a bundle of genome index files](https://www.dropbox.com/sh/1aikvdszjlvncic/AADzL8G55ayI3lRfzZ6LYjvPa?dl=0)
+ hg38 genome index is saved as default.
+ If you need another genome index files, download and move the files inside tReasure package folder (~Documents/R/win-library/tReasure/WholeGenomeFasta).


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

