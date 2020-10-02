
<br/>
<br/>
<br/>
<br/>

# tReasure
***
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
***
tReasure is a package for the R computing environment and it is assumed that you have already installed R according to your operating system. See the R project at (http://www.r-project.org). tReasure requires the gwidget2 graphical library to run and a few additional packages for the analysis of RNA-seq. 

<br/>
<br/>

Install source : https://github.com/jinoklee/tReasure_install.git

***

### **STEP 1.** Download and unzip a file
   > ##### Window : [tReasrue_v1_win.zip](https://www.dropbox.com/s/gnq144mputz4fdm/easy_int_win_v1.0.zip?dl=0)
   > ##### Linux/Mac OS : [tReasure_v1_mac.zip]
   
<br/>

### **STEP 2.** Double-click an install-file 
   > ##### Window : int_win_v1.0.bat
   > ##### Linux/Mac OS : int_mac_v1.0.sh.command
Installation may take several minutes. 
+ It automatically installs **tReasure package**.
+ It creates a folder named **tReasure_v1** on Documents. 
+ It creates **shortcut for tReasure** on Desktop (only Window).

<br/>

### **STEP 3.** Download the genome indices as needed
   > ##### Donwload : [a bundle of genome indics files](https://www.dropbox.com/sh/1aikvdszjlvncic/AADzL8G55ayI3lRfzZ6LYjvPa?dl=0)
+ hg38 genome indices are saved as default.
+ If you need another genome indices files, download and move the files inside tReasure packages folder (~Documents/R/win-library/tReasure/WholeGenomeFasta).


<br/>
<br/>

## Start
***
#### Using Rscript 
+ For Window users  : Double-click the icon of tReasure on the desktop.

+ For Linux or MacOS users : Double-Click the files named **run.sh.command** on Document/tReasure_v1


#### Using R or Rstudio

Open R or Rstudio and type on console as below.
<pre>
<code>
library(tReasure)
tReasure()
</code>
</pre>



