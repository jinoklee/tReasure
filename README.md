[header](https://github.com/jinoklee/tReasure/blob/master/docs/header.png?raw=true)

# tReasure
__________
tReasure (tRna Expression Analysis Software Utilizing R for Easy use) is user-friendly tool for the analysis of deep-sequencing experiments for small RNAs using R packages. The tool requires a simple input file containing a list of unique reads. Using these data, tReasure (i) detects all known tRNA sequences annotated in GtRNAdb, (ii) predicts the specific expressed tRNA compare to control.

## Installation guide
tReasure is a package for the R computing environment and it is assumed that you have already installed R. See the R project at (http://www.r-project.org). tReasure requires the gwidget2 graphical library to run and a few additional packages for the analysis of RNA-seq. 
  
  
tReasure currently automates the following tasks:
1) preprocessing of the raw sequence reads
2) pre-processing 
3) alignment and quantification
4) statistical analysis for differentially expressed tRNAs gene
5) visualization
