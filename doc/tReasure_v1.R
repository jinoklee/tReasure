options(guiToolkit="RGtk2")
library(gWidgets2)
library(gWidgets2RGtk2)
library(cairoDevice)
library(QuasR)
library(edgeR)
library(ggplot2)
library(grid)
library(gridExtra)
library(plotrix)
library(tidyverse)
library(stringr)
library(Rsamtools)
library(seqinr)
library(ShortRead)
library(dplyr)
library(DESeq2)
library(future)
library(statmod)

# Sys.setlocale('LC_ALL','C')
#-------------------------------------------------------------------------------------
# setting the PATH for TEST : system.file('', package="tReasure")
#......................................................................................#
intro <- system.file("extdata", "intro.png", package = "tReasure", mustWork = TRUE)
cl_name <-  read.table(system.file("extdata", "class_name.txt", package = "tReasure",mustWork = TRUE), sep = "\t", fill = T,header = T, as.is = T)
#-------------------------------------------------------------------------------------
# function
#......................................................................................#
source(system.file("extdata", "func/stopFuture.R", package = "tReasure", mustWork = TRUE))
source(system.file("extdata", "func/mk_sample.R", package = "tReasure", mustWork = TRUE))
source(system.file("extdata", "func/sel_sample.R", package = "tReasure", mustWork = TRUE))
source(system.file("extdata", "func/anl_trim.R", package = "tReasure", mustWork = TRUE))
source(system.file("extdata", "func/anl_align.R", package = "tReasure", mustWork = TRUE))
source(system.file("extdata", "func/anl_rm.R", package = "tReasure", mustWork = TRUE))
source(system.file("extdata", "func/anl_rc.R", package = "tReasure", mustWork = TRUE))
source(system.file("extdata", "func/anl_filter.R", package = "tReasure", mustWork = TRUE))
source(system.file("extdata", "func/deseq.R", package = "tReasure", mustWork = TRUE))
source(system.file("extdata", "func/edger.R", package = "tReasure", mustWork = TRUE))
source(system.file("extdata", "func/anl_deseq.R", package = "tReasure", mustWork = TRUE))
source(system.file("extdata", "func/anl_edger.R", package = "tReasure", mustWork = TRUE))
source(system.file("extdata", "func/pyramid.R", package = "tReasure", mustWork = TRUE))
source(system.file("extdata", "func/volcano.R", package = "tReasure", mustWork = TRUE))
source(system.file("extdata", "func/barplot.R", package = "tReasure", mustWork = TRUE))

#-------------------------------------------------------------------------------------
#  tReasure
#......................................................................................#
window <- gwindow("tReasure")
addHandlerUnrealize(window, handler = function(h,...) {
  val <- gconfirm("Are you sure you want to exit tReausre?", parent=h$obj)
  if(as.logical(val)){
    if(dir.exists("apid")){
      load("apid")
      stopFuture(apid)}
    if(dir.exists("bpid")){
      load("bpid")
      stopFuture(bpid)}
    dispose(window)
    gtkMainQuit()
  }else{
    return(TRUE)
  }
})


mother <- ggroup(container = window, horizontal = FALSE)
size(mother) <- c(1000,740)

# Sub window (Notebook)------------
notebook <- gnotebook(container = mother)
gr0 <- ggroup(container = notebook, label ="  Introduction  ")
gr1 <- ggroup(container = notebook, label ="  Uploading Samples  ", horizontal = FALSE)
gr2 <- ggroup(container = notebook, label ="  Quality Control  ",horizontal = FALSE)
gr3 <- ggroup(container = notebook, label ="  Alignment & Counting  ",horizontal = FALSE)
gr4 <- ggroup(container = notebook, label ="  Filtering   ",horizontal = FALSE)
#gr5 <- ggroup(container = notebook, label ="  Exploration   ",horizontal = FALSE)
gr6 <- ggroup(container = notebook, label ="  DEtRNA Detection ",horizontal = FALSE)
gr7 <- ggroup(container = notebook, label ="  Visualization ",horizontal = FALSE)
child <- gframe("  Progress  ", container = mother, horizontal = FALSE)
size(child) <- c(1000,100)
st <- gtext("  ",container = child)
size(st) <- c(1000,100)

#-------------------------------------------------------------------------------------
#  gr0. Introduction
#......................................................................................
ggr1 <- ggroup(container = gr0, horizontal = TRUE, fill = "both", expand = TRUE)
gimage(intro, container = ggr1)
#-------------------------------------------------------------------------------------
#  gr1. Uploading Samples
#......................................................................................
addSpace(gr1, 20)
ggr1 <- ggroup(container = gr1,  spacing = 10) ; size(ggr1) <- c(896,532)
addSpace(ggr1, 10)
paned.1 <- gpanedgroup(container = ggr1, horizontal = FALSE, spacing = 10)
tmp.1 <- gframe(" Create the sample list ", container = paned.1, horizontal = FALSE, spacing = 10 ) ; size(tmp.1) <- c(250,284)
addSpace(tmp.1, 10)
glabel(" Select a directory of FASTQ files : ", container = tmp.1, anchor = c(-1,0))
addSpace(tmp.1, 10)
fq_dir <- gfilebrowse(text = " ", quote = FALSE, type = "selectdir", container = tmp.1)
addSpace(tmp.1, 20)
make_button <- gbutton("RUN", container = tmp.1, handler = mk_sample)

tmp.11 <- gframe(" [*OPTION] Upload the sample list ", container = paned.1, horizontal = FALSE, spacing = 10) ; size(tmp.11) <- c(250,200)
addSpace(tmp.11, 10)
glabel(" Select a file of \n the pre-made sample list : ", container = tmp.11, anchor = c(-1,0))
addSpace(tmp.11, 20)

sel_button <- gfilebrowse(text=" ", quote = FALSE, type="open", container = tmp.11,
                          filter=list("*.txt" = list(patterns = c("*.txt")), "*.*" = list(patterns = c("*"))),
                          handler = sel_sample)

tbl <- gdf(data.frame(FileName=I(list("Path/samplename.fastq", NA, NA)),
                      SampleName=I(list("samplename", NA, NA)),
                      Group=I(list("control", NA, NA)),
                      Batch=I(list("batch I", NA, NA))), container = ggr1, multiple = FALSE)


addHandlerChanged(fq_dir, handler = function(h,...){
  setwd(svalue(fq_dir))
  dir <- getwd()
  if(!dir.exists(paste(dir, "/pre", sep = ""))){
    dir.create(paste(dir, "/pre", sep = ""), recursive = TRUE)}

  if(!dir.exists(paste(dir, "/post", sep = ""))){
    dir.create(paste(dir, "/post", sep = ""), recursive = TRUE)}

  if(!dir.exists(paste(dir, "/stat/plot", sep = ""))){
    dir.create(paste(dir, "/stat/plot", sep = ""), recursive = TRUE)}
  if(!dir.exists(paste(dir, "/rc", sep = ""))){
    dir.create(paste(dir, "/rc", sep = ""), recursive = TRUE)}
})


#-------------------------------------------------------------------------------------
#  gr2. Quality Control
#......................................................................................#
addSpace(gr2, 20)
ggr21 <- ggroup(container = gr2,  spacing = 10); size(ggr21) <- c(896,532)
addSpace(ggr21, 10)
paned.2 <- gpanedgroup(container = ggr21, horizontal = FALSE, spacing = 10)
tmp.2 <- gframe("  Settings  ", container = paned.2, horizontal = FALSE, spacing = 10); size(tmp.2) <- c(250,532)
addSpace(tmp.2, 10)

glabel(" Information of adapter : ", container = tmp.2, anchor = c(-1,0))
addSpace(tmp.2, 10)

adapt <- gradio(c("Illumina smallRNA 3' adapter", "Illumina universal adapter","SOLiD adapter","No adapter"), container = tmp.2)
addSpace(tmp.2, 10)

glabel(" Minimum quality : ", container = tmp.2, anchor = c(-1,0))
addSpace(tmp.2, 10)
q <- gradio(c("25", "30"), container = tmp.2)
addSpace(tmp.2, 10)

glabel(" Minimum lenght : ", container = tmp.2, anchor = c(-1,0))
addSpace(tmp.2, 10)
min <- gradio(c("10"), container = tmp.2)
addSpace(tmp.2, 20)

trim_button <- gbutton("RUN", container = tmp.2, expand = FALSE, handler = anl_trim)


#-------------------------------------------------------------------------------------
#  gr3. Alignment & Counting
#......................................................................................
addSpace(gr3, 20)
ggr31 <- ggroup(container = gr3,  spacing = 10); size(ggr31) <- c(896,532)
addSpace(ggr31, 10)
paned.3 <- gpanedgroup(container = ggr31, horizontal = FALSE, spacing = 10)
tmp.3 <- gframe("  Settings  ", container = paned.3 , horizontal = FALSE, spacing = 10); size(tmp.3) <- c(250,260)
pangr <- ggroup(container = paned.3, horizontal =FALSE)
addSpace(tmp.3, 10)

RC_view <- gnotebook(container = ggr31 ); size(RC_view) <- c(700, 500)
RC_trna <- ggroup(container = RC_view, horizontal = FALSE, label = " tRNA level ")
RC_codon <- ggroup(container = RC_view, horizontal = FALSE, label = " Isodecoder level")
RC_aa <- ggroup(container = RC_view, horizontal = FALSE, label = "Isoacceptor level")

glabel(" Genome assembly : ", container = tmp.3, anchor = c(-1,0))
addSpace(tmp.3, 10)


### child_window--------------------
          child_window <- gwindow("Search genome assembly", parent = window, width = 600, height = 300, visible = FALSE )
          ggroup(container = child_window, spacing = 10)
          gr_child <- ggroup(container = child_window, spacing = 10)

          addSpace(gr_child, 20)
          gr_frame <- gframe(" ", container = gr_child, horizontal = FALSE, spacing = 10)
          size(gr_frame) <- c(550, 250)
          addSpace(gr_child, 20)

          addSpace(gr_frame, 10)
          c1 <- glabel("  Domain : ", container = gr_frame, anchor = c(-1,0))
          ref_P1 <- gcombobox(c(" ", "Eukaryota" ,"Archaea", "Bacteria"),container = gr_frame) # container = tmp.3
          addSpace(gr_frame, 10)

          c2 <- glabel("  Next 1 : ", container = gr_frame, anchor = c(-1,0))
          ref_P2 <- (v1names <- gcombobox("  ",container = gr_frame)) # container = tmp.3
          addSpace(gr_frame, 10)

          c3 <- glabel("  Next 2 : ", container = gr_frame, anchor = c(-1,0))
          ref_P3 <- (v2names <- gcombobox("  ",container = gr_frame)) # container = tmp.3
          addSpace(gr_frame, 10)

          c4 <- glabel("  Species : ", container = gr_frame, anchor = c(-1,0))
          ref_P4 <- (v3names <- gcombobox(" ",container = gr_frame)) # container = tmp.3

          p2Nms <- function(d, envir=.GlobalEnv)  unique(cl_name$P2[grep(svalue(ref_P1), cl_name$P1)], envir=envir)
          p3Nms <- function(d, envir=.GlobalEnv)  unique(cl_name$P3[grep(svalue(ref_P2), cl_name$P2)], envir=envir)
          p4Nms <- function(d, envir=.GlobalEnv)  unique(cl_name$P4[grep(svalue(ref_P3), cl_name$P3)], envir=envir)

          addHandlerChanged(ref_P1, handler = function(h,...){
            ref_P2[] <- p2Nms(svalue(ref_P2))
          })

          addHandlerChanged(ref_P2, handler = function(h,...){
            ref_P3[] <- p3Nms(svalue(ref_P3))
          })

          addHandlerChanged(ref_P3, handler = function(h,...){
            ref_P4[] <- p4Nms(svalue(ref_P4))
          })

          addSpace(gr_frame, 20)
          g_run <- gbutton("OK", container = gr_frame, handler = function(h, ...){
              svalue(gsel_button) <<- svalue(ref_P4)
              visible(child_window) <- FALSE
              insert(st, "Click the RUN", do.newline = TRUE)
          })
          addSpace(gr_child, 20)

gsel_button <- gbutton("Search genome assembly", container = tmp.3, handler = function(h, ...){
  visible(child_window) <- TRUE
})
addSpace(tmp.3, 20)

anl_button <- gbutton("RUN", container=tmp.3, handler = anl_align)

addSpace(pangr, 10)
# tmp.31 <- gframe(" [*OPTION] Working directory ", container = pangr , horizontal = FALSE, spacing = 10); size(tmp.31) <- c(250,70)
# addSpace(tmp.31, 10)
#
# glabel(" Select the directory of the FASTQ files : ", container = tmp.31, anchor = c(-1,0))
# addSpace(tmp.31, 10)

tmp.32 <- gframe(" [*OPTION] Load the readcounts files ", container = pangr , horizontal = FALSE, spacing = 10); size(tmp.32) <- c(250,120)
addSpace(tmp.32, 10)

glabel(" Select the file of readcounts : ", container = tmp.32, anchor = c(-1,0))
addSpace(tmp.32, 10)

# selfq_button <- gfilebrowse(text = "", quote = FALSE, type = "selectdir", container = tmp.31,
#                             filter=list("*.fa" = list(patterns = c("*.fa")), "*.*" = list(patterns = c("*"))))

seltrn_button <- gfilebrowse(text = "tRNA level", quote = FALSE, type = "open", container = tmp.32,
                             filter=list("*.txt" = list(patterns = c("*.txt")), "*.*" = list(patterns = c("*"))))

selco_button <- gfilebrowse(text = "Isodecoder level", quote = FALSE, type = "open", container = tmp.32,
                            filter=list("*.txt" = list(patterns = c("*.txt")), "*.*" = list(patterns = c("*"))))

selaa_button <- gfilebrowse(text = "Isoacceptor level", quote = FALSE, type = "open", container = tmp.32,
                            filter=list("*.txt" = list(patterns = c("*.txt")), "*.*" = list(patterns = c("*"))))


# addHandlerChanged(selfq_button, handler = function(h,...){
#   setwd(svalue(selfq_button))
#   dir <- getwd()
#   if(!dir.exists(paste(dir, "/pre", sep = ""))){
#     dir.create(paste(dir, "/pre", sep = ""), recursive = TRUE)}
#   if(!dir.exists(paste(dir, "/stat/plot", sep = ""))){
#     dir.create(paste(dir, "/stat/plot", sep = ""), recursive = TRUE)}
# })

addHandlerChanged(seltrn_button, handler = function(h,...){
  tcount <- read.delim(svalue(seltrn_button))
  gtable(tcount, container=RC_trna)
})

addHandlerChanged(selco_button, handler = function(h, ...){
  ccount = read.delim(svalue(selco_button))
  gtable(ccount, container=RC_codon)
})

addHandlerChanged(selaa_button,handler = function(h, ...){
  acount = read.delim(svalue(selaa_button))
  gtable(acount, container=RC_aa)
})


#-------------------------------------------------------------------------------------
#  gr4. Filtering
#......................................................................................
addSpace(gr4, 20)
ggr41 <- ggroup(container = gr4,  spacing = 10); size(ggr41) <- c(896,532)
addSpace(ggr41, 10)
paned.4 <- gpanedgroup(container = ggr41, horizontal = FALSE, spacing = 10)
tmp.4 <- gframe("  Settings   ", container = paned.4, horizontal = FALSE, spacing = 10)
size(tmp.4) <- c(250,284)
addSpace(tmp.4, 20)

prelyt <- gformlayout(container = tmp.4, spacing = 1.5)
cfilv <- gcombobox(seq(0,10,by=1),label = "cpm value >=  ", container = prelyt)
sfil <- gcombobox(seq(0,100, by=10), label = "sample (%) >=  ", container = prelyt)

addSpace(tmp.4, 20)
Pre_button <- gbutton("RUN", container = tmp.4, handler = anl_filter)

tmp.42 <- gframe(" [*OPTION] Working directory ", container = paned.4, horizontal = FALSE, spacing = 10)
size(tmp.42) <- c(250,200)
addSpace(tmp.42, 10)
glabel(" Select the directory of the readcounts files :", container = tmp.42, anchor = c(-1,0))
addSpace(tmp.42, 10)

selrc_button <- gfilebrowse(text = "", quote = FALSE, type = "selectdir", container = tmp.42,
                            filter=list("*.txt" = list(patterns = c("*.txt")), "*.*" = list(patterns = c("*"))))

f_kep <- gtable(data.frame(No.reads=c("total", "keep", "remove")),container = ggr41, multiple = FALSE)

#-------------------------------------------------------------------------------------
#  gr5. Exploration
#......................................................................................
# addSpace(gr5, 20)
# ggr51 <- ggroup(container = gr5, spacing = 10, fill=TRUE); #size(ggr51) <- c(996,550)
# addSpace(ggr51, 10)
# mds_view <- gnotebook(container = ggr51, spacing = 10, fill=TRUE, expand=TRUE); #size(mds_view) <- c(950, 550)
# mds_plot <- ggroup(container = mds_view, horizontal = FALSE, label = " MDS plot ")
# addSpace(ggr51, 10)
# addSpace(gr5, 20)


#-------------------------------------------------------------------------------------
#  gr6. DEtRNA Detection
#......................................................................................
addSpace(gr6, 20)
ggr61 <- ggroup(container = gr6,  spacing = 10); size(ggr61) <- c(896,532)
addSpace(ggr61, 10)
paned.6 <- gpanedgroup(container = ggr61, horizontal = FALSE, spacing = 10)
addSpace(ggr61, 10)

ggr62 <- gnotebook(container=paned.6, tab.pos  = 3) ; size(ggr62) <- c(250, 200)
stat2 <- ggroup(container = ggr62, label = " DEseq2 ")
stat1 <- ggroup(container = ggr62, label = " EdgeR ")

tmp.61 <- gframe("  Statics  ", container = stat1, horizontal = FALSE, spacing = 10); size(tmp.61) <- c(235,150)
addSpace(tmp.61, 10)
glabel(" Stat.Method : ", container = tmp.61, anchor = c(-1,0))
addSpace(tmp.61, 10)
method_edgeR <- gradio(c("Exact test", "Likelihood F-test", "Quasi-likelihood F-test"), container = tmp.61)
addSpace(tmp.61, 10)

tmp.62 <- gframe("  Statics  ", container = stat2, horizontal = FALSE, spacing = 10); size(tmp.62) <- c(235,150)
addSpace(tmp.62, 10)
glabel(" Stat.Method : ", container = tmp.62, anchor = c(-1,0))
addSpace(tmp.62, 10)
method_deseq <- gradio("Walt test", container = tmp.62)
addSpace(tmp.62, 10)


statedgeR_button <- gbutton ("RUN EdgeR", container = tmp.61, handler = anl_edger)
statdeseq_button <- gbutton("RUN DESeq2", container = tmp.62, handler = anl_deseq)


widget_list <- list()
ggr.63 <- gframe("  Threshold  ", container = paned.6, horizontal = FALSE, spacing = 10)
addSpace(ggr.63, 10)
fdr_label<- glabel(" Adjust.Method : ", container = ggr.63, anchor = c(-1,0))
addSpace(ggr.63, 10)
widget_list$fdr_s <- gcombobox(c("FDR", "Bonferroni","Benjamini-Hochberg"  ), container = ggr.63)
addSpace(ggr.63, 10)

statlyt <- gformlayout(container = ggr.63, spacing = 1.5)
widget_list$pval <- gcombobox(c(0,0.001, 0.05, 0.01),label = "Adj P.value  < : ", container = statlyt)
widget_list$FC <- gcombobox(c(0,1,1.5,2),label = "Fold change  > : ", container = statlyt)

addSpace(ggr.63, 20)

DE_view <- gnotebook(container = ggr61, expand=TRUE ); #size(DE_view) <- c(700, 532)
DE_trna <- ggroup(container = DE_view, horizontal = FALSE, label=" tRNA level ")
DE_codon <- ggroup(container = DE_view, horizontal = FALSE, label=" Isodecoder level")
DE_aa <- ggroup(container=DE_view, horizontal = FALSE, label= "Isoacceptor level")

stat_trna <- gtable(data.frame(), container = DE_trna)
stat_codon <- gtable(data.frame(), container = DE_codon)
stat_aa <- gtable(data.frame(), container = DE_aa)


# Threshold handler  --------------------
update_trna_frame <- function(...){
  data <- read.delim("./stat/stat_trna_list.txt")
  if(svalue(widget_list$fdr_s) == "Benjamini-Hochberg"){
    stat_trna[] <- filter(data, data$Benjamini < svalue(widget_list$pval) & abs(data$logFC) > svalue(widget_list$FC))
  }else if(svalue(widget_list$fdr_s) =="Bonferroni"){
    stat_trna[] <- filter(data, data$Bonferroni < svalue(widget_list$pval) & abs(data$logFC) > svalue(widget_list$FC))
  }else{stat_trna[] <- filter(data, data$FDR < svalue(widget_list$pval) & abs(data$logFC) > svalue(widget_list$FC))}
}

update_codon_frame <- function(...){
  data <- read.delim("./stat/stat_isodecoder_list.txt")
  if(svalue(widget_list$fdr_s) == "Benjamini-Hochberg"){
    stat_codon[] <- filter(data, data$Benjamini < svalue(widget_list$pval) & abs(data$logFC) > svalue(widget_list$FC))
  }else if(svalue(widget_list$fdr_s) =="Bonferroni"){
    stat_codon[] <- filter(data, data$Bonferroni < svalue(widget_list$pval) & abs(data$logFC) > svalue(widget_list$FC))
  }else{stat_codon[] <- filter(data, data$FDR < svalue(widget_list$pval) & abs(data$logFC) > svalue(widget_list$FC))}
}

update_aa_frame <- function(...){
  data <- read.delim("./stat/stat_isoacceptor_list.txt")
  if(svalue(widget_list$fdr_s) == "Benjamini-Hochberg"){
    stat_aa[] <- filter(data, data$Benjamini < svalue(widget_list$pval) & abs(data$logFC) > svalue(widget_list$FC))
  }else if(svalue(widget_list$fdr_s) =="Bonferroni"){
    stat_aa[] <- filter(data, data$Bonferroni < svalue(widget_list$pval) & abs(data$logFC) > svalue(widget_list$FC))
  }else{stat_aa[] <- filter(data, data$FDR < svalue(widget_list$pval) & abs(data$logFC) > svalue(widget_list$FC))}
}

callback <- function(h, ...){
  update_trna_frame()
  update_codon_frame()
  update_aa_frame()
}

sapply(widget_list, addHandlerChanged, handler = callback)
DE_button <- gbutton(" Plot ", container = ggr.63)

# Plot save  handler --------------------
addHandlerClicked(DE_button, handler = function(h, ...){
  DEtRNA <- stat_trna[]
  DEcodon <- stat_codon[]
  DEaa <- stat_aa[]

  write.table(DEtRNA, "./stat/DEtrna_list.txt", sep = "\t", quote = FALSE)
  write.table(DEcodon, "./stat/DEisodecoder_list.txt", sep = "\t", quote = FALSE)
  write.table(DEaa, "./stat/DEisoacceptor_list.txt", sep = "\t", quote = FALSE)

  volcano(650, 460, 80)
  barplot(850, 460,80)
  pyramid(650, 460,80)

  insert(st, " ", do.newline = TRUE)
  insert(st,"Done : The plots of DEtRNA Detection. Click! Next tab of Visualization.", do.newline = TRUE )
  insert(st, ".", do.newline = TRUE)
})


#-------------------------------------------------------------------------------------
#  gr7. Visualiztion
#......................................................................................
addSpace(gr7, 20)
ggr71 <- ggroup(container = gr7,  spacing = 10); size(ggr71) <- c(996,550)
addSpace(ggr71, 10)
Plot_view <- gnotebook(container = ggr71,spacing = 10, expand=TRUE ); #size(Plot_view) <- c(950, 550)
devs <- lapply(c("MDS plot","Plot_trnas", "Plot_isodecoders", "Plot_isoacceptors"), function(i) ggraphics(container = Plot_view, visible = TRUE, label = as.character(i)))

# Plot viewe handler--------------------
addSpace(ggr71, 10)
addHandlerChanged(Plot_view, handler = function(h,...){
  gg <- h$obj[h$page.no]
  visible(gg) <- TRUE
  #print(p())
  if(h$page.no == "1"){
    load("./stat/plot/MDS_plot.RData")
    print(p())
  }else if(h$page.no == "2" ){
    load("./stat/plot/Volcano_Plot.RData")
    print(p())
  }else if(h$page.no == "3"){
    load("./stat/plot/Bar_Plot.RData")
    print(p())
  }else{
    load("./stat/plot/Pyramid_Plot.RData")
    print(p())
  }
})

gtkMain()



