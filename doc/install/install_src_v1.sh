#/bin/bash

mkdir $HOME/Documents/tReasure
cp -r $HOME/Downloads/tReasure_src.zip $HOME/Documents/tReasure

unzip $HOME/Documents/tReasure/tReasure_src.zip
Rscript $HOME/Documents/tReasure/pkg_lib.R
Rscript $HOME/Documents/tReasure/loadR_v1.R
