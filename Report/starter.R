
setwd("E:/Post-doc_data/Action_foreperiod")

my.knit = knitr::knit("action_foreperiod_paper.Rnw")
## document.tex is the latex file that will be compiled by the two following command:

system(paste0("pdflatex ", "action_foreperiod_paper.tex")) 
system(paste0("pdflatex ", "action_foreperiod_paper.tex")) 
