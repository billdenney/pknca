rm(list=ls())

## options(repos="http://cran.r-project.org/")
## install.packages("devtools")
## install.packages("testthat")
## install.packages("roxygen2")

## library(nlme)
library(devtools)

load_all(".")
test(".")
document(".")
check(".")

plot(Theoph)
theoph <- as.data.frame(Theoph)

myconc <- PKNCAconc(as.data.frame(Theoph), conc~Time|Subject)
mydose <- PKNCAdose(data.frame(Subject=unique(Theoph$Subject),
                       Time=0),
                    ~Time|Subject)

tmpconc <- splitBy(parseFormula(myconc)$groupFormula, myconc$data)
tmpdose <- splitBy(parseFormula(mydose)$groupFormula, mydose$data)
tmpmerge2 <- merge(conc=tmpconc, dose=tmpdose)

mydat <- 
  PKNCAdata(as.data.frame(Theoph),
            conc~Time|Subject,
            data.frame(Subject=unique(Theoph$Subject),
                       Time=0),
            ~Time|Subject)

q(save="no")
