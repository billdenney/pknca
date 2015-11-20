rm(list=ls())

## options(repos="http://cran.r-project.org/")
## install.packages("devtools")
## install.packages("testthat")
## install.packages("roxygen2")
## install.packages("plyr")
## install.packages("covr")
## install.packages("knitr")
## install.packages("Rcpp")
## install.packages("shiny")
## install.packages("doBy")

## library(nlme)
library(devtools)
library(covr)
library(testthat)
library(shiny)

load_all(".")
results <- test()
document(".")
check(".")

mycoverage <- package_coverage()
mycoverage
shine(mycoverage)

plot(Theoph)
theoph <- as.data.frame(Theoph)

## myconc <- PKNCAconc(as.data.frame(Theoph), conc~Time|Subject)
## mydose <- PKNCAdose(data.frame(Subject=unique(Theoph$Subject),
##                        Time=0),
##                     ~Time|Subject)
##
## tmpconc <- splitBy(parseFormula(myconc)$groupFormula, myconc$data)
## tmpdose <- splitBy(parseFormula(mydose)$groupFormula, mydose$data)
## tmpmerge2 <- merge(conc=tmpconc, dose=tmpdose)

myconc <-
  PKNCAconc(
    data=as.data.frame(Theoph),
    formula=conc~Time|Subject,
    labels=list(
      conc="Plasma Concentration",
      Time="Time Since First Dose"),
    units=list(conc="mg/L", Time="hr"))
mydose <-
  PKNCAdose(
    data=data.frame(Subject=unique(Theoph$Subject),
      Time=0),
    formula=~Time|Subject,
    labels=list(Time="Time Since First Dose"),
    units=list(Dose="mg/kg"))

mydat <-
  PKNCAdata(data.conc=myconc,
            data.dose=mydose)

myres <- pk.nca(mydat)

summary(myres, simplify.start=FALSE)

myconc$data$conc[myconc$data$Time == 0] <- 0
tmp.superposition <- superposition(myconc, tau=24)

plot(mydat$conc)

q(save="no")
