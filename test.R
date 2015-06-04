rm(list=ls())

## options(repos="http://cran.r-project.org/")
## install.packages("devtools")
## install.packages("testthat")
## install.packages("roxygen2")
## install.packages("plyr")
## install.packages("covr")

## library(nlme)
library(devtools)
library(covr)
library(testthat)
library(shiny)

load_all(".")
test(".")
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

mydat <-
  PKNCAdata(data.conc=PKNCAconc(
              data=as.data.frame(Theoph),
              formula=conc~Time|Subject,
              labels=list(
                conc="Plasma Concentration",
                Time="Time Since First Dose"),
              units=list(conc="mg/L", Time="hr")),
            data.dose=PKNCAdose(
              data=data.frame(Subject=unique(Theoph$Subject),
                Time=0),
              formula=~Time|Subject,
              labels=list(Time="Time Since First Dose"),
              units=list(Dose="mg/kg")))

myres <- pk.nca(mydat)

plot(mydat$conc)

q(save="no")
