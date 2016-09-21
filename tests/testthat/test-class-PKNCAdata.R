context("Class generation-PKNCAdata")

library(dplyr)
source("generate.data.R")

test_that("PKNCAdata", {
  tmp.conc <- generate.conc(nsub=5, ntreat=2, time.points=0:24)
  tmp.conc.analyte <- generate.conc(nsub=5, ntreat=2, time.points=0:24,
                                    nanalytes=2)
  tmp.conc.study <- generate.conc(nsub=5, ntreat=2, time.points=0:24,
                                  nstudies=2)
  tmp.conc.analyte.study <- generate.conc(nsub=5, ntreat=2, time.points=0:24,
                                          nanalytes=2, nstudies=2)
  tmp.dose <- generate.dose(tmp.conc)
  tmp.dose.analyte <- generate.dose(tmp.conc.analyte)
  tmp.dose.study <- generate.dose(tmp.conc.study)
  tmp.dose.analyte.study <- generate.dose(tmp.conc.analyte.study)
  obj.conc <-
    PKNCAconc(tmp.conc, formula=conc~time|treatment+ID)
  obj.conc.analyte <-
    PKNCAconc(tmp.conc.analyte,
              formula=conc~time|treatment+ID/analyte)
  obj.conc.study <-
    PKNCAconc(tmp.conc.study,
              formula=conc~time|study+treatment+ID)
  obj.conc.analyte.study <-
    PKNCAconc(tmp.conc.analyte.study,
              formula=conc~time|study+treatment+ID/analyte)
  
  obj.dose <- PKNCAdose(tmp.dose, formula=dose~time|treatment+ID)
  obj.dose.analyte <- PKNCAdose(tmp.dose.analyte, formula=dose~time|treatment+ID)
  obj.dose.study <- PKNCAdose(tmp.dose.study, formula=dose~time|study+treatment+ID)
  obj.dose.analyte.study <- PKNCAdose(tmp.dose.analyte.study, formula=dose~time|study+treatment+ID)
  
  expect_equal(PKNCAdata(obj.conc, obj.dose),
               PKNCAdata(obj.dose, obj.conc),
               info="Input arguments are reversible")
  expect_equal(PKNCAdata(obj.conc.analyte, obj.dose),
               PKNCAdata(obj.dose, obj.conc.analyte),
               info="Combination of dose and analyte works")
  expect_equal(PKNCAdata(data.conc=tmp.conc, formula.conc=conc~time|treatment+ID,
                         data.dose=tmp.dose, formula.dose=dose~time|treatment+ID),
               PKNCAdata(obj.conc, obj.dose),
               info="Concentration and dose data can be created on the fly")
  
  ## Input checking
  expect_error(PKNCAdata(obj.conc, obj.dose, options="a"),
               regexp="options must be a list.",
               info="Option class")
  expect_error(PKNCAdata(obj.conc, obj.dose, options=list(1)),
               regexp="options must have names.",
               info="Option structure")
  expect_error(PKNCAdata(obj.conc, obj.dose, options=list(foo=1)),
               regexp="Invalid setting for PKNCA.*foo",
               info="Option names")
  
  ## Single dose AUCs are appropriately selected
  expect_equal(PKNCAdata(obj.conc, obj.dose),
               {tmp.intervals <- merge(PKNCA.options("single.dose.aucs"), tmp.dose)
               tmp.intervals$time <- NULL
               tmp.intervals$dose <- NULL
               tmp <- list(conc=obj.conc,
                           dose=obj.dose,
                           options=list(),
                           intervals=tmp.intervals)
               class(tmp) <- c("PKNCAdata", "list")
               tmp
               },
               info="Selection of single dose AUCs")
})
