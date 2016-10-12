context("Class generation-PKNCAdose")

library(dplyr)
source("generate.data.R")

test_that("PKNCAdose", {
  tmp.conc <- generate.conc(nsub=5, ntreat=2, time.points=0:24)
  tmp.conc.analyte <- generate.conc(nsub=5, ntreat=2, time.points=0:24,
                                    nanalytes=2)
  tmp.conc.study <- generate.conc(nsub=5, ntreat=2, time.points=0:24,
                                  nstudies=2)
  tmp.conc.analyte.study <- generate.conc(nsub=5, ntreat=2, time.points=0:24,
                                          nanalytes=2, nstudies=2)
  tmp.dose <- generate.dose(tmp.conc)
  rownames(tmp.dose) <- NULL
  tmp.dose.analyte <- generate.dose(tmp.conc.analyte)
  tmp.dose.study <- generate.dose(tmp.conc.study)
  tmp.dose.analyte.study <- generate.dose(tmp.conc.analyte.study)
  
  ## Variables present
  expect_error(PKNCAdose(tmp.dose, formula=dosea~time|treatment+ID),
               regexp="The left side formula must be a variable in the data, empty, or '.'.",
               info="All formula parameters must be in the data (LHS)")
  expect_error(PKNCAdose(tmp.dose, formula=dose~timea|treatment+ID),
               regexp="The right side formula must be a variable in the data or '.'.",
               info="All formula parameters must be in the data (RHS)")
  expect_error(PKNCAdose(tmp.dose, formula=dose~time|treatmenta+ID),
               regexp="All of the variables in the groups must be in the data",
               info="All formula parameters must be in the data (groups)")
  
  ## Number of variables
  expect_error(PKNCAdose(tmp.dose, formula=dose+ID~time|treatment+ID),
               regexp="The left side of the formula must have zero or one variable",
               info="The right number of parameters in the formula (LHS)")
  expect_error(PKNCAdose(tmp.dose, formula=dose~time+ID|treatment+ID),
               regexp="The right side of the formula \\(excluding groups\\) must have exactly one variable",
               info="The right number of parameters in the formula (RHS)")
  expect_error(PKNCAdose(tmp.dose, formula=~time+ID|treatment+ID),
               regexp="The right side of the formula \\(excluding groups\\) must have exactly one variable",
               info="The right number of parameters in the formula (RHS)")
  
  ## Accept "." on either side of the ~
  expect_equal(PKNCAdose(tmp.dose, formula=.~time|treatment+ID),
               structure(list(
                 data=cbind(tmp.dose,
                            data.frame(route="extravascular",
                                       duration=0,
                                       stringsAsFactors=FALSE)),
                 formula = . ~ time | treatment + ID,
                 route="route",
                 duration="duration"),
                 .Names = c("data", "formula", "route", "duration"),
                 class = c("PKNCAdose", "list")),
               info="PKNCAdose accepts . on the left side of the formula")
  expect_equal(PKNCAdose(tmp.dose, formula=dose~.|treatment+ID),
               structure(list(
                 data=cbind(tmp.dose,
                            data.frame(route="extravascular",
                                       duration=0,
                                       stringsAsFactors=FALSE)),
                 formula = dose ~ . | treatment + ID,
                 route="route",
                 duration="duration"),
                 .Names = c("data", "formula", "route", "duration"),
                 class = c("PKNCAdose", "list")),
               info="PKNCAdose accepts . on the right side of the formula")
  
  tmp.dose.na <- tmp.dose
  tmp.dose.na$time[1] <- NA
  expect_error(PKNCAdose(tmp.dose.na, formula=dose~time|treatment+ID),
               regex="Some but not all values are missing for the independent variable",
               info="Dose time must either all or none be NA.")
  ## Keys must be unique
  bad.dose.analyte <- unique(tmp.conc.analyte[,c("treatment", "ID", "analyte")])
  bad.dose.analyte$dose <- 1
  bad.dose.analyte$time <- 0
  expect_error(PKNCAdose(bad.dose.analyte, formula=dose~time|treatment+ID),
               regexp="Rows that are not unique per group and time",
               info="Duplicated key rows")
  
  expect_equal(PKNCAdose(tmp.dose,
                         formula=dose~time|treatment+ID),
               PKNCAdose(tbl_df(tmp.dose),
                         formula=dose~time|treatment+ID),
               info="tbl_df and data.frame classes both work and create identical objects")
})

test_that("PKNCAdose model.frame", {
  tmp.conc <- generate.conc(nsub=5, ntreat=2, time.points=0:24)
  tmp.conc.analyte <- generate.conc(nsub=5, ntreat=2, time.points=0:24,
                                    nanalytes=2)
  tmp.conc.study <- generate.conc(nsub=5, ntreat=2, time.points=0:24,
                                  nstudies=2)
  tmp.conc.analyte.study <- generate.conc(nsub=5, ntreat=2, time.points=0:24,
                                          nanalytes=2, nstudies=2)
  tmp.dose <- generate.dose(tmp.conc)

  mydose1 <- PKNCAdose(formula=dose~time|treatment+ID, data=tmp.dose)
  expect_equal(getDepVar.PKNCAdose(mydose1),
               rep(1:2, each=5),
               info="getDepVar.PKNCAdose works with two-sided formula")
  expect_equal(getIndepVar.PKNCAdose(mydose1),
               rep(0, 10),
               info="getIndepVar.PKNCAdose works with two-sided formula")
  expect_equal(model.frame.PKNCAdose(mydose1),
               data.frame("getDepVar.PKNCAdose(formula)"=rep(1:2, each=5),
                          "getIndepVar.PKNCAdose(formula)"=0,
                          treatment=rep(c("Trt 1", "Trt 2"), each=5),
                          ID=rep(1:5, 2),
                          stringsAsFactors=FALSE),
               check.attributes=FALSE,
               info="model.frame.PKNCAdose works with two-sided formula")
  
  mydose2 <- PKNCAdose(formula=~time|treatment+ID, data=tmp.dose)
  expect_equal(getDepVar.PKNCAdose(mydose2),
               rep(NA_integer_, 10),
               info="getDepVar.PKNCAdose works with one-sided formula")
  expect_equal(getIndepVar.PKNCAdose(mydose2),
               rep(0, 10),
               info="getIndepVar.PKNCAdose works with one-sided formula")
  expect_equal(model.frame.PKNCAdose(mydose2),
               data.frame("getDepVar.PKNCAdose(formula)"=NA_integer_,
                          "getIndepVar.PKNCAdose(formula)"=0,
                          treatment=rep(c("Trt 1", "Trt 2"), each=5),
                          ID=rep(1:5, 2),
                          stringsAsFactors=FALSE),
               check.attributes=FALSE,
               info="model.frame.PKNCAdose works with one-sided formula")

  mydose3 <- PKNCAdose(formula=.~time|treatment+ID, data=tmp.dose)
  expect_equal(getDepVar.PKNCAdose(mydose3),
               rep(NA_integer_, 10),
               info="getDepVar.PKNCAdose works with one-sided formula ('.' on LHS)")
  expect_equal(getIndepVar.PKNCAdose(mydose3),
               rep(0, 10),
               info="getIndepVar.PKNCAdose works with one-sided formula ('.' on LHS)")
  expect_equal(model.frame.PKNCAdose(mydose3),
               data.frame("getDepVar.PKNCAdose(formula)"=NA_integer_,
                          "getIndepVar.PKNCAdose(formula)"=0,
                          treatment=rep(c("Trt 1", "Trt 2"), each=5),
                          ID=rep(1:5, 2),
                          stringsAsFactors=FALSE),
               check.attributes=FALSE,
               info="model.frame.PKNCAdose works with one-sided formula ('.' on LHS)")
  
  mydose4 <- PKNCAdose(formula=dose~.|treatment+ID, data=tmp.dose)
  expect_equal(getDepVar.PKNCAdose(mydose4),
               rep(1:2, each=5),
               info="getDepVar.PKNCAdose works with one-sided formula ('.' on RHS)")
  expect_equal(getIndepVar.PKNCAdose(mydose4),
               rep(NA_integer_, 10),
               info="getIndepVar.PKNCAdose works with one-sided formula ('.' on RHS)")
  expect_equal(model.frame.PKNCAdose(mydose4),
               data.frame("getDepVar.PKNCAdose(formula)"=rep(1:2, each=5),
                          "getIndepVar.PKNCAdose(formula)"=NA_integer_,
                          treatment=rep(c("Trt 1", "Trt 2"), each=5),
                          ID=rep(1:5, 2),
                          stringsAsFactors=FALSE),
               check.attributes=FALSE,
               info="model.frame.PKNCAdose works with one-sided formula ('.' on RHS)")
  
  ## You can't give multiple rows per group if you don't give time.
  expect_error(PKNCAdose(formula=dose~.|treatment+ID, data=rbind(tmp.dose, tmp.dose)),
               regexp="Rows that are not unique per group and time.*found within dosing data",
               info="Dosing must have unique values with time and group")
})