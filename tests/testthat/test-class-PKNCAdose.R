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
                 data=tmp.dose,
                 formula = . ~ time | treatment + ID),
                 .Names = c("data", "formula"),
                 class = c("PKNCAdose", "list")),
               info="PKNCAdose accepts . on the left side of the formula")
  expect_equal(PKNCAdose(tmp.dose, formula=dose~.|treatment+ID),
               structure(list(
                 data=tmp.dose,
                 formula = dose ~ . | treatment + ID),
                 .Names = c("data", "formula"),
                 class = c("PKNCAdose", "list")),
               info="PKNCAdose accepts . on the right side of the formula")
  
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
