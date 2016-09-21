context("Class generation-PKNCAconc")

library(dplyr)
source("generate.data.R")

test_that("PKNCAconc", {
  tmp.conc <- generate.conc(nsub=5, ntreat=2, time.points=0:24)
  tmp.conc.analyte <- generate.conc(nsub=5, ntreat=2, time.points=0:24,
                                    nanalytes=2)
  tmp.conc.study <- generate.conc(nsub=5, ntreat=2, time.points=0:24,
                                  nstudies=2)
  tmp.conc.analyte.study <- generate.conc(nsub=5, ntreat=2, time.points=0:24,
                                          nanalytes=2, nstudies=2)
  
  ## Variables present
  expect_error(PKNCAconc(tmp.conc, formula=conca~time|treatment+ID),
               regexp="All of the variables in the formula must be in the data",
               info="All formula parameters must be in the data (LHS)")
  expect_error(PKNCAconc(tmp.conc, formula=conc~timea|treatment+ID),
               regexp="All of the variables in the formula must be in the data",
               info="All formula parameters must be in the data (RHS)")
  expect_error(PKNCAconc(tmp.conc, formula=conc~time|treatmenta+ID),
               regexp="All of the variables in the formula must be in the data",
               info="All formula parameters must be in the data (groups)")
  
  ## Number of variables
  expect_error(PKNCAconc(tmp.conc, formula=conc+ID~time|treatment+ID),
               regexp="The left hand side of the formula must have exactly one variable",
               info="The right number of parameters in the formula (LHS)")
  expect_error(PKNCAconc(tmp.conc, formula=conc~time+ID|treatment+ID),
               regexp="The right hand side of the formula \\(excluding groups\\) must have exactly one variable",
               info="The right number of parameters in the formula (RHS)")
  
  ## Subject assignment
  expect_equal(PKNCAconc(tmp.conc.analyte, formula=conc~time|treatment+ID/analyte),
               PKNCAconc(tmp.conc.analyte, formula=conc~time|treatment+ID/analyte, subject="ID"))
  expect_error(PKNCAconc(tmp.conc.analyte, formula=conc~time|treatment+ID/analyte, subject=5),
               regexp="subject must be a character string")
  expect_error(PKNCAconc(tmp.conc.analyte, formula=conc~time|treatment+ID/analyte, subject=c("", "foo")),
               regexp="subject must be a scalar")
  expect_error(PKNCAconc(tmp.conc.analyte, formula=conc~time|treatment+ID/analyte, subject="foo"),
               regexp="The subject parameter must map to a name in the data")
  
  ## Keys must be unique
  expect_error(PKNCAconc(tmp.conc.analyte, formula=conc~time|treatment+ID),
               regexp="Rows that are not unique per group and time",
               info="Duplicated key rows")
  
  expect_equal(PKNCAconc(tmp.conc.analyte,
                         formula=conc~time|treatment+ID/analyte),
               PKNCAconc(tbl_df(tmp.conc.analyte),
                         formula=conc~time|treatment+ID/analyte),
               info="tbl_df and data.frame classes both work and create identical objects")
})
