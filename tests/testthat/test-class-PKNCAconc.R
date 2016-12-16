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
  expect_error(PKNCAconc(tmp.conc, formula=XXX~time|treatment+ID),
               regexp="All of the variables in the formula must be in the data",
               info="All formula parameters must be in the data (LHS)")
  expect_error(PKNCAconc(tmp.conc, formula=conc~XXX|treatment+ID),
               regexp="All of the variables in the formula must be in the data",
               info="All formula parameters must be in the data (RHS)")
  expect_error(PKNCAconc(tmp.conc, formula=conc~time|XXX+ID),
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

test_that("model frame and parameter extractions", {
  tmp.conc <- generate.conc(nsub=5, ntreat=2, time.points=0:24)
  expect_equal(model.frame(PKNCAconc(tmp.conc, formula=conc~time|treatment+ID)),
               tmp.conc[,c("conc", "time", "treatment", "ID")],
               info="model.frame.PKNCAconc extracts the correct components")
  expect_equal(getDepVar(PKNCAconc(tmp.conc, formula=conc~time|treatment+ID)),
               tmp.conc$conc,
               info="getDepVar.PKNCAconc extracts the correct component")
  expect_equal(getIndepVar(PKNCAconc(tmp.conc, formula=conc~time|treatment+ID)),
               tmp.conc$time,
               info="getIndepVar.PKNCAconc extracts the correct component")
  expect_equal(getGroups.PKNCAconc(PKNCAconc(tmp.conc, formula=conc~time|treatment+ID)),
               tmp.conc[,c("treatment", "ID")],
               info="getGroups.PKNCAconc extracts the correct components")

  tmp.conc <- generate.conc(nsub=5, ntreat=1, time.points=0:24)
  expect_equal(getGroups.PKNCAconc(PKNCAconc(tmp.conc, formula=conc~time|ID)),
               tmp.conc[,"ID", drop=FALSE],
               info="getGroups.PKNCAconc returns a data.frame even with a single grouping level")
  # Levels
  tmp.conc <- generate.conc(nsub=5, ntreat=2, time.points=0:24)
  expect_equal(getGroups.PKNCAconc(PKNCAconc(tmp.conc, formula=conc~time|treatment+ID), level=1),
               tmp.conc[,"treatment", drop=FALSE],
               info="getGroups.PKNCAconc the correct level (numeric scalar)")
  expect_equal(getGroups.PKNCAconc(PKNCAconc(tmp.conc, formula=conc~time|treatment+ID), level=-1),
               tmp.conc[,"ID", drop=FALSE],
               info="getGroups.PKNCAconc the correct level (negative numeric scalar)")
  expect_equal(getGroups.PKNCAconc(PKNCAconc(tmp.conc, formula=conc~time|treatment+ID), level=1:2),
               tmp.conc[,c("treatment", "ID"), drop=FALSE],
               info="getGroups.PKNCAconc the correct level (numeric vector)")
  expect_equal(getGroups.PKNCAconc(PKNCAconc(tmp.conc, formula=conc~time|treatment+ID), level=2:1),
               tmp.conc[,c("ID", "treatment"), drop=FALSE],
               info="getGroups.PKNCAconc the correct level (numeric vector)")
  expect_equal(getGroups.PKNCAconc(PKNCAconc(tmp.conc, formula=conc~time|treatment+ID), level="ID"),
               tmp.conc[,"ID", drop=FALSE],
               info="getGroups.PKNCAconc the correct level (character string scalar)")
  expect_equal(getGroups.PKNCAconc(PKNCAconc(tmp.conc, formula=conc~time|treatment+ID), level=c("ID", "treatment")),
               tmp.conc[,c("ID", "treatment"), drop=FALSE],
               info="getGroups.PKNCAconc the correct level (character string vector)")
  expect_error(getGroups.PKNCAconc(PKNCAconc(tmp.conc, formula=conc~time|treatment+ID), level="foo"),
               regexp="Not all levels are listed in the group names",
               info="getGroups.PKNCAconc gives an error if a group name is not present")
})

test_that("split.PKNCAconc", {
  tmp.conc <- generate.conc(nsub=2, ntreat=2, time.points=0:24)
  myconc <- PKNCAconc(tmp.conc, formula=conc~time|treatment+ID)
  expect_equal(base::split(myconc), split.PKNCAconc(myconc),
               info="The generic is correctly called")
  tmpsplit <- split.PKNCAconc(myconc)
  expect_true(all(sapply(tmpsplit,
                         function(x) {
                           all(names(x) == names(myconc))
                         })),
               info="All parameter names are accurately transferred")
  expect_true(all(sapply(tmpsplit,
                         function(x) {
                           ret <- TRUE
                           for (n in setdiff(names(x), "data")) {
                             ret <- ret & x[[n]] == myconc[[n]]
                           }
                           ret
                         })),
              info="All values (other than data) are accurately transferred.")
})

test_that("print.PKNCAconc", {
  tmp.conc <- generate.conc(nsub=2, ntreat=2, time.points=0:24)
  myconc <- PKNCAconc(tmp.conc, formula=conc~time|treatment+ID)

  expect_output(print.PKNCAconc(myconc),
                regexp="Formula for concentration:
 conc ~ time | treatment + ID
                With 2 subjects defined in the 'ID' column.
                Nominal time column is not specified.
                
                First 6 rows of concentration data:
                treatment ID time      conc
                Trt 1  1    0 0.0000000
                Trt 1  1    1 0.7052248
                Trt 1  1    2 0.7144320
                Trt 1  1    3 0.8596094
                Trt 1  1    4 0.9998126
                Trt 1  1    5 0.7651474",
                info="Generic print.PKNCAconc works")
  expect_output(print.PKNCAconc(myconc, n=0),
                regexp="Formula for concentration:
 conc ~ time | treatment + ID
With 2 subjects defined in the 'ID' column.
Nominal time column is not specified.",
                info="print.PKNCAconc respects the n argument.")
  expect_output(print.PKNCAconc(myconc, n=-98),
                regexp="Formula for concentration:
 conc ~ time | treatment + ID
                With 2 subjects defined in the 'ID' column.
                Nominal time column is not specified.
                
                First 2 rows of concentration data:
                treatment ID time      conc
                Trt 1  1    0 0.0000000
                Trt 1  1    1 0.7052248",
                info="print.PKNCAconc accurately uses negative n argument.")
})

test_that("summary.PKNCAconc", {
  tmp.conc <- generate.conc(nsub=2, ntreat=2, time.points=0:24)
  myconc <- PKNCAconc(tmp.conc, formula=conc~time|treatment+ID)
  
  expect_output(summary(myconc),
                regexp="Formula for concentration:
 conc ~ time | treatment + ID
                With 2 subjects defined in the 'ID' column.
                Nominal time column is not specified.
                
                Group summary:
                Group Name Count
                treatment     2
                ID     4",
                info="Generic summary.PKNCAconc works.")
})
