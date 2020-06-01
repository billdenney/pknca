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
  ## Data exists
  expect_error(PKNCAconc(data.frame()),
               regexp="data must have at least one row.",
               info="PKNCAconc requires data")
  
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
  
  expect_equal(
    PKNCAconc(
      tmp.conc.analyte,
      formula=conc~time|treatment+ID/analyte
    ),
    PKNCAconc(
      dplyr::as_tibble(tmp.conc.analyte),
      formula=conc~time|treatment+ID/analyte
    ),
    info="tibble and data.frame classes both work and create identical objects"
  )
})

test_that("PKNCAconc with input other than data.frames", {
  tmp <- structure(list(), class="foo")
  captured_message <- tryCatch(as.data.frame(tmp), error=function(e) e$message)
  expect_error(PKNCAconc(tmp, formula=conc~time|treatment+ID),
               regexp=captured_message,
               fixed=TRUE,
               info="Attempt to coerce into a data.frame.")
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
                           x_nodata <- x
                           x_nodata$data <- NULL
                           mc_nodata <- myconc
                           mc_nodata$data <- NULL
                           identical(x_nodata, mc_nodata)
                         })),
              info="All values (other than data) are accurately transferred.")
  expect_equal(split.PKNCAconc(NA),
               {
                 tmp <- list(NA)
                 attr(tmp, "groupid") <- data.frame(NA)[,c()]
                 tmp
               },
               info="NA split returns an effectively null split.")
  
  # There is a "feature" of base R split where NA values are ignored as
  # levels of the factor.  PKNCA works around this "feature".
  # This has 2 not 4 groups
  #
  # mydata <- data.frame(A=rep(c(NA_character_, "A"), each=4),
  #                      B=rep(1:2, 4),
  #                      C=11:18,
  #                      stringsAsFactors=FALSE)
  # split(mydata, f=mydata[,c("A", "B")])
  
  tmp_conc_na <- generate.conc(nsub=2, ntreat=2, time.points=0:24)
  tmp_conc_na$treatment[tmp_conc_na$treatment %in% "Trt 1"] <- NA_character_
  myconc_na <- PKNCAconc(tmp_conc_na, formula=conc~time|treatment+ID)
  tmp_myconc_na_split <- split(myconc_na)
  expect_equal(length(tmp_myconc_na_split),
               4,
               info="NA values in groups are kept not dropped")
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

test_that("PKNCAconc with exclusions", {
  tmp.conc <- generate.conc(nsub=2, ntreat=2, time.points=0:24)
  tmp.conc$excl <- NA_character_
  myconc <- PKNCAconc(tmp.conc, formula=conc~time|treatment+ID, exclude="excl")
  expect_equal(myconc,
               structure(
                 list(data=cbind(tmp.conc,
                                 volume=NA_real_,
                                 duration=0),
                      formula=conc~time|treatment+ID,
                      subject="ID", 
                      exclude="excl",
                      columns=list(volume="volume",
                                   duration="duration")),
                 class=c("PKNCAconc", "list")))
})

test_that("PKNCAconc with duration", {
  tmp.conc <- generate.conc(nsub=2, ntreat=2, time.points=0:24)
  tmp.conc$duration_test <- 0.1
  myconc <- PKNCAconc(tmp.conc,
                      formula=conc~time|treatment+ID,
                      duration="duration_test")
  expect_equal(myconc,
               structure(
                 list(data=cbind(tmp.conc,
                                 data.frame(exclude=NA_character_,
                                            volume=NA_real_,
                                            stringsAsFactors=FALSE)),
                      formula=conc~time|treatment+ID,
                      subject="ID", 
                      exclude="exclude",
                      columns=list(volume="volume",
                                   duration="duration_test")),
                 class=c("PKNCAconc", "list")))
})

test_that("PKNCAconc with nominal time added", {
  tmp.conc <- generate.conc(nsub=2, ntreat=2, time.points=0:24)
  tmp.conc$tnom <- tmp.conc$time
  myconc <- PKNCAconc(tmp.conc, formula=conc~time|treatment+ID, time.nominal="tnom")
  expect_equal(myconc,
               structure(
                 list(data=cbind(tmp.conc,
                                 data.frame(exclude=NA_character_,
                                            volume=NA_real_,
                                            duration=0,
                                            stringsAsFactors=FALSE)),
                      formula=conc~time|treatment+ID,
                      subject="ID",
                      exclude="exclude",
                      columns=list(volume="volume",
                                   duration="duration",
                                   time.nominal="tnom")),
                 class=c("PKNCAconc", "list")))
  expect_equal(PKNCAconc(tmp.conc, formula=conc~time|treatment+ID, time.nominal="foo"),
               structure(
                 list(data=cbind(tmp.conc,
                                 data.frame(exclude=NA_character_,
                                            volume=NA_real_,
                                            duration=0,
                                            foo=NA,
                                            stringsAsFactors=FALSE)),
                      formula=conc~time|treatment+ID,
                      subject="ID",
                      exclude="exclude",
                      columns=list(volume="volume",
                                   duration="duration",
                                   time.nominal="foo")),
                 class=c("PKNCAconc", "list")))
})

test_that("PKNCAconc with volume added", {
  tmp.conc <- generate.conc(nsub=2, ntreat=2, time.points=0:24)
  tmp.conc$vol <- 1:nrow(tmp.conc)
  myconc <- PKNCAconc(tmp.conc, formula=conc~time|treatment+ID, volume="vol")
  expect_equal(myconc,
               structure(
                 list(data=cbind(tmp.conc,
                                 data.frame(exclude=NA_character_,
                                            duration=0,
                                            stringsAsFactors=FALSE)),
                      formula=conc~time|treatment+ID,
                      subject="ID",
                      exclude="exclude",
                      columns=list(volume="vol",
                                   duration="duration")),
                 class=c("PKNCAconc", "list")))
  myconc_manual_vol <- PKNCAconc(tmp.conc, formula=conc~time|treatment+ID, volume=2)
  expect_equal(myconc_manual_vol,
               structure(
                 list(data=cbind(tmp.conc,
                                 data.frame(exclude=NA_character_,
                                            volume=2,
                                            duration=0,
                                            stringsAsFactors=FALSE)),
                      formula=conc~time|treatment+ID,
                      subject="ID",
                      exclude="exclude",
                      columns=list(volume="volume",
                                   duration="duration")),
                 class=c("PKNCAconc", "list")))
  myconc_manual_vol_vector <- PKNCAconc(tmp.conc, formula=conc~time|treatment+ID, volume=1:nrow(tmp.conc))
  expect_equal(myconc_manual_vol_vector,
               structure(
                 list(data=cbind(tmp.conc,
                                 data.frame(exclude=NA_character_,
                                            volume=1:nrow(tmp.conc),
                                            duration=0,
                                            stringsAsFactors=FALSE)),
                      formula=conc~time|treatment+ID,
                      subject="ID",
                      exclude="exclude",
                      columns=list(volume="volume",
                                   duration="duration")),
                 class=c("PKNCAconc", "list")))
})
