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
               {
                 tmp.intervals <- merge(PKNCA.options("single.dose.aucs"), tmp.dose)
                 tmp.intervals <- tmp.intervals[order(tmp.intervals$treatment, tmp.intervals$ID),]
                 tmp.intervals$time <- NULL
                 tmp.intervals$dose <- NULL
                 tmp <- list(conc=obj.conc,
                             dose=obj.dose,
                             options=list(),
                             intervals=tmp.intervals)
                 class(tmp) <- c("PKNCAdata", "list")
                 tmp
               }, check.attributes=FALSE,
               info="Selection of single dose AUCs")
  
  tmp.conc <- generate.conc(nsub=5, ntreat=2, time.points=0:24)
  tmp.dose <- generate.dose(tmp.conc)
  tmp.conc <- tmp.conc[!(tmp.conc$ID %in% 1),]
  obj.conc <-
    PKNCAconc(tmp.conc, formula=conc~time|treatment+ID)
  obj.dose <- PKNCAdose(tmp.dose, formula=dose~time|treatment+ID)
  expect_warning(PKNCAdata(obj.conc, obj.dose),
                 regexp="No intervals generated due to no concentration data for treatment=Trt 1, ID=1",
                 info="Missing concentration data with dose data gives a warning.")
  
  expect_warning(PKNCAdata(obj.conc, obj.dose, formula.conc=a~b),
                 regexp="data.conc was given as a PKNCAconc object.  Ignoring formula.conc")
  expect_warning(PKNCAdata(obj.conc, obj.dose, formula.dose=a~b),
                 regexp="data.dose was given as a PKNCAdose object.  Ignoring formula.dose")
})

test_that("PKNCAdata with no or limited dose information", {
  tmp.conc <- generate.conc(nsub=5, ntreat=2, time.points=0:24)
  obj.conc <-
    PKNCAconc(tmp.conc, formula=conc~time|treatment+ID)

  expect_error(PKNCAdata(obj.conc),
               regexp="If data.dose is not given, intervals must be given",
               info="One of dose and intervals is required (no dose)")
  expect_error(PKNCAdata(obj.conc, data.dose=NA),
               regexp="If data.dose is not given, intervals must be given",
               info="One of dose and intervals is required (NA dose)")
  expect_equal(PKNCAdata(obj.conc, intervals=data.frame(start=0, end=24, aucinf.obs=TRUE)),
               {
                 tmp <- list(conc=obj.conc,
                             dose=NA,
                             options=list(),
                             intervals=check.interval.specification(
                               data.frame(start=0, end=24, aucinf.obs=TRUE)))
                 class(tmp) <- c("PKNCAdata", "list")
                 tmp
               })
  
  tmp.conc <- generate.conc(nsub=5, ntreat=2, time.points=0:24)
  tmp.dose <- generate.dose(tmp.conc)
  obj.conc <-
    PKNCAconc(tmp.conc, formula=conc~time|treatment+ID)
  obj.dose <- PKNCAdose(tmp.dose, formula=dose~.|treatment+ID)
  expect_error(PKNCAdata(obj.conc, obj.dose),
               regexp="Dose times were not given, so intervals must be manually specified.",
               info="No dose times requires intervals.")
})

test_that("print.PKNCAdata", {
  tmp.conc <- generate.conc(nsub=2, ntreat=2, time.points=0:24)
  tmp.dose <- generate.dose(tmp.conc)
  obj.conc <- PKNCAconc(tmp.conc, formula=conc~time|treatment+ID)
  obj.dose <- PKNCAdose(tmp.dose, formula=dose~time|treatment+ID)
  obj.data.nodose <- PKNCAdata(obj.conc,
                               intervals=data.frame(start=0, end=24, aucinf.obs=TRUE))
  obj.data.nodose.opt <-
    PKNCAdata(obj.conc,
              intervals=data.frame(start=0, end=24, aucinf.obs=TRUE),
              options=list(min.hl.r.squared=0.95))
  obj.data.dose <- PKNCAdata(obj.conc, data.dose=obj.dose)

  expect_output(print.PKNCAdata(obj.data.nodose),
                regexp="Formula for concentration:
 conc ~ time | treatment + ID
With 2 subjects defined in the 'ID' column.
Nominal time column is not specified.
                
First 6 rows of concentration data:
 treatment ID time      conc exclude
     Trt 1  1    0 0.0000000    <NA>
     Trt 1  1    1 0.7052248    <NA>
     Trt 1  1    2 0.7144320    <NA>
     Trt 1  1    3 0.8596094    <NA>
     Trt 1  1    4 0.9998126    <NA>
     Trt 1  1    5 0.7651474    <NA>
No dosing information.

With 1 rows of AUC specifications.
No options are set differently than default.",
                info="Generic print.PKNCAdata works with no dosing")
  expect_output(print.PKNCAdata(obj.data.dose),
                regexp="Formula for concentration:
 conc ~ time | treatment + ID
With 2 subjects defined in the 'ID' column.
Nominal time column is not specified.
                
First 6 rows of concentration data:
 treatment ID time      conc exclude
     Trt 1  1    0 0.0000000    <NA>
     Trt 1  1    1 0.7052248    <NA>
     Trt 1  1    2 0.7144320    <NA>
     Trt 1  1    3 0.8596094    <NA>
     Trt 1  1    4 0.9998126    <NA>
     Trt 1  1    5 0.7651474    <NA>
Formula for dosing:
 dose ~ time | treatment + ID
Nominal time column is not specified.
                
Data for dosing:
 treatment ID dose time exclude         route duration
     Trt 1  1    1    0    <NA> extravascular        0
     Trt 1  2    1    0    <NA> extravascular        0
     Trt 2  1    2    0    <NA> extravascular        0
     Trt 2  2    2    0    <NA> extravascular        0
With 1 rows of AUC specifications.
No options are set differently than default.",
                info="Generic print.PKNCAdata works with dosing")

  expect_output(print.PKNCAdata(obj.data.nodose.opt),
                regexp="Formula for concentration:
 conc ~ time | treatment + ID
With 2 subjects defined in the 'ID' column.
Nominal time column is not specified.
                
First 6 rows of concentration data:
 treatment ID time      conc exclude
     Trt 1  1    0 0.0000000    <NA>
     Trt 1  1    1 0.7052248    <NA>
     Trt 1  1    2 0.7144320    <NA>
     Trt 1  1    3 0.8596094    <NA>
     Trt 1  1    4 0.9998126    <NA>
     Trt 1  1    5 0.7651474    <NA>
No dosing information.

With 1 rows of AUC specifications.
Options changed from default are:
$min.hl.r.squared
[1] 0.95",
                info="Generic print.PKNCAdata works with no dosing and with options changed")
})

test_that("summary.PKNCAdata", {
  tmp.conc <- generate.conc(nsub=2, ntreat=2, time.points=0:24)
  tmp.dose <- generate.dose(tmp.conc)
  obj.conc <- PKNCAconc(tmp.conc, formula=conc~time|treatment+ID)
  obj.dose <- PKNCAdose(tmp.dose, formula=dose~time|treatment+ID)
  obj.data.nodose <- PKNCAdata(obj.conc,
                               intervals=data.frame(start=0, end=24, aucinf.obs=TRUE))
  
  expect_output(summary(obj.data.nodose),
                regexp="Formula for concentration:
 conc ~ time | treatment + ID
With 2 subjects defined in the 'ID' column.
Nominal time column is not specified.

Group summary:
 Group Name Count
  treatment     2
         ID     4

First 6 rows of concentration data:
 treatment ID time      conc exclude
     Trt 1  1    0 0.0000000    <NA>
     Trt 1  1    1 0.7052248    <NA>
     Trt 1  1    2 0.7144320    <NA>
     Trt 1  1    3 0.8596094    <NA>
     Trt 1  1    4 0.9998126    <NA>
     Trt 1  1    5 0.7651474    <NA>
No dosing information.

With 1 rows of AUC specifications.
No options are set differently than default.",
                info="Generic summary.PKNCAdata works.")
})

test_that("splitting PKNCAdata", {
  tmp.conc <- generate.conc(nsub=2, ntreat=1, time.points=0:24)
  obj.conc <-
    PKNCAconc(tmp.conc, formula=conc~time|treatment+ID)
  intervals <- data.frame(start=0, end=24, aucinf.obs=TRUE)
  mydata <- PKNCAdata(obj.conc, intervals=intervals)

  splitconc <- split(obj.conc)
  expect_equal(split(mydata),
               list(
                 "Trt 1\n1"=list(
                   conc=splitconc[[1]],
                   dose=NA,
                   intervals=check.interval.specification(data.frame(start=0, end=24, aucinf.obs=TRUE)),
                   options=list()),
                 "Trt 1\n2"=list(
                   conc=splitconc[[2]],
                   dose=NA,
                   intervals=check.interval.specification(data.frame(start=0, end=24, aucinf.obs=TRUE)),
                   options=list())))
})

test_that("no intervals auto-determined (Fix GitHub issue #84)", {
  tmp_conc <-
    data.frame(
      Subject=1,
      Treatment=c(1, rep(2, 6)),
      Time=c(0, 1:6),
      Conc=1
    )
  tmp_dose <-
    data.frame(
      Subject=1,
      Treatment=c(1, 2),
      Time=c(0, 1),
      Dose=1
    )
  
  interval_1 <- PKNCA.options("single.dose.aucs")[c(1:2, 1:2),]
  interval_1$start <- rep(0:1, each=2)
  interval_1$end <- c(interval_1$end[1:2], interval_1$end[3:4] + 1)
  interval_1 <- cbind(interval_1, data.frame(Treatment=rep(1:2, each=2), Subject=1))
  two_single_dose_treatments <-
    PKNCAdata(
      PKNCAconc(data=tmp_conc, Conc~Time|Treatment+Subject),
      PKNCAdose(data=tmp_dose, Dose~Time|Treatment+Subject)
    )
  expect_equal(
    two_single_dose_treatments$intervals,
    interval_1,
    check.attributes=FALSE
  )
  interval_2 <-
    check.interval.specification(
      data.frame(start=1, end=c(2, Inf),
                 auclast=c(TRUE, FALSE),
                 cmax=c(TRUE, FALSE),
                 tmax=c(TRUE, FALSE),
                 half.life=c(FALSE, TRUE),
                 Treatment=2,
                 Subject=1
      )
    )
  expect_warning(
    two_multiple_dose_treatments <-
      PKNCAdata(
        PKNCAconc(data=tmp_conc, Conc~Time|Treatment+Subject),
        PKNCAdose(data=tmp_dose, Dose~Time|Subject)
      ),
    regexp="No intervals generated"
  )
  expect_equal(
    two_multiple_dose_treatments$intervals,
    interval_2
  )
})

test_that("Ensure that unexpected arguments to PKNCAdata give an error (related to issue #83)", {
  tmp.conc <- generate.conc(nsub=2, ntreat=1, time.points=0:24)
  tmp.dose <- generate.dose(tmp.conc)
  obj.conc <-
    PKNCAconc(tmp.conc, formula=conc~time|treatment+ID)
  obj.dose <-
    PKNCAdose(tmp.dose, formula=dose~time|treatment+ID)
  expect_error(mydata <- PKNCAdata(obj.conc, obj.dose, 1),
               regexp="Unknown argument")
})

test_that("intervals may be a tibble", {
  tmp.conc <- generate.conc(nsub=2, ntreat=1, time.points=0:24)
  tmp.dose <- generate.dose(tmp.conc)
  obj.conc <-
    PKNCAconc(tmp.conc, formula=conc~time|treatment+ID)
  obj.dose <-
    PKNCAdose(tmp.dose, formula=dose~time|treatment+ID)
  intervals <- data.frame(start=0, end=24, aucinf.obs=TRUE)
  mydata_tibble <- PKNCAdata(obj.conc, obj.dose, intervals=dplyr::as_tibble(intervals))
  mydata <- PKNCAdata(obj.conc, obj.dose, intervals=intervals)
  expect_equal(
    as.data.frame(pk.nca(mydata_tibble)),
    as.data.frame(pk.nca(mydata))
  )
})
