context("Class generation")

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
})

test_that("PKNCAdose", {
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
  
  ## Variables present
  expect_error(PKNCAdose(tmp.dose, formula=dosea~time|treatment+ID),
               regexp="All of the variables in the formula must be in the data",
               info="All formula parameters must be in the data (LHS)")
  expect_error(PKNCAdose(tmp.dose, formula=dose~timea|treatment+ID),
               regexp="All of the variables in the formula must be in the data",
               info="All formula parameters must be in the data (RHS)")
  expect_error(PKNCAdose(tmp.dose, formula=dose~time|treatmenta+ID),
               regexp="All of the variables in the formula must be in the data",
               info="All formula parameters must be in the data (groups)")

  ## Number of variables
  expect_error(PKNCAdose(tmp.dose, formula=dose+ID~time|treatment+ID),
               regexp="The left hand side of the formula must have zero or one variable",
               info="The right number of parameters in the formula (LHS)")
  expect_error(PKNCAdose(tmp.dose, formula=dose~time+ID|treatment+ID),
               regexp="The right hand side of the formula \\(excluding groups\\) must have exactly one variable",
               info="The right number of parameters in the formula (RHS)")
  expect_error(PKNCAdose(tmp.dose, formula=~time+ID|treatment+ID),
               regexp="The right hand side of the formula \\(excluding groups\\) must have exactly one variable",
               info="The right number of parameters in the formula (RHS)")

  ## Keys must be unique
  bad.dose.analyte <- unique(tmp.conc.analyte[,c("treatment", "ID", "analyte")])
  bad.dose.analyte$dose <- 1
  bad.dose.analyte$time <- 0
  expect_error(PKNCAdose(bad.dose.analyte, formula=dose~time|treatment+ID),
               regexp="Rows that are not unique per group and time",
               info="Duplicated key rows")
})

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

test_that("findOverlap", {
  ## Calculation checks
  expect_equal(findOverlap(c(0, 1, 2), c(1, 2, 3)),
               c(1, 2, 3),
               info="no overlap")
  expect_equal(findOverlap(c(0, 0, 1), c(0.5, 1, 2)),
               c(1, 1, 2),
               info="some overlap")
  expect_equal(findOverlap(c(0, 1, 2), c(1.1, 2.1, 3)),
               c(1, 1, 1),
               info="chaining overlap")
  expect_equal(findOverlap(as.numeric(c(NA, NA, NA)),
                           as.numeric(c(NA, NA, NA))),
               rep(NA, 3),
               info="all missing")
  expect_equal(findOverlap(c(NA, 1, 2), c(1, NA, NA)),
               rep(NA, 3),
               info="all missing but split by start and end")

  ## Input checks
  expect_error(findOverlap(1, 1:2),
               regexp="start and end must be the same length")
  expect_error(findOverlap("a", 1),
               regexp="start must be numeric")
  expect_error(findOverlap(1, "a"),
               regexp="end must be numeric")
})
