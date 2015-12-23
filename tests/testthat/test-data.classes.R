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

test_that("PKNCAresults and summary", {
  ## Note that generate.conc sets the random seed, so it doesn't have
  ## to happen here.
  tmpconc <- generate.conc(2, 1, 0:24)
  tmpdose <- generate.dose(tmpconc)
  myconc <- PKNCAconc(tmpconc, conc~time|treatment+ID)
  mydose <- PKNCAdose(tmpdose, dose~time|treatment+ID)
  mydata <- PKNCAdata(myconc, mydose)
  myresult <- pk.nca(mydata)

  expect_equal(names(myresult),
               c("result", "data", "provenance"),
               info="Make sure that the result has the expected names (and only the expected names) in it.")
  
  ## Test each of the pieces for myresult for accuracy

  expect_equal(myresult$data, {
    tmp <- mydata
    ## The options should be the default options after the
    ## calculations are done.
    tmp$options <- PKNCA.options()
    tmp
  }, info="The data is just a copy of the input data plus an instantiation of the PKNCA.options")

  ## The specific order of the levels isn't important-- the fact that
  ## they are factors and that the set doesn't change is important.
  test.code.levels <- levels(myresult$result$PPTESTCD)
  verify.result <-
    data.frame(
      start=0,
      end=c(24, rep(Inf, 12),
            24, rep(Inf, 12)),
      treatment="Trt 1",
      ID=rep(c(1, 2), each=13),
      PPTESTCD=factor(rep(c("auclast", "cmax", "tmax", "tlast",
                            "lambda.z", "r.squared", "adj.r.squared",
                            "lambda.z.time.first", "lambda.z.n.points",
                            "clast.pred", "half.life", "span.ratio",
                            "aucinf"),
                          times=2),
                      levels=test.code.levels),
      PPORRES=c(13.54, 0.9998, 4.000, 24.00,
                0.04297, 0.9072, 0.9021, 5.000,
                20.00, 0.3356, 16.13, 1.178,
                21.55, 14.03, 0.9410, 2.000,
                24.00, 0.05689, 0.9000, 0.8944,
                5.000, 20.00, 0.3011, 12.18,
                1.560, 19.56),
      stringsAsFactors=FALSE)
  expect_equal(myresult$result, verify.result,
               tol=0.001)

  ## Testing the summarization
  mysummary <- summary(myresult)
  expect_true(is.data.frame(mysummary))
  expect_equal(mysummary,
               data.frame(start=0,
                          end=c(24, Inf),
                          treatment="Trt 1",
                          auclast=c("13.8 [2.51]", "."),
                          cmax=c(".", "0.970 [4.29]"),
                          tmax=c(".", "3.00 [2.00, 4.00]"),
                          half.life=c(".", "14.2 [2.79]"),
                          aucinf=c(".", "20.5 [6.84]"),
                          stringsAsFactors=FALSE),
               info="simple summary of PKNCAresults performs as expected")

  tmpconc <- generate.conc(2, 1, 0:24)
  tmpconc$conc[tmpconc$ID %in% 2] <- 0
  tmpdose <- generate.dose(tmpconc)
  myconc <- PKNCAconc(tmpconc, conc~time|treatment+ID)
  mydose <- PKNCAdose(tmpdose, dose~time|treatment+ID)
  mydata <- PKNCAdata(myconc, mydose)
  expect_warning(myresult <- pk.nca(mydata),
                 regexp="Too few points for half-life calculation")
  mysummary <- summary(myresult)
  expect_equal(mysummary,
               data.frame(start=0,
                          end=c(24, Inf),
                          treatment="Trt 1",
                          auclast=c("13.5 [NC]", "."),
                          cmax=c(".", "1.00 [NC]"),
                          tmax=c(".", "4.00 [4.00, 4.00]"),
                          half.life=c(".", "16.1 [NC]"),
                          aucinf=c(".", "21.5 [NC]"),
                          stringsAsFactors=FALSE),
               info="summary of PKNCAresults with some missing values results in NA for spread")

  tmpconc <- generate.conc(2, 1, 0:24)
  tmpconc$conc <- 0
  tmpdose <- generate.dose(tmpconc)
  myconc <- PKNCAconc(tmpconc, conc~time|treatment+ID)
  mydose <- PKNCAdose(tmpdose, dose~time|treatment+ID)
  mydata <- PKNCAdata(myconc, mydose)
  expect_warning(myresult <- pk.nca(mydata),
                 regexp="Too few points for half-life calculation")
  mysummary <- summary(myresult)
  expect_equal(mysummary,
               data.frame(start=0,
                          end=c(24, Inf),
                          treatment="Trt 1",
                          auclast=c("NC", "."),
                          cmax=c(".", "NC"),
                          tmax=c(".", "NC"),
                          half.life=c(".", "NC"),
                          aucinf=c(".", "NC"),
                          stringsAsFactors=FALSE),
               info="summary of PKNCAresults without most results gives NC")

  mysummary <- summary(myresult,
                       not.requested.string="NR",
                       not.calculated.string="NoCalc")
  expect_equal(mysummary,
               data.frame(start=0,
                          end=c(24, Inf),
                          treatment="Trt 1",
                          auclast=c("NoCalc", "NR"),
                          cmax=c("NR", "NoCalc"),
                          tmax=c("NR", "NoCalc"),
                          half.life=c("NR", "NoCalc"),
                          aucinf=c("NR", "NoCalc"),
                          stringsAsFactors=FALSE),
               info="Summary respects the not.requested.string and not.calculated.string")
})
