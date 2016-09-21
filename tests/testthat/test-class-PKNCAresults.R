context("Class generation-PKNCAresults")

library(dplyr)
source("generate.data.R")

test_that("PKNCAresults generation", {
  ## Note that generate.conc sets the random seed, so it doesn't have
  ## to happen here.
  tmpconc <- generate.conc(2, 1, 0:24)
  tmpdose <- generate.dose(tmpconc)
  myconc <- PKNCAconc(tmpconc, formula=conc~time|treatment+ID)
  mydose <- PKNCAdose(tmpdose, formula=dose~time|treatment+ID)
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
  
  test.code.levels <- levels(myresult$result$PPTESTCD)
  verify.result <-
    data.frame(
      start=0,
      end=c(24, rep(Inf, 13),
            24, rep(Inf, 13)),
      treatment="Trt 1",
      ID=rep(c(1, 2), each=14),
      PPTESTCD=factor(rep(c("auclast", "cmax", "tmax", "tlast", "clast.obs",
                            "lambda.z", "r.squared", "adj.r.squared",
                            "lambda.z.time.first", "lambda.z.n.points",
                            "clast.pred", "half.life", "span.ratio",
                            "aucinf.obs"),
                          times=2),
                      levels=test.code.levels),
      PPORRES=c(13.54, 0.9998, 4.000, 24.00, 0.3441,
                0.04297, 0.9072, 0.9021, 5.000,
                20.00, 0.3356, 16.13, 1.178,
                21.55, 14.03, 0.9410, 2.000,
                24.00, 0.3148, 0.05689, 0.9000, 0.8944,
                5.000, 20.00, 0.3011, 12.18,
                1.560, 19.56),
      stringsAsFactors=FALSE)
  expect_equal(myresult$result, verify.result,
               tol=0.001,
               info="The specific order of the levels isn't important-- the fact that they are factors and that the set doesn't change is important.")
  
  ## Test conversion to a data.frame
  expect_equal(as.data.frame(myresult), verify.result, tol=0.001,
               info="Conversion of PKNCAresults to a data.frame in long format (default long format)")
  expect_equal(as.data.frame(myresult), verify.result, tol=0.001,
               info="Conversion of PKNCAresults to a data.frame in long format (specifying long format)")
  expect_equal(as.data.frame(myresult, out.format="wide"),
               tidyr::spread_(verify.result, "PPTESTCD", "PPORRES"),
               tol=0.001,
               info="Conversion of PKNCAresults to a data.frame in wide format (specifying wide format)")
})

test_that("PKNCAresults summary", {
  ## Note that generate.conc sets the random seed, so it doesn't have
  ## to happen here.
  tmpconc <- generate.conc(2, 1, 0:24)
  tmpdose <- generate.dose(tmpconc)
  myconc <- PKNCAconc(tmpconc, formula=conc~time|treatment+ID)
  mydose <- PKNCAdose(tmpdose, formula=dose~time|treatment+ID)
  mydata <- PKNCAdata(myconc, mydose)
  myresult <- pk.nca(mydata)
  
  ## Testing the summarization
  mysummary <- summary(myresult)
  expect_true(is.data.frame(mysummary))
  expect_equal(mysummary,
               data.frame(start=0,
                          end=c(24, Inf),
                          treatment="Trt 1",
                          N="2",
                          auclast=c("13.8 [2.51]", "."),
                          cmax=c(".", "0.970 [4.29]"),
                          tmax=c(".", "3.00 [2.00, 4.00]"),
                          half.life=c(".", "14.2 [2.79]"),
                          aucinf.obs=c(".", "20.5 [6.84]"),
                          stringsAsFactors=FALSE),
               info="simple summary of PKNCAresults performs as expected")
  
  tmpconc <- generate.conc(2, 1, 0:24)
  tmpconc$conc[tmpconc$ID %in% 2] <- 0
  tmpdose <- generate.dose(tmpconc)
  myconc <- PKNCAconc(tmpconc, conc~time|treatment+ID)
  mydose <- PKNCAdose(tmpdose, dose~time|treatment+ID)
  mydata <- PKNCAdata(myconc, mydose)
  # Not capturing the warning due to R bug
  # https://bugs.r-project.org/bugzilla3/show_bug.cgi?id=17122
  #expect_warning(myresult <- pk.nca(mydata),
  #               regexp="Too few points for half-life calculation")
  myresult <- pk.nca(mydata)
  mysummary <- summary(myresult)
  expect_equal(mysummary,
               data.frame(start=0,
                          end=c(24, Inf),
                          treatment="Trt 1",
                          N="2",
                          auclast=c("13.5 [NC]", "."),
                          cmax=c(".", "1.00 [NC]"),
                          tmax=c(".", "4.00 [4.00, 4.00]"),
                          half.life=c(".", "16.1 [NC]"),
                          aucinf.obs=c(".", "21.5 [NC]"),
                          stringsAsFactors=FALSE),
               info="summary of PKNCAresults with some missing values results in NA for spread")
  
  tmpconc <- generate.conc(2, 1, 0:24)
  tmpconc$conc <- 0
  tmpdose <- generate.dose(tmpconc)
  myconc <- PKNCAconc(tmpconc, conc~time|treatment+ID)
  mydose <- PKNCAdose(tmpdose, dose~time|treatment+ID)
  mydata <- PKNCAdata(myconc, mydose)
  # Not capturing the warning due to R bug
  # https://bugs.r-project.org/bugzilla3/show_bug.cgi?id=17122
  #expect_warning(myresult <- pk.nca(mydata),
  #               regexp="Too few points for half-life calculation")
  myresult <- pk.nca(mydata)
  mysummary <- summary(myresult)
  expect_equal(mysummary,
               data.frame(start=0,
                          end=c(24, Inf),
                          treatment="Trt 1",
                          N="2",
                          auclast=c("NC", "."),
                          cmax=c(".", "NC"),
                          tmax=c(".", "NC"),
                          half.life=c(".", "NC"),
                          aucinf.obs=c(".", "NC"),
                          stringsAsFactors=FALSE),
               info="summary of PKNCAresults without most results gives NC")
  
  mysummary <- summary(myresult,
                       not.requested.string="NR",
                       not.calculated.string="NoCalc")
  expect_equal(mysummary,
               data.frame(start=0,
                          end=c(24, Inf),
                          treatment="Trt 1",
                          N="2",
                          auclast=c("NoCalc", "NR"),
                          cmax=c("NR", "NoCalc"),
                          tmax=c("NR", "NoCalc"),
                          half.life=c("NR", "NoCalc"),
                          aucinf.obs=c("NR", "NoCalc"),
                          stringsAsFactors=FALSE),
               info="Summary respects the not.requested.string and not.calculated.string")
  
  mysummary <- summary(myresult,
                       summarize.n.per.group=FALSE,
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
                          aucinf.obs=c("NR", "NoCalc"),
                          stringsAsFactors=FALSE),
               info="N is optionally omitted")
})
