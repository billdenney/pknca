context("All NCA calculations")

library(dplyr)
source("generate.data.R")

test_that("pk.nca", {
  ## Note that generate.conc sets the random seed, so it doesn't have
  ## to happen here.
  tmpconc <- generate.conc(2, 1, 0:24)
  tmpdose <- generate.dose(tmpconc)
  myconc <- PKNCAconc(tmpconc, formula=conc~time|treatment+ID)
  mydose <- PKNCAdose(tmpdose, formula=dose~time|treatment+ID)
  mydata <- PKNCAdata(myconc, mydose)
  myresult <- pk.nca(mydata)

  expect_equal(names(myresult),
               c("result", "data", "exclude"),
               info="Make sure that the result has the expected names (and only the expected names) in it.")
  expect_true(checkProvenance(myresult),
              info="Provenance works on results")

  mydata.failure <- mydata
  ## There's no way to automatically make a PKNCAdata object with no
  ## intervals, but we want to ensure that users cannot cause this error
  ## by playing in the internals.
  mydata.failure$intervals <- data.frame()
  expect_warning(myresult.failure <- pk.nca(mydata.failure),
                 regexp="No intervals given; no calculations done.",
                 info="An empty result is returned if there are no intervals")
  
  tmpconc <- generate.conc(2, 1, 0:24)
  tmpdose <- generate.dose(tmpconc)
  myconc <- PKNCAconc(tmpconc, formula=conc~time|treatment+ID)
  mydose.nodose <- PKNCAdose(tmpdose, formula=~time|treatment+ID)
  mydata.nodose <- PKNCAdata(myconc, mydose.nodose)
  expect_equal(pk.nca(mydata.nodose)$result,
               myresult$result,
               info="missing dose information is handled without an issue")
  
  
  ## Test each of the pieces for myresult for accuracy

  expect_equal(myresult$data, {
    tmp <- mydata
    ## The options should be the default options after the
    ## calculations are done.
    tmp$options <- PKNCA.options()
    tmp
  }, info="The data is just a copy of the input data plus an instantiation of the PKNCA.options")

  verify.result <-
    data.frame(
      start=0,
      end=c(24, rep(Inf, 13),
            24, rep(Inf, 13)),
      treatment="Trt 1",
      ID=rep(c(1, 2), each=14),
      PPTESTCD=rep(c("auclast", "cmax", "tmax", "tlast", "clast.obs",
                     "lambda.z", "r.squared", "adj.r.squared",
                     "lambda.z.time.first", "lambda.z.n.points",
                     "clast.pred", "half.life", "span.ratio",
                     "aucinf.obs"),
                   times=2),
      PPORRES=c(13.54, 0.9998, 4.000, 24.00, 0.3441,
                0.04297, 0.9072, 0.9021, 5.000,
                20.00, 0.3356, 16.13, 1.178,
                21.55, 14.03, 0.9410, 2.000,
                24.00, 0.3148, 0.05689, 0.9000, 0.8944,
                5.000, 20.00, 0.3011, 12.18,
                1.560, 19.56),
      exclude=NA_character_,
      stringsAsFactors=FALSE)
  expect_equal(myresult$result, verify.result,
               tol=0.001,
               info=paste("The specific order of the levels isn't important--",
                          "the fact that they are factors and that the set",
                          "doesn't change is important."))

  ## Specifying new intervals
  mydata.newinterval <-
      PKNCAdata(myconc, mydose,
                intervals=data.frame(start=0, end=c(24, Inf),
                                     auclast=c(TRUE, FALSE),
                                     aucinf.obs=c(FALSE, TRUE),
                                     cmax=c(FALSE, TRUE),
                                     tmax=c(FALSE, TRUE),
                                     half.life=c(FALSE, TRUE)))
  myresult.newinterval <- pk.nca(mydata)
  expect_equal(myresult.newinterval$result,
               myresult$result,
               info="Intervals can be specified manually, and will apply across appropriate parts of the grouping variables.")
  
  
  ## Dosing not at time 0
  tmpconc.multi <- generate.conc(2, 1, 0:24)
  tmpdose.multi <- generate.dose(tmpconc.multi)
  tmpconc.multi$time <- tmpconc.multi$time + 2
  tmpdose.multi$time <- tmpdose.multi$time + 2
  myconc.multi <- PKNCAconc(tmpconc.multi, conc~time|treatment+ID)
  mydose.multi <- PKNCAdose(tmpdose.multi, dose~time|treatment+ID)
  mydata.multi <- PKNCAdata(myconc.multi, mydose.multi)
  myresult.multi <- pk.nca(mydata.multi)

  verify.result.multi <- verify.result
  verify.result.multi$start <- verify.result.multi$start + 2
  verify.result.multi$end <- verify.result.multi$end + 2
  expect_equal(myresult.multi$result, verify.result.multi,
               tol=0.001,
               info="Shifted dosing works the same as un-shifted where time parameters like tmax and tlast are reported relative to the start of the interval")

  tmpconc <- generate.conc(2, 1, 0:24)
  tmpdose <- generate.dose(tmpconc)
  myconc <- PKNCAconc(tmpconc, conc~time|treatment+ID)
  mydose <- PKNCAdose(tmpdose, dose~time|treatment+ID)
  mydata <- PKNCAdata(myconc, mydose,
                      intervals=data.frame(start=0, end=Inf, cmax=TRUE))
  myresult <- pk.nca(mydata)
  expect_equal(myresult$result$PPORRES,
               c(0.99981, 0.94097), tol=0.00001,
               info="Calculations work with a single row of intervals and a single parameter requested")

  mydata <- PKNCAdata(myconc, mydose,
                      intervals=data.frame(start=0, end=Inf, cl.obs=TRUE))
  myresult <- pk.nca(mydata)
  expect_equal(subset(myresult$result, PPTESTCD %in% "cl.obs")$PPORRES,
               c(0.04640, 0.05111), tol=0.0001,
               info="PK intervals work with passing in dose as a parameter")

  tmpconc <- generate.conc(2, 1, 0:24)
  tmpdose <- generate.dose(tmpconc)
  myconc <- PKNCAconc(tmpconc, conc~time|treatment+ID)
  mydose <- PKNCAdose(tmpdose, dose~time|treatment+ID)
  mydata <- PKNCAdata(myconc, mydose,
                      intervals=data.frame(start=0, end=24, cmax=TRUE, cav=TRUE))
  myresult <- pk.nca(mydata)
  expect_equal(subset(myresult$result, PPTESTCD %in% "cav")$PPORRES,
               c(0.5642, 0.5846), tol=0.0001,
               info="PK intervals work with passing in start and end as parameters")
  
  ## Ensure that the correct number of doses are included in parameters that use dosing.
  tmpconc <- generate.conc(2, 1, 0:24)
  tmpdose <- generate.dose(tmpconc)
  tmpdose$time <- NULL
  tmpdose <- merge(tmpdose, data.frame(time=c(0, 6, 12, 18, 24)))
  myconc <- PKNCAconc(tmpconc, conc~time|treatment+ID)
  mydose <- PKNCAdose(tmpdose, dose~time|treatment+ID)
  mydata <- PKNCAdata(myconc, mydose,
                      intervals=data.frame(start=0, end=24, cl.obs=TRUE))
  myresult <- pk.nca(mydata)
  expect_equal(myresult$result$PPORRES[myresult$result$PPTESTCD %in% "cl.obs"],
               4/myresult$result$PPORRES[myresult$result$PPTESTCD %in% "aucinf.obs"],
               tol=0.0001,
               info="The correct number of doses is selected for an interval (>=start and <end), 4 doses and not 5")

  mydata <- PKNCAdata(myconc, mydose,
                      intervals=data.frame(start=0, end=6, cl.last=TRUE))
  myresult <- pk.nca(mydata)
  expect_equal(myresult$result$PPORRES[myresult$result$PPTESTCD %in% "cl.last"],
               1/myresult$result$PPORRES[myresult$result$PPTESTCD %in% "auclast"],
               tol=0.0001,
               info="The correct number of doses is selected for an interval (>=start and <end), 1 dose and not 5")

  mydata <- PKNCAdata(myconc, mydose,
                      intervals=data.frame(start=1, end=6, cl.last=TRUE))
  myresult <- pk.nca(mydata)
  expect_equal(myresult$result$PPORRES[myresult$result$PPTESTCD %in% "cl.last"],
               NA/myresult$result$PPORRES[myresult$result$PPTESTCD %in% "auclast"],
               tol=0.0001,
               info="The correct number of doses is selected for an interval (>=start and <end), no doses selected")
  
})

test_that("Calculations when dose time is missing", {
  ## Ensure that the correct number of doses are included in parameters that use dosing.
  tmpconc <- generate.conc(2, 1, 0:24)
  tmpdose <- generate.dose(tmpconc)
  myconc <- PKNCAconc(tmpconc, conc~time|treatment+ID)
  mydose <- PKNCAdose(tmpdose, dose~.|treatment+ID)
  mydata <- PKNCAdata(myconc, mydose,
                      intervals=data.frame(start=0, end=24, cl.obs=TRUE))
  myresult <- pk.nca(mydata)
  expect_equal(myresult$result$PPORRES[myresult$result$PPTESTCD %in% "cl.obs"],
               1/myresult$result$PPORRES[myresult$result$PPTESTCD %in% "aucinf.obs"],
               info="The correct number of doses is selected for an interval (>=start and <end), 4 doses and not 5")

  tmpdose$time <- NULL
  tmpdose <- merge(tmpdose, data.frame(time=c(0, 6, 12, 18, 24)))
})

test_that("Calculations when no dose info is given", {
  tmpconc <- generate.conc(2, 1, 0:24)
  myconc <- PKNCAconc(tmpconc, formula=conc~time|treatment+ID)
  mydata <- PKNCAdata(myconc, intervals=data.frame(start=0, end=24, cmax=TRUE, cl.last=TRUE))
  expect_message(myresult <- pk.nca(mydata),
                 regexp="No dose information provided, assuming default dosing information.",
                 info="Default dosing information is assumed if no dosing information is given.")
  expect_equal(myresult$result,
               data.frame(start=0,
                          end=24,
                          treatment="Trt 1",
                          ID=rep(1:2, each=3),
                          PPTESTCD=rep(c("auclast", "cmax", "cl.last"), 2),
                          PPORRES=c(13.5417297156528, 0.999812606062292, NA,
                                    14.0305397438242, 0.94097296083447, NA),
                          exclude=NA_character_,
                          stringsAsFactors=FALSE))
})

test_that("pk.nca with exclusions", {
  ## Note that generate.conc sets the random seed, so it doesn't have
  ## to happen here.
  tmpconc <- generate.conc(2, 1, 0:24)
  tmpdose <- generate.dose(tmpconc)
  myconc <- PKNCAconc(tmpconc, formula=conc~time|treatment+ID)
  mydose <- PKNCAdose(tmpdose, formula=dose~time|treatment+ID)
  mydata <- PKNCAdata(myconc, mydose)
  myresult <- pk.nca(mydata)
  tmpconc.excl <- tmpconc
  tmpconc.excl$excl <- NA_character_
  tmpconc.excl$excl[5] <- "test exclusion"
  myconc.excl <- PKNCAconc(tmpconc.excl,
                           formula=conc~time|treatment+ID,
                           exclude="excl")
  mydata.excl <- PKNCAdata(myconc.excl, mydose)
  myresult.excl <- pk.nca(mydata.excl)
  expect_true(identical(myresult$result[myresult$result$ID %in% 2,],
                        myresult.excl$result[myresult.excl$result$ID %in% 2,]),
               info="Results are unchanged for the subject who has the same data")
  expect_false(identical(myresult$result[myresult$result$ID %in% 1,],
                         myresult.excl$result[myresult.excl$result$ID %in% 1,]),
               info="Results are changed for the subject who has the same data")
  expect_equal(myresult.excl$result$PPORRES[myresult.excl$result$ID %in% 1 & myresult.excl$result$PPTESTCD %in% "cmax"],
               max(tmpconc.excl$conc[tmpconc.excl$ID %in% 1 & is.na(tmpconc.excl$excl)]),
               info="Cmax is affected by the exclusion")
})