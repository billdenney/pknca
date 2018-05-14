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
                 regexp="No dose information provided, calculations requiring dose will return NA.",
                 info="Dosing information not required.")
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

test_that("pk.calc.all with duration.dose required", {
  tmpconc <- generate.conc(2, 1, 0:24)
  tmpdose <- generate.dose(tmpconc)
  tmpdose$duration_dose <- 0.1
  myconc <- PKNCAconc(tmpconc, formula=conc~time|treatment+ID)
  mydose <- PKNCAdose(tmpdose, formula=dose~time|treatment+ID, duration="duration_dose", route="intravascular")
  mydata <- PKNCAdata(myconc, mydose,
                      intervals=data.frame(start=0, end=24,
                                           mrt.iv.last=TRUE))
  myresult <- pk.nca(mydata)
  expect_equal(myresult$result$PPORRES[myresult$result$PPTESTCD %in% "mrt.iv.last"],
               c(10.36263, 10.12515),
               tol=1e-5,
               info="duration.dose is used when requested")
})

test_that("half life inclusion and exclusion", {
  tmpconc <- generate.conc(2, 1, 0:24)
  tmpdose <- generate.dose(tmpconc)
  tmpconc$include_hl <- tmpconc$time <= 22
  tmpconc$exclude_hl <- tmpconc$time == 22
  myconc <- PKNCAconc(tmpconc, formula=conc~time|treatment+ID)
  myconc_incl <- PKNCAconc(tmpconc, formula=conc~time|treatment+ID,
                           include_half.life="include_hl")
  myconc_excl <- PKNCAconc(tmpconc, formula=conc~time|treatment+ID,
                           exclude_half.life="exclude_hl")
  mydose <- PKNCAdose(tmpdose, formula=dose~time|treatment+ID)
  mydata <-  PKNCAdata(myconc, mydose,
                       intervals=data.frame(start=0, end=24, half.life=TRUE))
  mydata_incl <- PKNCAdata(myconc_incl, mydose,
                           intervals=data.frame(start=0, end=24, half.life=TRUE))
  mydata_excl <- PKNCAdata(myconc_excl, mydose,
                           intervals=data.frame(start=0, end=24, half.life=TRUE))
  myresult <- pk.nca(mydata)
  myresult_incl <- pk.nca(mydata_incl)
  myresult_excl <- pk.nca(mydata_excl)
  expect_false(identical(myresult$result, myresult_excl$result))
  expect_false(identical(myresult$result, myresult_incl$result))
})

test_that("No interval requested (e.g. for placebo)", {
  tmpconc <- generate.conc(2, 1, 0:24)
  tmpdose <- generate.dose(tmpconc)
  myconc <- PKNCAconc(tmpconc, formula=conc~time|treatment+ID)
  mydose <- PKNCAdose(tmpdose, formula=dose~time|treatment+ID)
  mydata <-  PKNCAdata(myconc, mydose,
                       intervals=data.frame(treatment="Trt 3", start=0, end=24, cmax=TRUE,
                                            stringsAsFactors=FALSE))
  expect_message(myresult <- pk.nca(mydata),
                 regexp="2 groups have no interval calculations requested.",
                 info="No intervals apply to a group provides a message.")
  expect_equal(nrow(as.data.frame(myresult)), 0,
               info="No rows were generated when no intervals applied")
})

test_that("Volume-related calculations", {
  tmpconc <- generate.conc(2, 1, c(4, 12, 24))
  tmpconc$conc <- 1:nrow(tmpconc)
  tmpconc$vol <- 2
  tmpdose <- generate.dose(tmpconc)
  myconc <- PKNCAconc(tmpconc, formula=conc~time|treatment+ID, volume="vol")
  mydose <- PKNCAdose(tmpdose, formula=dose~time|treatment+ID)
  mydata <-  PKNCAdata(myconc, mydose,
                       intervals=data.frame(treatment="Trt 1", start=0, end=24,
                                            ae=TRUE, fe=TRUE,
                                            stringsAsFactors=FALSE))
  myresult <- pk.nca(mydata)
  expect_equal(as.data.frame(myresult)[["PPORRES"]], c(12, 12, 30, 30),
              info="ae and fe are correctly calculated")
  tmpdose2 <- tmpdose
  tmpdose2$dose <- 2
  mydose2 <- PKNCAdose(tmpdose2, formula=dose~time|treatment+ID)
  mydata2 <-  PKNCAdata(myconc, mydose2,
                       intervals=data.frame(treatment="Trt 1", start=0, end=24,
                                            ae=TRUE, fe=TRUE,
                                            stringsAsFactors=FALSE))
  myresult2 <- pk.nca(mydata2)
  expect_equal(as.data.frame(myresult2)[["PPORRES"]], c(12, 6, 30, 15),
               info="fe respects dose")
})

test_that("pk.nca can calculate values with group-level data", {
  tmpconc_impute <- generate.conc(2, 1, 0:24)
  # This is what will happen in the imputation
  tmpconc_observe_05 <- tmpconc_impute[tmpconc_impute$time %in% 0,]
  tmpconc_observe_05$time <- 0.5
  tmpconc_observe <- rbind(tmpconc_impute, tmpconc_observe_05)
  tmpconc_observe <- tmpconc_observe[order(tmpconc_observe$treatment, tmpconc_observe$ID, tmpconc_observe$time),]
  tmpdose <- generate.dose(tmpconc_impute)
  tmpdose$time <- 0.5
  
  myconc_impute <- PKNCAconc(tmpconc_impute, formula=conc~time|treatment+ID)
  myconc_observe <- PKNCAconc(tmpconc_observe, formula=conc~time|treatment+ID)
  mydose <- PKNCAdose(tmpdose, formula=dose~time|treatment+ID)
  mydata_impute <-
    PKNCAdata(myconc_impute, mydose,
              intervals=data.frame(treatment="Trt 1", start=0, end=24,
                                   aucint.last.dose=TRUE,
                                   stringsAsFactors=FALSE))
  mydata_observe <-
    PKNCAdata(myconc_observe, mydose,
              intervals=data.frame(treatment="Trt 1", start=0, end=24,
                                   auclast=TRUE,
                                   stringsAsFactors=FALSE))
  myres_impute <- pk.nca(mydata_impute)
  myres_observe <- pk.nca(mydata_observe)
  expect_equal(as.data.frame(myres_impute)$PPORRES,
               as.data.frame(myres_observe)$PPORRES,
               info="Manually imputing values gives the same result as aucint")
})

test_that("Missing dose info for some subjects gives a warning, not a difficult-to-interpret error", {
  tmpconc <- generate.conc(2, 1, 0:24)
  tmpdose <- generate.dose(tmpconc)[1,]
  myconc <- PKNCAconc(tmpconc, formula=conc~time|treatment+ID)
  mydose <- PKNCAdose(tmpdose, formula=dose~time|treatment+ID)
  mydata <-  PKNCAdata(myconc, mydose,
                       intervals=data.frame(start=0, end=24,
                                            cl.last=TRUE))
  expect_warning(myresult <- pk.nca(mydata),
                 regexp="The following intervals are missing dosing data:",
                 fixed=TRUE,
                 info="Warning is issued when dose data are missing for an interval but provided for some data.")
  expect_true(all(is.na(myresult$result[["PPORRES"]]) == c(FALSE, FALSE, FALSE, TRUE)) &
                all(myresult$result[["PPTESTCD"]] == rep(c("auclast", "cl.last"), 2)),
              info="cl.last is not calculated when dose information is missing, but only for the subject where dose info is missing.")
})

# Fix issue #68
test_that("Ensure that options are respected during pk.nca call", {
  doses <- data.frame(ID=1:2, Time=0, Dose=0.5)
  
  conc.data <- c(0, 1, 2, 1.3, 0.4, 0.35, 0.125)
  time.data <- c(0, 1, 2, 4,   8,   24,   48)
  concs <- merge(doses[c("ID")], data.frame(Conc=conc.data, Time=time.data))
  
  myconc <- PKNCA::PKNCAconc(concs, formula=Conc~Time|ID)
  mydose <- PKNCA::PKNCAdose(doses, formula=Dose~Time|ID)
  
  myintervals <- data.frame(start=c(0,0,0),
                            end=c(24,48,Inf),
                            auclast=TRUE,
                            aucinf.obs=TRUE,
                            aucinf.pred=TRUE,
                            aumclast=TRUE,
                            aumcall=TRUE,
                            half.life=TRUE)
  
  linear.mydata <- PKNCA::PKNCAdata(myconc, mydose, intervals = myintervals,
                                    options = list(auc.method = "linear"))
  linear.results <- PKNCA::pk.nca(linear.mydata)
  
  linlog.mydata <- PKNCA::PKNCAdata(myconc, mydose, intervals = myintervals,
                                    options = list(auc.method = "lin up/log down"))
  linlog.results <- PKNCA::pk.nca(linlog.mydata)
  expect_true(all.equal(linear.results$result$PPORRES[linear.results$result$PPTESTCD %in% "aucinf.obs" &
                                                         linear.results$result$ID %in% 1 &
                                                         linear.results$result$end %in% Inf],
                         24.54319,
                         tolerance=0.0001) &
              all.equal(linlog.results$result$PPORRES[linlog.results$result$PPTESTCD %in% "aucinf.obs" &
                                                        linlog.results$result$ID %in% 1 &
                                                        linlog.results$result$end %in% Inf],
                        23.68317,
                        tolerance=0.0001),
              info="linear and loglinear effects are calculated differently.")
})
