test_that("pk.nca", {
  # Note that generate.conc sets the random seed, so it doesn't have to happen
  # here.
  tmpconc <- generate.conc(2, 1, 0:24)
  tmpdose <- generate.dose(tmpconc)
  myconc <- PKNCAconc(tmpconc, formula=conc~time|treatment+ID)
  mydose <- PKNCAdose(tmpdose, formula=dose~time|treatment+ID)
  mydata <- PKNCAdata(myconc, mydose)
  myresult <- pk.nca(mydata)

  expect_equal(names(myresult),
               c("result", "data", "columns"),
               info="Make sure that the result has the expected names (and only the expected names) in it.")
  expect_true(checkProvenance(myresult),
              info="Provenance works on results")

  mydata.failure <- mydata
  # There's no way to automatically make a PKNCAdata object with no intervals,
  # but we want to ensure that users cannot cause this error by playing in the
  # internals.
  mydata.failure$intervals <- data.frame()
  expect_warning(myresult.failure <- pk.nca(mydata.failure),
                 regexp="No intervals given; no calculations will be done.",
                 info="An empty result is returned if there are no intervals")

  tmpconc <- generate.conc(2, 1, 0:24)
  tmpdose <- generate.dose(tmpconc)
  myconc <- PKNCAconc(tmpconc, formula=conc~time|treatment+ID)
  mydose.nodose <- PKNCAdose(tmpdose, formula=~time|treatment+ID)
  mydata.nodose <- PKNCAdata(myconc, mydose.nodose)
  expect_equal(
    pk.nca(mydata.nodose)$result,
    myresult$result,
    info="missing dose information is handled without an issue"
  )

  # Test each of the pieces for myresult for accuracy

  expect_equal(myresult$data, {
    tmp <- mydata
    # The options should be the default options after the calculations are done.
    tmp$options <- PKNCA.options()
    tmp
  }, info="The data is just a copy of the input data plus an instantiation of the PKNCA.options")

  verify.result <-
    tibble::tibble(
      treatment="Trt 1",
      ID=rep(c(1, 2), each=15),
      start=0,
      end=c(24, rep(Inf, 14),
            24, rep(Inf, 14)),
      PPTESTCD=rep(c("auclast", "cmax", "tmax", "tlast", "clast.obs",
                     "lambda.z", "r.squared", "adj.r.squared",
                     "lambda.z.time.first", "lambda.z.time.last",
                     "lambda.z.n.points", "clast.pred", "half.life",
                     "span.ratio", "aucinf.obs"),
                   times=2),
      PPORRES=c(13.54, 0.9998, 4.000, 24.00, 0.3441,
                0.04297, 0.9072, 0.9021, 5.000, 24.00,
                20.00, 0.3356, 16.13, 1.178,
                21.55, 14.03, 0.9410, 2.000,
                24.00, 0.3148, 0.05689, 0.9000, 0.8944,
                5.000, 24.00, 20.00, 0.3011, 12.18,
                1.560, 19.56),
      exclude=NA_character_
    )
  expect_equal(
    myresult$result,
    verify.result,
    tolerance=0.001,
    info=paste("The specific order of the levels isn't important--",
               "the fact that they are factors and that the set",
               "doesn't change is important.")
  )

  # Specifying new intervals
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


  # Dosing not at time 0
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
               tolerance=0.001,
               info="Shifted dosing works the same as un-shifted where time parameters like tmax and tlast are reported relative to the start of the interval")

  tmpconc <- generate.conc(2, 1, 0:24)
  tmpdose <- generate.dose(tmpconc)
  myconc <- PKNCAconc(tmpconc, conc~time|treatment+ID)
  mydose <- PKNCAdose(tmpdose, dose~time|treatment+ID)
  mydata <- PKNCAdata(myconc, mydose,
                      intervals=data.frame(start=0, end=Inf, cmax=TRUE))
  myresult <- pk.nca(mydata)
  expect_equal(myresult$result$PPORRES,
               c(0.99981, 0.94097), tolerance=0.00001,
               info="Calculations work with a single row of intervals and a single parameter requested")

  mydata <- PKNCAdata(myconc, mydose,
                      intervals=data.frame(start=0, end=Inf, cl.obs=TRUE))
  myresult <- pk.nca(mydata)
  expect_equal(subset(myresult$result, PPTESTCD %in% "cl.obs")$PPORRES,
               c(0.04640, 0.05111), tolerance=0.0001,
               info="PK intervals work with passing in dose as a parameter")

  tmpconc <- generate.conc(2, 1, 0:24)
  tmpdose <- generate.dose(tmpconc)
  myconc <- PKNCAconc(tmpconc, conc~time|treatment+ID)
  mydose <- PKNCAdose(tmpdose, dose~time|treatment+ID)
  mydata <- PKNCAdata(myconc, mydose,
                      intervals=data.frame(start=0, end=24, cmax=TRUE, cav=TRUE))
  myresult <- pk.nca(mydata)
  expect_equal(subset(myresult$result, PPTESTCD %in% "cav")$PPORRES,
               c(0.5642, 0.5846), tolerance=0.0001,
               info="PK intervals work with passing in start and end as parameters")

  # Ensure that the correct number of doses are included in parameters that use dosing.
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
               tolerance=0.0001,
               info="The correct number of doses is selected for an interval (>=start and <end), 4 doses and not 5")

  mydata <- PKNCAdata(myconc, mydose,
                      intervals=data.frame(start=0, end=6, cl.last=TRUE))
  myresult <- pk.nca(mydata)
  expect_equal(myresult$result$PPORRES[myresult$result$PPTESTCD %in% "cl.last"],
               1/myresult$result$PPORRES[myresult$result$PPTESTCD %in% "auclast"],
               tolerance=0.0001,
               info="The correct number of doses is selected for an interval (>=start and <end), 1 dose and not 5")

  mydata <- PKNCAdata(myconc, mydose,
                      intervals=data.frame(start=1, end=6, cl.last=TRUE))
  myresult <- pk.nca(mydata)
  expect_equal(myresult$result$PPORRES[myresult$result$PPTESTCD %in% "cl.last"],
               NA/myresult$result$PPORRES[myresult$result$PPTESTCD %in% "auclast"],
               tolerance=0.0001,
               info="The correct number of doses is selected for an interval (>=start and <end), no doses selected")

})

test_that("verbose pk.nca", {
  tmpconc <- generate.conc(nsub=1, ntreat=1, 0:4)
  tmpdose <- generate.dose(tmpconc)
  myconc <- PKNCAconc(tmpconc, formula=conc~time)
  mydose <- PKNCAdose(tmpdose, formula=dose~time)
  mydata <- PKNCAdata(myconc, mydose)
  expect_message(expect_message(expect_message(
    suppressWarnings(pk.nca(mydata, verbose=TRUE)),
    regexp = "Setting up options"),
    regexp = "Starting dense PK NCA calculations"),
    regexp = "Combining completed dense PK calculation results"
  )
  expect_message(
    suppressWarnings(pk.nca(mydata, verbose=FALSE)),
    NA
  )
})

test_that("pk.nca warnings", {
  tmpconc <- generate.conc(nsub=1, ntreat=1, 0:4)
  tmpdose <- generate.dose(tmpconc)
  myconc <- PKNCAconc(tmpconc, formula=conc~time)
  mydose <- PKNCAdose(tmpdose, formula=dose~time)
  mydata <- PKNCAdata(myconc, mydose, intervals=data.frame(start=24, end=48, cmax=TRUE))
  expect_warning(
    pk.nca(mydata),
    regexp="No data for interval"
  )
})

test_that("pk.nca.interval errors", {
  expect_error(
    pk.nca.interval(interval="A"),
    regexp="Interval must be a data.frame"
  )
  expect_error(
    pk.nca.interval(interval=data.frame()),
    regexp="Interval must be a one-row data.frame"
  )
})

test_that("Calculations when dose time is missing", {
  # Ensure that the correct number of doses are included in parameters that use dosing.
  tmpconc <- generate.conc(2, 1, 0:24)
  tmpdose <- generate.dose(tmpconc)
  myconc <- PKNCAconc(tmpconc, conc~time|treatment+ID)
  mydose <- PKNCAdose(tmpdose, dose~.|treatment+ID)
  mydata <- PKNCAdata(myconc, mydose,
                      intervals=data.frame(start=0, end=24, cl.obs=TRUE))
  myresult <- pk.nca(mydata)
  expect_equal(
    myresult$result$PPORRES[myresult$result$PPTESTCD %in% "cl.obs"],
    1/myresult$result$PPORRES[myresult$result$PPTESTCD %in% "aucinf.obs"],
    info="The correct number of doses is selected for an interval (>=start and <end), 4 doses and not 5"
  )
})

test_that("Calculations when no dose info is given", {
  tmpconc <- generate.conc(2, 1, 0:24)
  myconc <- PKNCAconc(tmpconc, formula=conc~time|treatment+ID)
  mydata <- PKNCAdata(myconc, intervals=data.frame(start=0, end=24, cmax=TRUE, cl.last=TRUE))
  expect_message(
    myresult <- pk.nca(mydata),
    regexp="No dose information provided, calculations requiring dose will return NA.",
    info="Dosing information not required."
  )
  expect_equal(
    myresult$result,
    tibble::tibble(
      treatment="Trt 1",
      ID=rep(1:2, each=3),
      start=0,
      end=24,
      PPTESTCD=rep(c("auclast", "cmax", "cl.last"), 2),
      PPORRES=c(13.5417297156528, 0.999812606062292, NA,
                14.0305397438242, 0.94097296083447, NA),
      exclude=NA_character_
    )
  )
})

test_that("pk.nca with exclusions", {
  # Note that generate.conc sets the random seed, so it doesn't have to happen
  # here.
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
               tolerance=1e-5,
               info="duration.dose is used when requested")
})

test_that("pk.calc.all with duration.conc required", {
  tmpconc <- generate.conc(2, 1, 0:24)
  tmpdose <- generate.dose(tmpconc)
  tmpconc$duration_conc <- 0.1
  myconc <- PKNCAconc(tmpconc, formula=conc~time|treatment+ID, duration="duration_conc")
  mydose <- PKNCAdose(tmpdose, formula=dose~time|treatment+ID, route="intravascular")
  mydata <- PKNCAdata(myconc, mydose,
                      intervals=data.frame(start=0, end=24,
                                           mrt.iv.last=TRUE))
  myresult <- pk.nca(mydata)
  expect_equal(myresult$result$PPORRES[myresult$result$PPTESTCD %in% "mrt.iv.last"],
               c(10.41263, 10.17515),
               tolerance=1e-5,
               info="duration.conc is used when requested")
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

test_that("include_half.life and exclude_half.life work with NAs treated as missing for all NA and as FALSE for partial NA (#372)", {
  # Partial NA include_hl is used
  d_conc_incl <- data.frame(conc = c(1, 0.6, 0.3, 0.25, 0.15, 0.1), time = 0:5, include_hl = c(FALSE, NA, TRUE, TRUE, TRUE, TRUE))
  o_conc_incl <- PKNCAconc(d_conc_incl, conc~time, include_half.life = "include_hl")
  o_data_incl <- PKNCAdata(o_conc_incl, intervals = data.frame(start = 0, end = Inf, half.life = TRUE))
  suppressMessages(o_nca_incl <- pk.nca(o_data_incl))
  expect_equal(as.data.frame(o_nca_incl, out_format = "wide")$half.life, 1.820879, tolerance = 0.00001)

  # All FALSE include_hl is used
  d_conc_false <- data.frame(conc = c(1, 0.6, 0.3, 0.25, 0.15, 0.1), time = 0:5, include_hl = FALSE)
  o_conc_false <- PKNCAconc(d_conc_false, conc~time, include_half.life = "include_hl")
  o_data_false <- PKNCAdata(o_conc_false, intervals = data.frame(start = 0, end = Inf, half.life = TRUE))
  suppressWarnings(suppressMessages(o_nca_false <- pk.nca(o_data_false)))
  d_nca_false <- as.data.frame(o_nca_false)
  expect_equal(d_nca_false$PPORRES[d_nca_false$PPTESTCD %in% "half.life"], NA_real_)

  # All NA include_hl is ignored
  d_conc <- data.frame(conc = c(1, 0.6, 0.3, 0.25, 0.15, 0.1), time = 0:5, include_hl = NA)
  o_conc <- PKNCAconc(d_conc, conc~time, include_half.life = "include_hl")
  o_data <- PKNCAdata(o_conc, intervals = data.frame(start = 0, end = Inf, half.life = TRUE))
  suppressMessages(o_nca <- pk.nca(o_data))
  expect_equal(as.data.frame(o_nca, out_format = "wide")$half.life, 1.512942, tolerance = 0.00001)

  # Partial NA include_hl is used
  d_conc_excl <- data.frame(conc = c(1, 0.6, 0.3, 0.25, 0.15, 0.1), time = 0:5, exclude_hl = c(FALSE, NA, TRUE, TRUE, TRUE, TRUE))
  o_conc_excl <- PKNCAconc(d_conc_excl, conc~time, exclude_half.life = "exclude_hl")
  o_data_excl <- PKNCAdata(o_conc_excl, intervals = data.frame(start = 0, end = Inf, half.life = TRUE))
  suppressWarnings(suppressMessages(o_nca_excl <- pk.nca(o_data_excl)))
  d_nca_excl <- as.data.frame(o_nca_excl)
  expect_equal(d_nca_excl$PPORRES[d_nca_excl$PPTESTCD %in% "half.life"], NA_real_)

  # All NA exclude_hl is ignored
  d_conc <- data.frame(conc = c(1, 0.6, 0.3, 0.25, 0.15, 0.1), time = 0:5, exclude_hl = NA)
  o_conc <- PKNCAconc(d_conc, conc~time, exclude_half.life = "exclude_hl")
  o_data <- PKNCAdata(o_conc, intervals = data.frame(start = 0, end = Inf, half.life = TRUE))
  suppressMessages(o_nca <- pk.nca(o_data))
  expect_equal(as.data.frame(o_nca, out_format = "wide")$half.life, 1.512942, tolerance = 0.00001)

  # All FALSE exclude_hl is used
  d_conc_false <- data.frame(conc = c(1, 0.6, 0.3, 0.25, 0.15, 0.1), time = 0:5, exclude_hl = FALSE)
  o_conc_false <- PKNCAconc(d_conc_false, conc~time, exclude_half.life = "exclude_hl")
  o_data_false <- PKNCAdata(o_conc_false, intervals = data.frame(start = 0, end = Inf, half.life = TRUE))
  suppressWarnings(suppressMessages(o_nca_false <- pk.nca(o_data_false)))
  d_nca_false <- as.data.frame(o_nca_false)
  expect_equal(d_nca_false$PPORRES[d_nca_false$PPTESTCD %in% "half.life"], 1.512942, tolerance = 0.00001)
})

test_that("No interval requested (e.g. for placebo)", {
  tmpconc <- generate.conc(2, 1, 0:24)
  tmpdose <- generate.dose(tmpconc)
  myconc <- PKNCAconc(tmpconc, formula=conc~time|treatment+ID)
  mydose <- PKNCAdose(tmpdose, formula=dose~time|treatment+ID)
  mydata <-
    PKNCAdata(
      myconc, mydose,
      intervals=
        data.frame(
          treatment="Trt 3", start=0, end=24, cmax=TRUE,
          stringsAsFactors=FALSE
        )
    )
  expect_warning(expect_warning(expect_warning(expect_warning(
    myresult <- pk.nca(mydata),
    class = "pknca_no_intervals"),
    class = "pknca_no_intervals"),
    class = "pknca_no_conc_data"),
    class = "pknca_all_warnings_no_results"
  )
  expect_equal(
    nrow(as.data.frame(myresult)),
    0,
    info="No rows were generated when no intervals applied"
  )
})

test_that("Volume-related calculations", {
  tmpconc <- generate.conc(2, 1, c(4, 12, 24))
  tmpconc$conc <- seq_len(nrow(tmpconc))
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
  myresult <- pk.nca(mydata)
  expect_true(all(is.na(myresult$result[["PPORRES"]]) == c(FALSE, FALSE, FALSE, TRUE)) &
                all(myresult$result[["PPTESTCD"]] == rep(c("auclast", "cl.last"), 2)),
              info="cl.last is not calculated when dose information is missing, but only for the subject where dose info is missing.")
})

# Fix issue #68
test_that("Ensure that options are respected during pk.nca call", {
  doses <- data.frame(ID=1:2, Time=0, Dose=0.5)

  conc.data <- c(0, 1, 2, 1.3, 0.4, 0.35, 0.125)
  time.data <- c(0, 1, 2, 4,   8,   24,   48)
  concs <- merge(doses["ID"], data.frame(Conc=conc.data, Time=time.data))

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
  expect_true(
    all.equal(
      linear.results$result$PPORRES[linear.results$result$PPTESTCD %in% "aucinf.obs" &
                                      linear.results$result$ID %in% 1 &
                                      linear.results$result$end %in% Inf],
      24.54319,
      tolerance=0.0001) &
      all.equal(
        linlog.results$result$PPORRES[linlog.results$result$PPTESTCD %in% "aucinf.obs" &
                                        linlog.results$result$ID %in% 1 &
                                        linlog.results$result$end %in% Inf],
        23.68317,
        tolerance=0.0001),
    info="linear and loglinear effects are calculated differently."
  )
})

test_that("Can calculate parameters requiring extra arguments", {
  o_conc <- PKNCAconc(conc~time, data=data.frame(conc=c(1:3, 2:1), time=0:4))
  d_intervals <- data.frame(start=0, end=4, time_above=TRUE, conc_above=2)
  o_data <- PKNCAdata(o_conc, intervals=d_intervals, options=list(auc.method="linear"))
  o_nca <- suppressMessages(pk.nca(o_data))
  expect_equal(as.data.frame(o_nca)$PPORRES, 2)
})

test_that("calculate with sparse data", {
  d_sparse <-
    data.frame(
      id = c(1L, 2L, 3L, 1L, 2L, 3L, 1L, 2L, 3L, 4L, 5L, 6L, 4L, 5L, 6L, 7L, 8L, 9L, 7L, 8L, 9L),
      conc = c(0, 0, 0,  1.75, 2.2, 1.58, 4.63, 2.99, 1.52, 3.03, 1.98, 2.22, 3.34, 1.3, 1.22, 3.54, 2.84, 2.55, 0.3, 0.0421, 0.231),
      time = c(0, 0, 0, 1, 1, 1, 6, 6, 6, 2, 2, 2, 10, 10, 10, 4, 4, 4, 24, 24, 24),
      dose = c(100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100)
    )
  o_conc_sparse <- PKNCAconc(d_sparse, conc~time|id, sparse=TRUE)

  d_intervals <-
    data.frame(
      start=0,
      end=24,
      aucinf.obs=TRUE,
      cmax=TRUE,
      sparse_auclast=TRUE
    )
  o_data_sparse <- PKNCAdata(o_conc_sparse, intervals=d_intervals)
  suppressMessages(
    expect_warning(expect_warning(
      o_nca <- pk.nca(o_data_sparse),
      class = "pknca_sparse_df_multi"),
      class = "pknca_halflife_too_few_points"
    )
  )
  df_result <- as.data.frame(o_nca)
  expect_true("sparse_auclast" %in% df_result$PPTESTCD)
  expect_equal(df_result$PPORRES[df_result$PPTESTCD %in% "sparse_auclast"], 39.4689)
  sum_o_nca <- summary(o_nca)
  expect_s3_class(sum_o_nca, "summary_PKNCAresults")

  # Mixed sparse and dense calculations when only one type is requested in an
  # interval works The example below has dense-only; sparse and dense; and and
  # sparse-only.
  d_intervals_mixed <-
    data.frame(
      start=0,
      end=c(23, 24, 25),
      cmax=c(TRUE, TRUE, FALSE),
      sparse_auclast=c(FALSE, TRUE, TRUE)
    )
  o_data_sparse_mixed <- PKNCAdata(o_conc_sparse, intervals=d_intervals_mixed)
  suppressMessages(
    expect_warning(expect_warning(
      o_nca_sparse_mixed <- pk.nca(o_data_sparse_mixed),
      class = "pknca_sparse_df_multi"),
      class = "pknca_sparse_df_multi"
    )
  )
  df_result_sparse_mixed <- as.data.frame(o_nca_sparse_mixed)
  expect_true("sparse_auclast" %in% df_result_sparse_mixed$PPTESTCD)
  expect_equal(df_result_sparse_mixed$PPORRES[df_result_sparse_mixed$PPTESTCD %in% "sparse_auclast"], rep(39.4689, 2))
  suppressMessages(
    expect_message(
      expect_warning(expect_warning(
        o_nca_sparse_mixed <- pk.nca(o_data_sparse_mixed, verbose=TRUE),
        class = "pknca_sparse_df_multi"),
        class = "pknca_sparse_df_multi"
      ),
      regexp="No sparse calculations requested for an interval"
    )
  )

  # Sparse data with multiple treatments, confirm the correct number of rows of
  # outputs are created.
  d_sparse_200 <- d_sparse
  d_sparse_200$dose <- 200
  d_sparse_multi_trt <- rbind(d_sparse, d_sparse_200)
  d_sparse_multi_trt$dose_grp <- d_sparse_multi_trt$dose
  o_conc_sparse_multi_trt <- PKNCAconc(d_sparse_multi_trt, conc~time|dose_grp+id, sparse=TRUE)
  d_intervals_mixed <-
    data.frame(
      start=0,
      end=c(23, 24, 25),
      cmax=c(TRUE, TRUE, FALSE),
      sparse_auclast=c(FALSE, TRUE, TRUE)
    )
  d_dose_sparse_multi_trt <- unique(d_sparse_multi_trt[, c("id", "dose")])
  d_dose_sparse_multi_trt$time <- 0
  d_dose_sparse_multi_trt$dose_grp <- d_dose_sparse_multi_trt$dose
  o_dose_sparse_multi_trt <- PKNCAdose(d_dose_sparse_multi_trt, dose~time|dose_grp+id)
  o_data_sparse_multi_trt <- PKNCAdata(o_conc_sparse_multi_trt, o_dose_sparse_multi_trt, intervals=d_intervals_mixed)
  suppressMessages(
    expect_warning(expect_warning(expect_warning(expect_warning(
      o_nca_sparse_multi_trt <- pk.nca(o_data_sparse_multi_trt),
      class = "pknca_sparse_df_multi"),
      class = "pknca_sparse_df_multi"),
      class = "pknca_sparse_df_multi"),
      class = "pknca_sparse_df_multi"
    )
  )
  expect_equal(nrow(as.data.frame(o_nca_sparse_multi_trt)), 16)

  # Correct detection of mixed doses within a sparse dose group when there are
  # no grouping columns other than subject
  d_sparse_multi_trt_bad_dose_single <- d_sparse_multi_trt[d_sparse_multi_trt$dose == 100, ]
  d_dose_sparse_multi_trt_bad_dose_single <- unique(d_sparse_multi_trt_bad_dose_single[, c("id", "dose")])
  d_dose_sparse_multi_trt_bad_dose_single$time <- 0
  d_dose_sparse_multi_trt_bad_dose_single$dose[1] <- d_dose_sparse_multi_trt_bad_dose_single$dose[1] + 1
  o_conc_sparse_multi_trt_bad_dose_single <- PKNCAconc(d_sparse_multi_trt_bad_dose_single, conc~time|id, sparse=TRUE)
  o_dose_sparse_multi_trt_bad_dose_single <- PKNCAdose(d_dose_sparse_multi_trt_bad_dose_single, dose~time|id)
  o_data_sparse_multi_trt_bad_dose_single <- PKNCAdata(o_conc_sparse_multi_trt_bad_dose_single, o_dose_sparse_multi_trt_bad_dose_single, intervals=d_intervals_mixed)
  expect_error(
    pk.nca(o_data_sparse_multi_trt_bad_dose_single),
    regexp="With sparse PK, all subjects in a group must have the same dosing information.*Not all subjects have the same dosing information"
  )

  # Correct detection of mixed doses within a sparse dose group
  d_dose_sparse_multi_trt_bad_dose <- d_dose_sparse_multi_trt
  d_dose_sparse_multi_trt_bad_dose$dose[1] <- d_dose_sparse_multi_trt$dose[1] + 1
  o_dose_sparse_multi_trt_bad_dose <- PKNCAdose(d_dose_sparse_multi_trt_bad_dose, dose~time|dose_grp+id)
  o_data_sparse_multi_trt_bad_dose <- PKNCAdata(o_conc_sparse_multi_trt, o_dose_sparse_multi_trt_bad_dose, intervals=d_intervals_mixed)
  expect_error(
    pk.nca(o_data_sparse_multi_trt_bad_dose),
    regexp="With sparse PK, all subjects in a group must have the same dosing information.*Not all subjects have the same dosing information for this group: +dose_grp=100"
  )
  # Correct detection of mixed doses within a sparse dose group when there are no groups
})

test_that("Unexpected interval columns now not cause an error (#238)", {
  d_conc <-
    data.frame(
      ID = 1L,
      time = 0:6,
      conc = c(0, 0.7, 0.71, 0.85, 1, 0.76, 0.74)
    )
  d_dose <- data.frame(dose = 1)
  d_intervals <- data.frame(start = 0, end = 6, cmax = TRUE, aucinf = TRUE)
  o_conc <- PKNCAconc(d_conc, formula = conc~time|ID)
  o_dose <- PKNCAdose(d_dose, formula = dose~.)
  expect_error(PKNCAdata(o_conc, o_dose, intervals = d_intervals),
               "The following columns in 'intervals' are not allowed:"
  )
})

test_that("aucint works within pk.calc.all for all zero concentrations with interpolated or extrapolated concentrations", {
  # AUCint.inf.obs
  d_interval <- data.frame(start = 0, end = 4, aucint.inf.obs = TRUE)
  d_conctime <- data.frame(conc = c(0, 0, 0, 0), time = 0:3)
  o_conc <- PKNCAconc(d_conctime, conc~time)
  o_data <- PKNCAdata(o_conc, intervals = d_interval)
  suppressWarnings(suppressMessages(
    o_nca <- pk.nca(o_data)
  ))
  results <- setNames(as.data.frame(o_nca)$PPORRES, nm = as.data.frame(o_nca)$PPTESTCD)
  zero_names <- c("clast.obs", "aucint.inf.obs")
  na_names <- setdiff(names(results), zero_names)
  expect_equal(
    results[zero_names],
    setNames(rep(0, length(zero_names)), zero_names)
  )
  expect_equal(
    results[na_names],
    setNames(rep(NA_real_, length(na_names)), na_names)
  )

  # AUCint.inf.pred
  d_interval <- data.frame(start = 0, end = 4, aucint.inf.pred = TRUE)
  d_conctime <- data.frame(conc = c(0, 0, 0, 0), time = 0:3)
  o_conc <- PKNCAconc(d_conctime, conc~time)
  o_data <- PKNCAdata(o_conc, intervals = d_interval)
  suppressWarnings(suppressMessages(
    o_nca <- pk.nca(o_data)
  ))
  expect_equal(
    as.data.frame(o_nca)$PPORRES,
    c(rep(NA_real_, 11), 0)
  )
})

test_that("The option keep_interval_cols is respected", {
  d_interval <- data.frame(start = 0, end = 4, cmax = TRUE, foo = "A")
  d_conctime <- data.frame(conc = c(0, 0, 0, 0), time = 0:3)
  o_conc <- PKNCAconc(d_conctime, conc~time)
 expect_error(PKNCAdata(o_conc, intervals = d_interval),
              "The following columns in 'intervals' are not allowed:")

  o_data <- PKNCAdata(o_conc, intervals = d_interval, options = list(keep_interval_cols = "foo"))
  suppressWarnings(suppressMessages(
    o_nca <- pk.nca(o_data)
  ))
  expect_equal(o_nca$result$foo, "A")
  expect_true("foo" %in% names(summary(o_nca)))
})

test_that("dose is calculable", {
  tmpconc <- generate.conc(2, 1, 0:24)
  tmpdose <- generate.dose(tmpconc)
  myconc <- PKNCAconc(tmpconc, formula=conc~time|treatment+ID)
  mydose <- PKNCAdose(tmpdose, formula=dose~time|treatment+ID)
  mydata <- PKNCAdata(myconc, mydose, intervals = data.frame(start = 0, end = Inf, totdose = TRUE))
  myresult <- pk.nca(mydata)

  # One dose in the interval
  expect_equal(as.data.frame(myresult)$PPORRES, rep(1, 2))
  expect_equal(as.data.frame(myresult)$PPTESTCD, rep("totdose", 2))

  # Don't give dose data
  mydata <- PKNCAdata(myconc, intervals = data.frame(start = 0, end = Inf, totdose = TRUE))
  suppressMessages(myresult <- pk.nca(mydata))
  expect_equal(as.data.frame(myresult)$PPORRES, rep(NA_real_, 2))
  expect_equal(as.data.frame(myresult)$PPTESTCD, rep("totdose", 2))

  # Multiple doses in the interval
  tmpdose_second <- tmpdose
  tmpdose_second$time <- 1
  mydose <- PKNCAdose(rbind(tmpdose, tmpdose_second), formula=dose~time|treatment+ID)
  mydata <- PKNCAdata(myconc, mydose, intervals = data.frame(start = 0, end = Inf, totdose = TRUE))
  myresult <- pk.nca(mydata)
  expect_equal(as.data.frame(myresult)$PPORRES, rep(2, 2))
  expect_equal(as.data.frame(myresult)$PPTESTCD, rep("totdose", 2))
})

test_that("do not give rbind error when interval columns have attributes (#381)", {
  o_conc <- PKNCAconc(data = data.frame(conc = 1, time = 0), conc~time)
  d_interval <- data.frame(start = 0, end = Inf, cmax = TRUE)

  d_interval <- data.frame(start = 0, end = Inf, cmax = TRUE, tmax = TRUE)
  attr(d_interval$start, "label") <- "start"
  o_data <- PKNCAdata(o_conc, intervals = d_interval)
  suppressMessages(o_nca <- pk.nca(o_data))
  # interval attributes are preserved
  expect_equal(
    attributes(as.data.frame(o_nca)$start),
    list(label = "start")
  )
})
