test_that("PKNCAresults summary", {
  # Note that generate.conc sets the random seed, so it doesn't have
  # to happen here.
  tmpconc <- generate.conc(2, 1, 0:24)
  tmpdose <- generate.dose(tmpconc)
  myconc <- PKNCAconc(tmpconc, formula=conc~time|treatment+ID)
  mydose <- PKNCAdose(tmpdose, formula=dose~time|treatment+ID)
  mydata <- PKNCAdata(myconc, mydose)
  myresult <- pk.nca(mydata)

  # Testing the summarization
  mysummary <- summary(myresult)
  expect_true(is.data.frame(mysummary))
  expect_equal(
    mysummary,
    as_summary_PKNCAresults(
      data.frame(
        start=0,
        end=c(24, Inf),
        treatment="Trt 1",
        N="2",
        auclast=c("13.8 [2.51]", "."),
        cmax=c(".", "0.970 [4.29]"),
        tmax=c(".", "3.00 [2.00, 4.00]"),
        half.life=c(".", "14.2 [2.79]"),
        aucinf.obs=c(".", "20.5 [6.84]"),
        stringsAsFactors=FALSE
      ),
      caption="auclast, cmax, aucinf.obs: geometric mean and geometric coefficient of variation; tmax: median and range; half.life: arithmetic mean and standard deviation"
    ),
    info="simple summary of PKNCAresults performs as expected"
  )

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
  suppressWarnings(myresult <- pk.nca(mydata))
  mysummary <- summary(myresult)
  expect_equal(
    mysummary,
    as_summary_PKNCAresults(
      data.frame(
        start=0,
        end=c(24, Inf),
        treatment="Trt 1",
        N="2",
        auclast=c("13.5 [NC]", "."),
        cmax=c(".", "1.00 [NC]"),
        tmax=c(".", "4.00 [4.00, 4.00]"),
        half.life=c(".", "16.1 [NC]"),
        aucinf.obs=c(".", "21.5 [NC]"),
        stringsAsFactors=FALSE
      ),
      caption="auclast, cmax, aucinf.obs: geometric mean and geometric coefficient of variation; tmax: median and range; half.life: arithmetic mean and standard deviation"
    ),
    info="summary of PKNCAresults with some missing values results in NA for spread"
  )

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
  suppressWarnings(myresult <- pk.nca(mydata))
  mysummary <- summary(myresult)
  expect_equal(
    mysummary,
    as_summary_PKNCAresults(
      data.frame(
        start=0,
        end=c(24, Inf),
        treatment="Trt 1",
        N="2",
        auclast=c("NC", "."),
        cmax=c(".", "NC"),
        tmax=c(".", "NC"),
        half.life=c(".", "NC"),
        aucinf.obs=c(".", "NC"),
        stringsAsFactors=FALSE),
      caption="auclast, cmax, aucinf.obs: geometric mean and geometric coefficient of variation; tmax: median and range; half.life: arithmetic mean and standard deviation"
    ),
    info="summary of PKNCAresults without most results gives NC"
  )

  mysummary <- summary(myresult,
                       not_requested="NR",
                       not_calculated="NoCalc")
  expect_equal(
    mysummary,
    as_summary_PKNCAresults(
      data.frame(
        start=0,
        end=c(24, Inf),
        treatment="Trt 1",
        N="2",
        auclast=c("NoCalc", "NR"),
        cmax=c("NR", "NoCalc"),
        tmax=c("NR", "NoCalc"),
        half.life=c("NR", "NoCalc"),
        aucinf.obs=c("NR", "NoCalc"),
        stringsAsFactors=FALSE
      ),
      caption="auclast, cmax, aucinf.obs: geometric mean and geometric coefficient of variation; tmax: median and range; half.life: arithmetic mean and standard deviation"
    ),
    info="Summary respects the not.requested.string and not.calculated.string"
  )

  mysummary <- summary(myresult,
                       summarize_n=FALSE,
                       not_requested="NR",
                       not_calculated="NoCalc")
  expect_equal(
    mysummary,
    as_summary_PKNCAresults(
      data.frame(
        start=0,
        end=c(24, Inf),
        treatment="Trt 1",
        auclast=c("NoCalc", "NR"),
        cmax=c("NR", "NoCalc"),
        tmax=c("NR", "NoCalc"),
        half.life=c("NR", "NoCalc"),
        aucinf.obs=c("NR", "NoCalc"),
        stringsAsFactors=FALSE),
      caption="auclast, cmax, aucinf.obs: geometric mean and geometric coefficient of variation; tmax: median and range; half.life: arithmetic mean and standard deviation"
    ),
    info="N is optionally omitted"
  )
})

test_that("PKNCAresults summary counts N correctly", {
  tmpconc <- generate.conc(2, 1, 0:24)
  tmpdose <- generate.dose(tmpconc)
  myconc <- PKNCAconc(tmpconc, formula=conc~time|treatment+ID)
  mydose <- PKNCAdose(tmpdose, formula=dose~time|treatment+ID)
  mydata <- PKNCAdata(myconc, mydose, intervals = data.frame(start = 0, end = c(24, Inf), cmax = TRUE))
  myresult <- pk.nca(mydata)

  # Testing the summarization
  mysummary_two_row <- summary(myresult)
  expect_warning(expect_warning(
    mysummary_one_row <- summary(myresult, drop_group = c("ID", "end")),
    "Some subjects may have more than one result for cmax"),
    "drop.group including start or end may result in incorrect groupings"
  )
  expect_equal(mysummary_two_row$N, c("2", "2"))
  expect_equal(mysummary_one_row$N, "2")

  # No subject identifier
  tmpconc <- generate.conc(1, 1, 0:24)
  tmpdose <- generate.dose(tmpconc)
  myconc <- PKNCAconc(tmpconc, formula=conc~time)
  mydose <- PKNCAdose(tmpdose, formula=dose~time)
  mydata <- PKNCAdata(myconc, mydose, intervals = data.frame(start = 0, end = c(24, Inf), cmax = TRUE))
  myresult <- pk.nca(mydata)

  mysummary_one_subject <- summary(myresult)
  expect_false("N" %in% names(mysummary_one_subject))
  expect_warning(
    mysummary_one_subject_askn <- summary(myresult, summarize_n = TRUE),
    "summarize_n was requested, but no subject column exists"
  )
  expect_false("N" %in% names(mysummary_one_subject_askn))
})

test_that("dropping `start` and `end` from groups is allowed with a warning.", {
  tmpconc <- generate.conc(2, 1, 0:24)
  tmpdose <- generate.dose(tmpconc)
  myconc <- PKNCAconc(tmpconc, formula=conc~time|treatment+ID)
  mydose <- PKNCAdose(tmpdose, formula=dose~time|treatment+ID)
  mydata <- PKNCAdata(myconc, mydose)
  myresult <- pk.nca(mydata)

  expect_warning(
    current_summary <- summary(myresult, drop_group=c("ID", "start")),
    regexp="drop.group including start or end may result", fixed=TRUE
  )
  expect_false("start" %in% names(current_summary))
})

test_that("summary.PKNCAresults manages exclusions as missing not as non-existent.", {
  # Note that generate.conc sets the random seed, so it doesn't have
  # to happen here.
  tmpconc <- generate.conc(2, 1, 0:24)
  tmpdose <- generate.dose(tmpconc)
  myconc <- PKNCAconc(tmpconc, formula=conc~time|treatment+ID)
  mydose <- PKNCAdose(tmpdose, formula=dose~time|treatment+ID)
  mydata <- PKNCAdata(myconc, mydose)
  myresult <- pk.nca(mydata)
  myresult_excluded <-
    exclude(
      myresult,
      reason="testing",
      mask=with(as.data.frame(myresult),
                PPTESTCD %in% "auclast" & ID %in% 1)
    )
  myresult_excluded2 <-
    exclude(
      myresult,
      reason="testing",
      mask=with(as.data.frame(myresult),
                PPTESTCD %in% "auclast")
    )
  # Testing the summarization
  mysummary <- summary(myresult)
  mysummary_excluded <- summary(myresult_excluded)
  mysummary_excluded2 <- summary(myresult_excluded2)
  expect_equal(
    mysummary,
    as_summary_PKNCAresults(
      data.frame(
        start=0,
        end=c(24, Inf),
        treatment="Trt 1",
        N="2",
        auclast=c("13.8 [2.51]", "."),
        cmax=c(".", "0.970 [4.29]"),
        tmax=c(".", "3.00 [2.00, 4.00]"),
        half.life=c(".", "14.2 [2.79]"),
        aucinf.obs=c(".", "20.5 [6.84]"),
        stringsAsFactors=FALSE
      ),
      caption="auclast, cmax, aucinf.obs: geometric mean and geometric coefficient of variation; tmax: median and range; half.life: arithmetic mean and standard deviation"
    ),
    info="simple summary of PKNCAresults performs as expected"
  )
  expect_equal(
    mysummary_excluded,
    as_summary_PKNCAresults(
      data.frame(
        start=0,
        end=c(24, Inf),
        treatment="Trt 1",
        N="2",
        auclast=c("14.0 [NC]", "."),
        cmax=c(".", "0.970 [4.29]"),
        tmax=c(".", "3.00 [2.00, 4.00]"),
        half.life=c(".", "14.2 [2.79]"),
        aucinf.obs=c(".", "20.5 [6.84]"),
        stringsAsFactors=FALSE
      ),
      caption="auclast, cmax, aucinf.obs: geometric mean and geometric coefficient of variation; tmax: median and range; half.life: arithmetic mean and standard deviation"
    ),
    info="summary of PKNCAresults correctly excludes auclast when requested")
  expect_equal(
    mysummary_excluded2,
    as_summary_PKNCAresults(
      data.frame(
        start=0,
        end=c(24, Inf),
        treatment="Trt 1",
        N="2",
        auclast=c("NC", "."),
        cmax=c(".", "0.970 [4.29]"),
        tmax=c(".", "3.00 [2.00, 4.00]"),
        half.life=c(".", "14.2 [2.79]"),
        aucinf.obs=c(".", "20.5 [6.84]"),
        stringsAsFactors=FALSE
      ),
      caption="auclast, cmax, aucinf.obs: geometric mean and geometric coefficient of variation; tmax: median and range; half.life: arithmetic mean and standard deviation"
    ),
    info="summary of PKNCAresults correctly excludes all of auclast when requested"
  )
})

test_that("print.summary_PKNCAresults works", {
  tmpconc <- generate.conc(2, 1, 0:24)
  tmpdose <- generate.dose(tmpconc)
  myconc <- PKNCAconc(tmpconc, formula=conc~time|treatment+ID)
  mydose <- PKNCAdose(tmpdose, formula=dose~time|treatment+ID)
  mydata <- PKNCAdata(myconc, mydose)
  myresult <- pk.nca(mydata)

  expect_output(
    print(summary(myresult)),
    paste(
      " start end treatment N     auclast         cmax              tmax   half.life.*",
      "     0  24     Trt 1 2 13.8 \\[2.51\\]            .                 .           ..*",
      "     0 Inf     Trt 1 2           . 0.970 \\[4.29\\] 3.00 \\[2.00, 4.00\\] 14.2 \\[2.79\\].*",
      "",
      "Caption: auclast, cmax, aucinf.obs: geometric mean and geometric coefficient of variation; tmax: median and range; half.life: arithmetic mean and standard deviation",
      sep="\n"
    )
  )
})

test_that("summary pretty_name control", {
  tmpconc <- generate.conc(2, 1, 0:24)
  tmpdose <- generate.dose(tmpconc)
  myconc <- PKNCAconc(tmpconc, formula=conc~time|treatment+ID)
  mydose <- PKNCAdose(tmpdose, formula=dose~time|treatment+ID)
  mydata <- PKNCAdata(myconc, mydose)
  myresult <- pk.nca(mydata)

  d_units_orig <- pknca_units_table(concu="ng/mL", doseu="mg", amountu="mg", timeu="hr")
  d_units_std <-
    pknca_units_table(
      concu="ng/mL", doseu="mg", amountu="mg", timeu="hr",
      conversions=data.frame(PPORRESU="ng/mL", PPSTRESU="mg/mL")
    )
  mydata_orig <- PKNCAdata(myconc, mydose, units=d_units_orig)
  myresult_units_orig <- pk.nca(mydata_orig)

  s_plain <- summary(myresult)
  s_pretty <- summary(myresult, pretty_names=TRUE)
  s_plain_units <- summary(myresult_units_orig, pretty_names=FALSE)
  s_pretty_units <- summary(myresult_units_orig)
  expect_equal(
    names(s_plain),
    c("start", "end", "treatment", "N", "auclast", "cmax", "tmax", "half.life", "aucinf.obs")
  )
  expect_equal(
    names(s_pretty),
    c("Interval Start", "Interval End", "treatment", "N", "AUClast",
      "Cmax", "Tmax", "Half-life", "AUCinf,obs")
  )
  expect_equal(
    names(s_plain_units),
    c("start", "end", "treatment", "N", "auclast (hr*ng/mL)", "cmax (ng/mL)",
      "tmax (hr)", "half.life (hr)", "aucinf.obs (hr*ng/mL)")
  )
  expect_equal(
    names(s_pretty_units),
    c(
      "Interval Start", "Interval End", "treatment", "N", "AUClast (hr*ng/mL)",
      "Cmax (ng/mL)", "Tmax (hr)", "Half-life (hr)", "AUCinf,obs (hr*ng/mL)"
    )
  )
  # Captions use the pretty_names, if requested
  expect_equal(
    attr(s_plain, "caption"),
    "auclast, cmax, aucinf.obs: geometric mean and geometric coefficient of variation; tmax: median and range; half.life: arithmetic mean and standard deviation"
  )
  expect_equal(
    attr(s_pretty, "caption"),
    "AUClast, Cmax, AUCinf,obs: geometric mean and geometric coefficient of variation; Tmax: median and range; Half-life: arithmetic mean and standard deviation"
  )
  # Default for pretty_names are kept
  expect_equal(
    names(s_plain),
    names(summary(myresult, pretty_names=FALSE))
  )
  expect_equal(
    names(s_pretty_units),
    names(summary(myresult_units_orig, pretty_names=TRUE))
  )
})

test_that("roundingSummarize", {
  expect_error(
    roundingSummarize(1, "foo"),
    regexp="foo is not in the summarization instructions from PKNCA.set.summary"
  )

  PKNCA.set.summary(name="lambda.z.n.points", description="not a real parameter", point=mean, spread=sd, rounding=function(x) round(x, 1))
  expect_equal(roundingSummarize(1.2345, "lambda.z.n.points"), "1.2")
  PKNCA.set.summary(name="lambda.z.n.points", description="not a real parameter", point=mean, spread=sd, rounding=list(round=1))
  expect_equal(roundingSummarize(1.2345, "lambda.z.n.points"), "1.2")

  # reset it
  PKNCA.set.summary(
    name="lambda.z.n.points",
    description="median and range",
    point=business.median,
    spread=business.range
  )
})
