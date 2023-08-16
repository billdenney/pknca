source("generate.data.R")

test_that("PKNCAresults object creation", {
  minimal_result <- PKNCAresults(data.frame(a=1), data=list())
  expect_equal(minimal_result$columns$exclude, "exclude")
  result_with_exclude_col <- PKNCAresults(data.frame(exclude=1), data=list())
  expect_equal(result_with_exclude_col$columns$exclude, "exclude.exclude")
})

test_that("PKNCAresults generation", {
  # Note that generate.conc sets the random seed, so it doesn't have
  # to happen here.
  tmpconc <- generate.conc(2, 1, 0:24)
  tmpdose <- generate.dose(tmpconc)
  myconc <- PKNCAconc(tmpconc, formula=conc~time|treatment+ID)
  mydose <- PKNCAdose(tmpdose, formula=dose~time|treatment+ID)
  mydata <- PKNCAdata(myconc, mydose)
  myresult <- pk.nca(mydata)

  expect_equal(
    names(myresult),
    c("result", "data", "columns"),
    info="Make sure that the result has the expected names (and only the expected names) in it."
  )
  expect_true(
    checkProvenance(myresult),
    info="Provenance exists and can be confirmed on results"
  )

  # Test each of the pieces for myresult for accuracy
  expect_equal(
    myresult$data, {
      tmp <- mydata
      # The options should be the default options after the
      # calculations are done.
      tmp$options <- PKNCA.options()
      tmp
    }, info="The data is just a copy of the input data plus an instantiation of the PKNCA.options"
  )

  verify.result <-
    tibble::tibble(
      treatment="Trt 1",
      ID=as.integer(rep(c(1, 2), each=14)),
      start=0,
      end=c(24, rep(Inf, 13),
            24, rep(Inf, 13)),
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
      exclude=NA_character_
    )
  expect_equal(
    myresult$result,
    verify.result,
    tolerance=0.001,
    info="The specific order of the levels isn't important-- the fact that they are factors and that the set doesn't change is important."
  )

  # Test conversion to a data.frame
  expect_equal(
    as.data.frame(myresult),
    verify.result,
    tolerance=0.001,
    info="Conversion of PKNCAresults to a data.frame in long format (default long format)"
  )
  expect_equal(
    as.data.frame(myresult, out.format="long"),
    verify.result,
    tolerance=0.001,
    info="Conversion of PKNCAresults to a data.frame in long format (specifying long format)"
  )
  expect_equal(
    as.data.frame(myresult, out.format="wide"),
    tidyr::spread(verify.result, key="PPTESTCD", value="PPORRES"),
    tolerance=0.001,
    info="Conversion of PKNCAresults to a data.frame in wide format (specifying wide format)"
  )

  tmpconc <- generate.conc(2, 1, 0:24)
  tmpdose <- generate.dose(tmpconc)
  myconc <- PKNCAconc(tmpconc, formula=conc~time|treatment+ID)
  mydose <- PKNCAdose(tmpdose, formula=dose~time|treatment+ID)
  mydata <- PKNCAdata(myconc, mydose, intervals=data.frame(start=0, end=12, aucint.inf.obs=TRUE))
  myresult <- pk.nca(mydata)

  tmpconc12 <- tmpconc
  tmpconc12$time <- tmpconc$time + 12
  tmpdose12 <- generate.dose(tmpconc12)
  myconc12 <- PKNCAconc(tmpconc12, formula=conc~time|treatment+ID)
  mydose12 <- PKNCAdose(tmpdose12, formula=dose~time|treatment+ID)
  mydata12 <- PKNCAdata(myconc12, mydose12, intervals=data.frame(start=12, end=24, aucint.inf.obs=TRUE))
  myresult12 <- pk.nca(mydata12)
  comparison_orig <- as.data.frame(myresult)
  comparison_12 <- as.data.frame(myresult12)
  expect_equal(
    comparison_orig$PPORRES[comparison_orig$PPTESTCD %in% "aucint.inf.obs"],
    comparison_12$PPORRES[comparison_12$PPTESTCD %in% "aucint.inf.obs"],
    info="Time shift does not affect aucint calculations."
  )
})

test_that("PKNCAresults has exclude, when applicable", {
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
  myresult_df <- as.data.frame(myresult)
  expect_true(
    all(myresult_df$PPTESTCD %in%
          c(
            "adj.r.squared", "aucinf.obs", "auclast", "clast.obs",
            "clast.pred", "cmax", "half.life", "lambda.z", "lambda.z.n.points",
            "lambda.z.time.first", "r.squared", "span.ratio", "tlast", "tmax"
          )
    ),
    info="verify that only expected results are present"
  )
  expect_equal(
    unique(
      myresult_df$exclude[
        myresult_df$ID == 2 &
          myresult_df$PPTESTCD %in%
          c("lambda.z", "r.squared", "adj.r.squared", "lambda.z.time.first",
            "lambda.z.n.points", "clast.pred", "half.life", "span.ratio")
        ]
    ),
    "Too few points for half-life calculation (min.hl.points=3 with only 0 points)",
    info="exclusions are propogated to results"
  )
  expect_equal(
    unique(
      myresult_df$exclude[
        !(myresult_df$ID == 2 &
            myresult_df$PPTESTCD %in%
            c("lambda.z", "r.squared", "adj.r.squared", "lambda.z.time.first",
              "lambda.z.n.points", "clast.pred", "half.life", "span.ratio")
        )
        ]
    ),
    NA_character_,
    info="exclusions are propogated to results only when applicable"
  )
})

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
                       not.requested.string="NR",
                       not.calculated.string="NoCalc")
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
                       summarize.n.per.group=FALSE,
                       not.requested.string="NR",
                       not.calculated.string="NoCalc")
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
    mysummary_one_row <- summary(myresult, drop.group = c("ID", "end")),
    "Some subjects may have more than one result for cmax"),
    "drop.group including start or end may result in incorrect groupings"
  )
  expect_equal(mysummary_two_row$N, c("2", "2"))
  expect_equal(mysummary_one_row$N, "2")
})

test_that("dropping `start` and `end` from groups is allowed with a warning.", {
  tmpconc <- generate.conc(2, 1, 0:24)
  tmpdose <- generate.dose(tmpconc)
  myconc <- PKNCAconc(tmpconc, formula=conc~time|treatment+ID)
  mydose <- PKNCAdose(tmpdose, formula=dose~time|treatment+ID)
  mydata <- PKNCAdata(myconc, mydose)
  myresult <- pk.nca(mydata)

  expect_warning(
    current_summary <- summary(myresult, drop.group=c("ID", "start")),
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

test_that("ptr works as a parameter", {
  tmpconc <- generate.conc(2, 1, 0:24)
  tmpdose <- generate.dose(tmpconc)
  myconc <- PKNCAconc(tmpconc, formula=conc~time|treatment+ID)
  mydose <- PKNCAdose(tmpdose, formula=dose~time|treatment+ID)
  myinterval <- data.frame(start=0, end=24, ptr=TRUE)
  mydata <- PKNCAdata(myconc, mydose, intervals=myinterval)
  myresult <- pk.nca(mydata)
  ptr_result <- as.data.frame(myresult)
  expect_equal(
    ptr_result$PPORRES[ptr_result$PPTESTCD %in% "ptr"],
    c(2.9055, 2.9885),
    tolerance=0.0001
  )
})

test_that("exclude values are maintained in derived parameters during automatic calculation (#112)", {
  my_conc <-
    data.frame(
      conc = c(0, 2.5, 3, 2.7, 2.3),
      time = 0:4,
      subject = 1
    )

  conc_obj <- PKNCAconc(my_conc, conc~time|subject)
  data_obj <-
    PKNCAdata(
      data.conc=conc_obj,
      intervals=
        data.frame(
          start=0,
          end=Inf,
          aucinf.obs=TRUE
        )
    )
  expect_message(
    expect_warning(
      results_obj <- pk.nca(data_obj),
      regexp="Too few points for half-life"
    ),
    regexp="No dose information provided"
  )
  d_results <- as.data.frame(results_obj)
  expect_equal(
    d_results$exclude[d_results$PPTESTCD == "aucinf.obs"],
    d_results$exclude[d_results$PPTESTCD == "half.life"]
  )
})

test_that("ctrough is correctly calculated", {
  my_conc <- data.frame(time=0:6, conc=2^(0:-6), subject=1)
  conc_obj <- PKNCAconc(my_conc, conc~time|subject)
  data_obj <-
    PKNCAdata(
      data.conc=conc_obj,
      intervals=
        data.frame(
          start=0,
          end=c(6, Inf),
          ctrough=TRUE
        )
    )
  expect_message(
    expect_equal(
      as.data.frame(pk.nca(data_obj))$PPORRES,
      c(2^-6, NA_real_)
    ),
    regexp="No dose information provided"
  )
})

test_that("single subject, ungrouped data works (#74)", {
  my_conc <- data.frame(time=0:6, conc=2^(0:-6))
  conc_obj <- PKNCAconc(my_conc, conc~time)
  data_obj <-
    PKNCAdata(
      data.conc=conc_obj,
      intervals=
        data.frame(
          start=0,
          end=Inf,
          cmax=TRUE
        )
    )
  expect_message(
    expect_equal(
      as.data.frame(pk.nca(data_obj))$PPORRES,
      1
    ),
    regexp="No dose information provided",
  )
})

test_that("units work for calculations and summaries with one set of units across all analytes", {
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
  mydata_std <- PKNCAdata(myconc, mydose, units=d_units_std)
  myresult_units_std <- pk.nca(mydata_std)

  # Summaries are the same except for the column names
  expect_equal(
    unname(summary(myresult)),
    unname(summary(myresult_units_orig)),
    # The caption attribute will differ
    ignore_attr = TRUE
  )
  expect_equal(
    summary(myresult_units_orig) %>% dplyr::select(-`Cmax (ng/mL)`),
    summary(myresult_units_std) %>% dplyr::select(-`Cmax (mg/mL)`)
  )
  # The units are converted to standard units, if requested
  expect_equal(
    summary(myresult_units_orig)$`Cmax (ng/mL)`,
    c(".", "0.970 [4.29]")
  )
  expect_equal(
    summary(myresult_units_std)$`Cmax (mg/mL)`,
    c(".", "9.70e-7 [4.29]")
  )
  # Wide conversion works for original and standardized units
  df_wide_orig <- as.data.frame(myresult_units_orig, out.format="wide")
  df_wide_std <- as.data.frame(myresult_units_std, out.format="wide")
  expect_equal(
    as.data.frame(myresult, out.format="wide"),
    # The difference is the addition of units to the column names
    df_wide_orig %>%
      dplyr::rename_with(.fn=gsub, pattern=" \\(.*$", replacement="")
  )
  expect_true(
    all(
      names(df_wide_orig) %in% c("treatment", "ID", "start", "end", "exclude") |
        grepl(x=names(df_wide_orig), pattern=" (", fixed=TRUE)
    )
  )
  # Everything is the same unless it is a concentration which has been converted
  expect_equal(
    df_wide_orig %>% dplyr::select(-`cmax (ng/mL)`, -`clast.obs (ng/mL)`, -`clast.pred (ng/mL)`),
    df_wide_std %>% dplyr::select(-`cmax (mg/mL)`, -`clast.obs (mg/mL)`, -`clast.pred (mg/mL)`)
  )
  # Concentration conversion works correctly
  expect_equal(
    df_wide_orig$`cmax (ng/mL)`,
    df_wide_std$`cmax (mg/mL)`*1e6
  )
})

test_that("units work for calculations and summaries with one set of units across all analytes", {
  tmpconc1 <- generate.conc(2, 1, 0:24)
  tmpconc1$analyte <- "drug1"
  tmpconc2 <- tmpconc1
  tmpconc2$analyte <- "drug2"
  tmpconc <- rbind(tmpconc1, tmpconc2)

  tmpdose <- generate.dose(tmpconc)
  myconc <- PKNCAconc(tmpconc, formula=conc~time|treatment+ID/analyte)
  mydose <- PKNCAdose(tmpdose, formula=dose~time|treatment+ID)
  mydata <- PKNCAdata(myconc, mydose)
  myresult <- pk.nca(mydata)

  d_units_std1 <-
    pknca_units_table(
      concu="ng/mL", doseu="mg", amountu="mg", timeu="hr",
      conversions=data.frame(PPORRESU="ng/mL", PPSTRESU="mg/mL")
    )
  d_units_std1$analyte <- "drug1"
  d_units_std2 <-
    pknca_units_table(
      concu="ng/mL", doseu="mg", amountu="mg", timeu="hr",
      conversions=data.frame(PPORRESU="ng/mL", PPSTRESU="mmol/L", conversion_factor=2)
    )
  d_units_std2$analyte <- "drug2"
  d_units_std <- rbind(d_units_std1, d_units_std2)
  mydata_std <- PKNCAdata(myconc, mydose, units=d_units_std)
  myresult_units_std <- pk.nca(mydata_std)
  summary_myresult_units_std <- summary(myresult_units_std)
  # Everything is the same between analytes except for "cmax"
  for (nm in setdiff(names(summary_myresult_units_std), c("analyte", "Cmax"))) {
    expect_equal(
      summary_myresult_units_std[[nm]][1:2],
      summary_myresult_units_std[[nm]][3:4]
    )
  }
  # Different units in the same column are shown in the cell
  expect_equal(
    summary_myresult_units_std$Cmax,
    c(".", "9.70e-7 [4.29] mg/mL", ".", "1.94 [4.29] mmol/L")
  )

  # I can't think of a way to trigger this error without explicit manipulation.
  myresult_units_manipulated <- myresult_units_std
  myresult_units_manipulated$result$PPSTRESU[myresult_units_manipulated$result$PPTESTCD %in% "auclast"][1] <- "foo"
  expect_error(
    summary(myresult_units_manipulated),
    regexp="Multiple units cannot be summarized together.  For auclast, trying to combine: foo, hr*ng/mL",
    fixed=TRUE
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

test_that("getGroups.PKNCAresults", {
  tmpconc <- generate.conc(2, 1, 0:24)
  tmpdose <- generate.dose(tmpconc)
  myconc <- PKNCAconc(tmpconc, formula=conc~time|treatment+ID)
  mydose <- PKNCAdose(tmpdose, formula=dose~time|treatment+ID)
  mydata <- PKNCAdata(myconc, mydose)
  myresult <- pk.nca(mydata)

  expect_equal(
    getGroups(myresult, level="treatment"),
    myresult$result[, "treatment", drop=FALSE]
  )
  expect_equal(
    getGroups(myresult, level=factor("treatment")),
    myresult$result[, "treatment", drop=FALSE]
  )
  expect_error(
    getGroups(myresult, level="foo"),
    regexp="Not all levels are listed in the group names.  Missing levels are: foo"
  )
  expect_equal(
    getGroups(myresult, level=2),
    myresult$result[, c("treatment", "ID")]
  )
  expect_equal(
    getGroups(myresult, level=2:3),
    myresult$result[, c("ID", "start")]
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
