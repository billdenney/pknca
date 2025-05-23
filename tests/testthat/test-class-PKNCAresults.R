test_that("PKNCAresults object creation", {
  minimal_result <- PKNCAresults(data.frame(a=1), data=list())
  expect_equal(minimal_result$columns$exclude, "exclude")
  result_with_exclude_col <- PKNCAresults(data.frame(exclude=1), data=list())
  expect_equal(result_with_exclude_col$columns$exclude, "exclude.exclude")
})

test_that("PKNCAresults generation", {
  # Note that generate.conc sets the random seed, so it doesn't have
  # to happen here.
  d_conc <- generate.conc(2, 1, 0:24)
  d_dose <- generate.dose(d_conc)
  o_conc <- PKNCAconc(d_conc, formula=conc~time|treatment+ID)
  o_dose <- PKNCAdose(d_dose, formula=dose~time|treatment+ID)
  o_data <- PKNCAdata(o_conc, o_dose)
  o_result <- pk.nca(o_data)

  expect_equal(
    names(o_result),
    c("result", "data", "columns"),
    info="Make sure that the result has the expected names (and only the expected names) in it."
  )
  expect_true(
    checkProvenance(o_result),
    info="Provenance exists and can be confirmed on results"
  )

  # Test each of the pieces for o_result for accuracy
  expect_equal(
    o_result$data, {
      tmp <- o_data
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
    o_result$result,
    verify.result,
    tolerance=0.001,
    info="The specific order of the levels isn't important-- the fact that they are factors and that the set doesn't change is important."
  )

  # Test conversion to a data.frame
  expect_equal(
    as.data.frame(o_result),
    verify.result,
    tolerance=0.001,
    info="Conversion of PKNCAresults to a data.frame in long format (default long format)"
  )
  expect_equal(
    as.data.frame(o_result, out_format="long"),
    verify.result,
    tolerance=0.001,
    info="Conversion of PKNCAresults to a data.frame in long format (specifying long format)"
  )
  expect_equal(
    as.data.frame(o_result, out_format="wide"),
    tidyr::spread(verify.result, key="PPTESTCD", value="PPORRES"),
    tolerance=0.001,
    info="Conversion of PKNCAresults to a data.frame in wide format (specifying wide format)"
  )

  d_conc <- generate.conc(2, 1, 0:24)
  d_dose <- generate.dose(d_conc)
  o_conc <- PKNCAconc(d_conc, formula=conc~time|treatment+ID)
  o_dose <- PKNCAdose(d_dose, formula=dose~time|treatment+ID)
  o_data <- PKNCAdata(o_conc, o_dose, intervals=data.frame(start=0, end=12, aucint.inf.obs=TRUE))
  o_result <- pk.nca(o_data)

  d_conc12 <- d_conc
  d_conc12$time <- d_conc$time + 12
  d_dose12 <- generate.dose(d_conc12)
  o_conc12 <- PKNCAconc(d_conc12, formula=conc~time|treatment+ID)
  o_dose12 <- PKNCAdose(d_dose12, formula=dose~time|treatment+ID)
  o_data12 <- PKNCAdata(o_conc12, o_dose12, intervals=data.frame(start=12, end=24, aucint.inf.obs=TRUE))
  o_result12 <- pk.nca(o_data12)
  comparison_orig <- as.data.frame(o_result)
  comparison_12 <- as.data.frame(o_result12)
  expect_equal(
    comparison_orig$PPORRES[comparison_orig$PPTESTCD %in% "aucint.inf.obs"],
    comparison_12$PPORRES[comparison_12$PPTESTCD %in% "aucint.inf.obs"],
    info="Time shift does not affect aucint calculations."
  )
})

test_that("PKNCAresults has exclude, when applicable", {
  d_conc <- generate.conc(2, 1, 0:24)
  d_conc$conc[d_conc$ID %in% 2] <- 0
  d_dose <- generate.dose(d_conc)
  o_conc <- PKNCAconc(d_conc, conc~time|treatment+ID)
  o_dose <- PKNCAdose(d_dose, dose~time|treatment+ID)
  o_data <- PKNCAdata(o_conc, o_dose)
  # Not capturing the warning due to R bug
  # https://bugs.r-project.org/bugzilla3/show_bug.cgi?id=17122
  #expect_warning(o_result <- pk.nca(o_data),
  #               regexp="Too few points for half-life calculation")
  suppressWarnings(o_result <- pk.nca(o_data))
  o_result_df <- as.data.frame(o_result)
  expect_true(
    all(o_result_df$PPTESTCD %in%
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
      o_result_df$exclude[
        o_result_df$ID == 2 &
          o_result_df$PPTESTCD %in%
          c("lambda.z", "r.squared", "adj.r.squared", "lambda.z.time.first",
            "lambda.z.n.points", "clast.pred", "half.life", "span.ratio")
        ]
    ),
    "Too few points for half-life calculation (min.hl.points=3 with only 0 points)",
    info="exclusions are propogated to results"
  )
  expect_equal(
    unique(
      o_result_df$exclude[
        !(o_result_df$ID == 2 &
            o_result_df$PPTESTCD %in%
            c("lambda.z", "r.squared", "adj.r.squared", "lambda.z.time.first",
              "lambda.z.n.points", "clast.pred", "half.life", "span.ratio")
        )
        ]
    ),
    NA_character_,
    info="exclusions are propogated to results only when applicable"
  )
})

test_that("ptr works as a parameter", {
  d_conc <- generate.conc(2, 1, 0:24)
  d_dose <- generate.dose(d_conc)
  o_conc <- PKNCAconc(d_conc, formula=conc~time|treatment+ID)
  o_dose <- PKNCAdose(d_dose, formula=dose~time|treatment+ID)
  myinterval <- data.frame(start=0, end=24, ptr=TRUE)
  o_data <- PKNCAdata(o_conc, o_dose, intervals=myinterval)
  o_result <- pk.nca(o_data)
  ptr_result <- as.data.frame(o_result)
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
  d_conc <- generate.conc(2, 1, 0:24)
  d_dose <- generate.dose(d_conc)
  o_conc <- PKNCAconc(d_conc, formula=conc~time|treatment+ID)
  o_dose <- PKNCAdose(d_dose, formula=dose~time|treatment+ID)
  o_data <- PKNCAdata(o_conc, o_dose)
  o_result <- pk.nca(o_data)

  d_units_orig <- pknca_units_table(concu="ng/mL", doseu="mg", amountu="mg", timeu="hr")
  d_units_std <-
    pknca_units_table(
      concu="ng/mL", doseu="mg", amountu="mg", timeu="hr",
      conversions=data.frame(PPORRESU="ng/mL", PPSTRESU="mg/mL")
    )
  o_data_orig <- PKNCAdata(o_conc, o_dose, units=d_units_orig)
  o_result_units_orig <- pk.nca(o_data_orig)
  o_data_std <- PKNCAdata(o_conc, o_dose, units=d_units_std)
  o_result_units_std <- pk.nca(o_data_std)

  # Summaries are the same except for the column names
  expect_equal(
    unname(summary(o_result)),
    unname(summary(o_result_units_orig)),
    # The caption attribute will differ
    ignore_attr = TRUE
  )
  expect_equal(
    summary(o_result_units_orig) %>% dplyr::select(-`Cmax (ng/mL)`),
    summary(o_result_units_std) %>% dplyr::select(-`Cmax (mg/mL)`)
  )
  # The units are converted to standard units, if requested
  expect_equal(
    summary(o_result_units_orig)$`Cmax (ng/mL)`,
    c(".", "0.970 [4.29]")
  )
  expect_equal(
    summary(o_result_units_std)$`Cmax (mg/mL)`,
    c(".", "9.70e-7 [4.29]")
  )
  # Wide conversion works for original and standardized units
  df_wide_orig <- as.data.frame(o_result_units_orig, out_format="wide")
  df_wide_std <- as.data.frame(o_result_units_std, out_format="wide")
  expect_equal(
    as.data.frame(o_result, out_format="wide"),
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
  d_conc1 <- generate.conc(2, 1, 0:24)
  d_conc1$analyte <- "drug1"
  d_conc2 <- d_conc1
  d_conc2$analyte <- "drug2"
  d_conc <- rbind(d_conc1, d_conc2)

  d_dose <- generate.dose(d_conc)
  o_conc <- PKNCAconc(d_conc, formula=conc~time|treatment+ID/analyte)
  o_dose <- PKNCAdose(d_dose, formula=dose~time|treatment+ID)
  o_data <- PKNCAdata(o_conc, o_dose)
  o_result <- pk.nca(o_data)

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
  o_data_std <- PKNCAdata(o_conc, o_dose, units=d_units_std)
  o_result_units_std <- pk.nca(o_data_std)
  summary_o_result_units_std <- summary(o_result_units_std)
  # Everything is the same between analytes except for "cmax"
  for (nm in setdiff(names(summary_o_result_units_std), c("analyte", "Cmax"))) {
    expect_equal(
      summary_o_result_units_std[[nm]][1:2],
      summary_o_result_units_std[[nm]][3:4]
    )
  }
  # Different units in the same column are shown in the cell
  expect_equal(
    summary_o_result_units_std$Cmax,
    c(".", "9.70e-7 [4.29] mg/mL", ".", "1.94 [4.29] mmol/L")
  )

  # I can't think of a way to trigger this error without explicit manipulation.
  o_result_units_manipulated <- o_result_units_std
  o_result_units_manipulated$result$PPSTRESU[o_result_units_manipulated$result$PPTESTCD %in% "auclast"][1] <- "foo"
  expect_error(
    summary(o_result_units_manipulated),
    regexp="Multiple units cannot be summarized together.  For auclast, trying to combine: foo, hr*ng/mL",
    fixed=TRUE
  )
})

test_that("getGroups.PKNCAresults", {
  d_conc <- generate.conc(2, 1, 0:24)
  d_dose <- generate.dose(d_conc)
  o_conc <- PKNCAconc(d_conc, formula=conc~time|treatment+ID)
  o_dose <- PKNCAdose(d_dose, formula=dose~time|treatment+ID)
  o_data <- PKNCAdata(o_conc, o_dose)
  o_result <- pk.nca(o_data)

  expect_equal(
    getGroups(o_result, level="treatment"),
    o_result$result[, "treatment", drop=FALSE]
  )
  expect_equal(
    getGroups(o_result, level=factor("treatment")),
    o_result$result[, "treatment", drop=FALSE]
  )
  expect_error(
    getGroups(o_result, level="foo"),
    regexp="Not all levels are listed in the group names.  Missing levels are: foo"
  )
  expect_equal(
    getGroups(o_result, level=2),
    o_result$result[, c("treatment", "ID")]
  )
  expect_equal(
    getGroups(o_result, level=2:3),
    o_result$result[, c("ID", "start")]
  )
})

test_that("group_vars.PKNCAresult", {
  o_conc_group <- PKNCAconc(as.data.frame(datasets::Theoph), conc~Time|Subject)
  o_data_group <- PKNCAdata(o_conc_group, intervals = data.frame(start = 0, end = 1, cmax = TRUE))
  suppressMessages(o_nca_group <- pk.nca(o_data_group))

  expect_equal(dplyr::group_vars(o_nca_group), "Subject")

  # Check that it works without groupings as expected [empty]
  o_conc_nongroup <- PKNCAconc(as.data.frame(datasets::Theoph)[datasets::Theoph$Subject == 1,], conc~Time)
  o_data_nogroup <- PKNCAdata(o_conc_nongroup, intervals = data.frame(start = 0, end = 1, cmax = TRUE))
  suppressMessages(o_nca_nogroup <- pk.nca(o_data_nogroup))

  expect_equal(dplyr::group_vars(o_nca_nogroup), character(0))
})

test_that("as.data.frame.PKNCAresults can filter for only requested parameters", {
  d_conc <- generate.conc(2, 1, 0:24)
  d_dose <- generate.dose(d_conc)
  o_conc <- PKNCAconc(d_conc, formula=conc~time|treatment+ID)
  o_dose <- PKNCAdose(d_dose, formula=dose~time|treatment+ID)
  o_data <- PKNCAdata(o_conc, o_dose, intervals = data.frame(start = 0, end = Inf, half.life = TRUE))
  o_result <- pk.nca(o_data)

  expect_equal(nrow(as.data.frame(o_result)), 20)
  expect_equal(nrow(as.data.frame(o_result, filter_requested = TRUE)), 2)
})

test_that("as.data.frame.PKNCAresults can filter to remove excluded parameters", {
  d_conc <- generate.conc(2, 1, c(0, 2, 6, 12, 24))
  d_dose <- generate.dose(d_conc)
  o_conc <- PKNCAconc(d_conc, formula=conc~time|treatment+ID)
  o_dose <- PKNCAdose(d_dose, formula=dose~time|treatment+ID)
  o_data <- PKNCAdata(o_conc, o_dose, intervals = data.frame(start = 0, end = Inf, half.life = TRUE))
  o_result <- exclude(pk.nca(o_data), FUN = exclude_nca_span.ratio(1))

  expect_equal(nrow(as.data.frame(o_result)), 20)
  expect_equal(nrow(as.data.frame(o_result, filter_excluded = TRUE)), 12)
})
