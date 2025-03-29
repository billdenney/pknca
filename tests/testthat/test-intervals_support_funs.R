# Create sample data for testing
d_conc <- data.frame(
  conc = c(1, 0.6, 0.2, 0.1, 0.9, 0.4, 1.2, 0.8, 0.3, 0.2, 1.1, 0.5),
  time = rep(0:5, 2),
  analyte = rep(c("Analyte1", "Analyte2"), each = 6),
  include_hl = c(FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, NA, TRUE, TRUE, TRUE, TRUE),
  ID = rep(1:2, each = 6)
)

d_dose <- data.frame(
  dose = c(100, 200),
  time = c(0, 0),
  ID = c(1, 2)
)

intervals <- data.frame(
  start = c(0, 0, 0),
  end = c(24, 48, Inf),
  half.life = c(TRUE, TRUE, TRUE),
  cmax = c(TRUE, TRUE, TRUE),
  impute = c("start_conc0,start_predose", "start_predose", "start_conc0"),
  analyte = c("Analyte1", "Analyte2", "Analyte1"),
  ID = c(1, 2, 1)
)

o_conc <- PKNCA::PKNCAconc(d_conc, conc ~ time | ID / analyte, include_half.life = "include_hl")
o_dose <- PKNCA::PKNCAdose(d_dose, dose ~ time | ID)
o_data <- PKNCA::PKNCAdata(o_conc, o_dose, intervals = intervals)

describe("interval_add_impute", {
  it("adds the impute method in the impute column of a dummy intervals dataframe", {
    simple_df <- data.frame(
      cmax = TRUE,
      impute =  c("", "m0", "m0,m1")
    )
    expected_res <- data.frame(
      cmax = TRUE,
      impute = c("mlast", "m0,mlast", "m0,m1,mlast")
    )
    res <- interval_add_impute(simple_df, target_impute = "mlast")
    expect_equal(res, expected_res)
  })

  it("throws an error if either data or target_impute is missing", {
    expect_error(interval_add_impute(o_data), "Both 'data' and 'target_impute' must be provided.")
  })

  it("throws an error for non-character target_impute", {
    expect_error(interval_add_impute(o_data, target_impute = 123),
                 "'target_impute' must be a character string.")
  })

  it("throws an error when input data is not a proper format object", {
    expect_error(interval_add_impute(data = o_conc, target_impute = "start_conc0"))
    expect_no_error(interval_add_impute(data = o_data, target_impute = "start_conc0"))
  })

  it("throws an error for unknown target_params", {
    expect_error(interval_add_impute(o_data,
                                     target_impute = "start_conc0",
                                     target_params = "unknown_param"))
  })

  it("handles impute column with FALSE values correctly", {
    o_data_with_na_impute <- o_data
    o_data_with_na_impute$intervals$impute <- NA_character_
    expected_result <- data.frame(analyte = c("Analyte1", "Analyte2", "Analyte1"),
                                  half.life = c(TRUE, TRUE, TRUE),
                                  cmax = c(TRUE, TRUE, TRUE),
                                  impute = rep("new_impute", 3))
    result <- interval_add_impute(o_data_with_na_impute, target_impute = "new_impute")
    expect_equal(result$intervals[, c("analyte", "half.life", "cmax", "impute")], expected_result)
  })

  it("reports an error when the impute column is not a character", {
    o_data_not_character_impute <- o_data
    o_data_not_character_impute$intervals$impute <- 1
    expect_error(interval_add_impute(o_data_not_character_impute, target_impute = "new_impute"),
                 "The 'impute' column in the intervals data.frame must be a character column.")
  })

  it("warns and makes no changes when target_impute is NA or empty", {
    expect_warning({
      result <- interval_add_impute(o_data, target_impute = NA_character_)
      expect_equal(result, o_data)
    },
    "No impute method specified. No changes made."
    )

    expect_warning({
      result <- interval_add_impute(o_data, target_impute = "")
      expect_equal(result, o_data)
    },
    "No impute method specified. No changes made."
    )
  })

  it("creates missing impute col as NA_char & adds impute", {
    d_no_imp <- o_data
    d_no_imp$intervals$impute <- NULL
    res <- interval_add_impute(d_no_imp, target_impute = "new_impute")
    expect_equal(res$intervals, transform(d_no_imp$intervals, impute = "new_impute"))
    res <- interval_add_impute(d_no_imp$intervals, target_impute = "new_impute")
    expect_equal(res, transform(d_no_imp$intervals, impute = "new_impute"))
  })

  it("with no optional parameters uses all, with new intervals below", {
    expected_result <- data.frame(analyte = c("Analyte1", "Analyte2", "Analyte1"),
                                  half.life = c(TRUE, TRUE, TRUE),
                                  cmax = c(TRUE, TRUE, TRUE),
                                  impute = c("start_conc0,start_predose,new_impute",
                                             "start_predose,new_impute",
                                             "start_conc0,new_impute"))
    result <- interval_add_impute(o_data, target_impute = "new_impute")
    expect_equal(result$intervals[, c("analyte", "half.life", "cmax", "impute")], expected_result)
  })

  it("handles specified target_params correctly", {
    expected_result_half_life <- data.frame(
      analyte = c("Analyte1", "Analyte2", "Analyte1"),
      half.life = c(TRUE, TRUE, TRUE),
      impute = c("start_conc0,start_predose,new_impute",
                 "start_predose,new_impute",
                 "start_conc0,new_impute")
    )
    expected_result_cmax <- o_data$intervals[o_data$intervals$cmax,
                                             c("analyte", "cmax", "impute")] |>
      `rownames<-`(NULL)
    result <- interval_add_impute(o_data, target_impute = "new_impute", target_params = "half.life")
    expect_equal(result$intervals[result$intervals$half.life & !is.na(result$intervals$half.life),
                                  c("analyte", "half.life", "impute")] |>
                   `rownames<-`(NULL), expected_result_half_life)
    expect_equal(result$intervals[result$intervals$cmax & !is.na(result$intervals$cmax),
                                  c("analyte", "cmax", "impute")] |>
                   `rownames<-`(NULL), expected_result_cmax)
  })

  it("handles target_groups correctly", {
    expected_result_analyte1 <- data.frame(analyte = c("Analyte1", "Analyte1"),
                                           half.life = c(TRUE, TRUE),
                                           cmax = c(TRUE, TRUE),
                                           impute = c("start_conc0,start_predose,new_impute",
                                                      "start_conc0,new_impute"))
    expected_result_analyte2 <- o_data$intervals[o_data$intervals$analyte == "Analyte2",
                                                 c("analyte", "half.life", "cmax", "impute")] |>
      `rownames<-`(NULL)
    result <- interval_add_impute(o_data, target_impute = "new_impute",
                                  target_groups = data.frame(analyte = "Analyte1"))
    expect_equal(result$intervals[result$intervals$analyte == "Analyte1",
                                  c("analyte", "half.life", "cmax", "impute")] |>
                   `rownames<-`(NULL), expected_result_analyte1)
    expect_equal(result$intervals[result$intervals$analyte == "Analyte2",
                                  c("analyte", "half.life", "cmax", "impute")] |>
                   `rownames<-`(NULL), expected_result_analyte2)
  })

  it("handles multiple target_params correctly", {
    expected_result <- data.frame(analyte = c("Analyte1", "Analyte2", "Analyte1"),
                                  half.life = c(TRUE, TRUE, TRUE),
                                  cmax = c(TRUE, TRUE, TRUE),
                                  impute = c("start_conc0,start_predose,new_impute",
                                             "start_predose,new_impute",
                                             "start_conc0,new_impute"))
    result <- interval_add_impute(o_data,
                                  target_impute = "new_impute",
                                  target_params = c("half.life", "cmax"))
    expect_equal(result$intervals[, c("analyte", "half.life", "cmax", "impute")], expected_result)
  })

  it("makes no changes and warns when no matching intervals are found", {
    expect_warning({
      res <- interval_remove_impute(o_data,
                                    target_impute = "start_conc0",
                                    target_groups = data.frame(analyte = "Analyte3"))
      expect_equal(res, o_data)
    },
    paste0("No intervals found with the specified target parameters,",
           " groups and/or impute method. No changes made.")
    )
  })

  it("handles mixed TRUE/FALSE for cmax and half.life correctly", {
    intervals_mixed <- data.frame(
      start = c(0, 0, 0, 0),
      end = c(24, 48, Inf, 72),
      half.life = c(TRUE, FALSE, TRUE, FALSE),
      cmax = c(FALSE, TRUE, FALSE, TRUE),
      impute = c("start_conc0,start_predose", "start_predose", "start_conc0", "start_predose"),
      analyte = c("Analyte1", "Analyte2", "Analyte1", "Analyte2"),
      ID = c(1, 2, 1, 2)
    )

    o_data_mixed <- PKNCA::PKNCAdata(o_conc, o_dose, intervals = intervals_mixed)
    expected_result <- data.frame(analyte = c("Analyte1", "Analyte2", "Analyte1", "Analyte2"),
                                  half.life = c(TRUE, FALSE, TRUE, FALSE),
                                  cmax = c(FALSE, TRUE, FALSE, TRUE),
                                  impute = c("start_conc0,start_predose,new_impute",
                                             "start_predose,new_impute",
                                             "start_conc0,new_impute",
                                             "start_predose,new_impute"))
    result <- interval_add_impute(o_data_mixed, target_impute = "new_impute")
    expect_equal(result$intervals[, c("analyte", "half.life", "cmax", "impute")], expected_result)
  })

  it("does not create duplicates but removes the originals & adds impute method based on after", {
    result <- interval_add_impute(o_data, target_impute = "start_conc0", after = Inf)
    expected_result <- data.frame(
      analyte = c("Analyte1", "Analyte2", "Analyte1"),
      half.life = c(TRUE, TRUE, TRUE),
      cmax = c(TRUE, TRUE, TRUE),
      impute = c("start_predose,start_conc0",
                 "start_predose,start_conc0",
                 "start_conc0")
    )
    expect_equal(result$intervals[, c("analyte", "half.life", "cmax", "impute")],
                 expected_result)
  })

  it("adds new rows with added imputations after the original ones", {
    result <- interval_add_impute(o_data, target_impute = "new_impute", target_param = "cmax")
    expected_result <- data.frame(
      analyte = c("Analyte1", "Analyte1", "Analyte2",
                  "Analyte2", "Analyte1", "Analyte1"),
      half.life = c(TRUE, NA, TRUE, NA, TRUE, NA),
      cmax = c(NA, TRUE, NA, TRUE, NA, TRUE),
      impute = c("start_conc0,start_predose",
                 "start_conc0,start_predose,new_impute",
                 "start_predose",
                 "start_predose,new_impute",
                 "start_conc0",
                 "start_conc0,new_impute")
    )
    expect_equal(result$intervals[, c("analyte", "half.life", "cmax", "impute")],
                 expected_result)
  })

  it("does not add new interval if non-target & target params share target impute", {
    intervals_mixed <- data.frame(
      start = c(0, 0),
      end = c(24, 48),
      half.life = c(TRUE, TRUE),
      cmax = c(TRUE, TRUE),
      impute = c("start_conc0,start_predose", "start_predose"),
      analyte = c("Analyte1", "Analyte2"),
      ID = 1
    )

    o_data_mixed <- PKNCA::PKNCAdata(o_conc, o_dose, intervals = intervals_mixed)
    result <- suppressWarnings(interval_add_impute(o_data_mixed,
                                                   target_impute = "start_predose",
                                                   target_param = "cmax",
                                                   after = Inf))
    expected_result <- data.frame(
      analyte = c("Analyte1", "Analyte2"),
      half.life = c(TRUE, TRUE),
      cmax = c(TRUE, TRUE),
      impute = c("start_conc0,start_predose", "start_predose")
    )
    expect_equal(result$intervals[, c("analyte", "half.life", "cmax", "impute")],
                 expected_result)
  })
})

describe("interval_remove_impute", {
  it("removes the impute method in the impute column of a dummy intervals dataframe", {
    simple_df <- data.frame(
      cmax = TRUE,
      impute =  c("", "m0", "m0,m1")
    )
    expected_res <- data.frame(
      cmax = TRUE,
      impute = c("", NA_character_, "m1")
    )
    res <- interval_remove_impute(simple_df, target_impute = "m0")
    expect_equal(res, expected_res)
  })

  it("throws an error if either data or target_impute is missing", {
    expect_error(interval_remove_impute(o_data),
                 "Both 'data' and 'target_impute' must be provided.")
  })

  it("throws an error for non-character target_impute", {
    expect_error(interval_remove_impute(o_data, target_impute = 123),
                 "'target_impute' must be a character string.")
  })

  it("throws an error when input data is not in correct format", {
    expect_error(interval_remove_impute(data = o_conc, target_impute = "start_conc0"))
    expect_no_error(interval_remove_impute(data = o_data, target_impute = "start_conc0"))
  })

  it("throws an error for unknown target_params", {
    expect_error(interval_remove_impute(o_data,
                                        target_impute = "start_conc0",
                                        target_params = "unknown_param"))
  })

  it("handles impute column with FALSE values correctly", {
    o_data_with_na_impute <- o_data
    o_data_with_na_impute$intervals$impute <- NA_character_
    expected_result <- data.frame(analyte = c("Analyte1", "Analyte2", "Analyte1"),
                                  half.life = c(TRUE, TRUE, TRUE),
                                  cmax = c(TRUE, TRUE, TRUE),
                                  impute = c(NA_character_, NA_character_, NA_character_))
    result <- suppressWarnings(
      interval_remove_impute(o_data_with_na_impute, target_impute = "start_conc0")
    )
    expect_equal(result$intervals[, c("analyte", "half.life", "cmax", "impute")], expected_result)
  })

  it("reports an error when impute column is not a character", {
    o_data_not_character_impute <- o_data
    o_data_not_character_impute$intervals$impute <- 1
    expect_error(interval_remove_impute(o_data_not_character_impute, target_impute = "start_conc0"),
                 "The 'impute' column in the intervals data.frame must be a character column.")
  })

  it("warns and makes no changes when target_impute is NA or empty", {
    expect_warning({
      result <- interval_remove_impute(o_data, target_impute = NA_character_)
      expect_equal(result, o_data)
    },
    "No impute method specified. No changes made."
    )

    expect_warning({
      result <- interval_remove_impute(o_data, target_impute = "")
      expect_equal(result, o_data)
    },
    "No impute method specified. No changes made."
    )
  })

  it("does not modify data if global impute & column are missing", {
    d_no_imp <- o_data
    d_no_imp$intervals <- d_no_imp$intervals[, !names(d_no_imp$intervals) %in% "impute"]
    d_no_imp$impute <- NA_character_
    expect_warning({
      res <- interval_remove_impute(d_no_imp, target_impute = "start_conc0")
      expect_equal(res, d_no_imp)
    },
    paste0("No default impute column or global method identified.",
           " No impute methods to remove")
    )

    expect_warning({
      res <- interval_remove_impute(d_no_imp$intervals, target_impute = "start_conc0")
      expect_equal(res, d_no_imp$intervals)
    }, "No default impute column identified. No impute methods to remove")
  })

  it("if impute col is missing uses global impute", {
    o_d_no_imp <- o_data
    o_d_no_imp$intervals <- o_d_no_imp$intervals[, !names(o_d_no_imp$intervals) %in% "impute"]
    o_d_no_imp$impute <- "start_conc0, start_predose"

    # When targets are all intervals, global method is changed
    res_no_target <- interval_remove_impute(o_d_no_imp, target_impute = "start_conc0")
    expect_equal(res_no_target$impute, "start_predose")

    # When targets are specific intervals, then a new column is created and the action handled
    res_target <- interval_remove_impute(o_d_no_imp, target_impute = "start_conc0",
                                         target_groups = data.frame(analyte = "Analyte1"))
    expect_equal(unique(res_target$intervals[res_target$intervals$analyte == "Analyte1", "impute"]),
                 "start_predose")
    expect_equal(res_target$intervals[res_target$intervals$analyte == "Analyte2", "impute"],
                 "start_conc0, start_predose")
  })

  it("with no optional parameters uses all relevant cases", {
    expected_result <- data.frame(analyte = c("Analyte1", "Analyte2", "Analyte1"),
                                  half.life = c(TRUE, TRUE, TRUE),
                                  cmax = c(TRUE, TRUE, TRUE),
                                  impute = c("start_predose", "start_predose", NA_character_))
    result <- interval_remove_impute(o_data, target_impute = "start_conc0")
    expect_equal(result$intervals[, c("analyte", "half.life", "cmax", "impute")], expected_result)
  })

  it("handles specified target_params correctly", {
    expected_result_half_life <- data.frame(
      analyte = c("Analyte1", "Analyte2", "Analyte1"),
      half.life = c(TRUE, TRUE, TRUE),
      impute = c("start_predose", "start_predose", NA_character_)
    )
    expected_result_cmax <- o_data$intervals[o_data$intervals$cmax,
                                             c("analyte", "cmax", "impute")] |>
      `rownames<-`(NULL)
    result <- interval_remove_impute(o_data,
                                     target_impute = "start_conc0",
                                     target_params = "half.life")
    expect_equal(result$intervals[result$intervals$half.life & !is.na(result$intervals$half.life),
                                  c("analyte", "half.life", "impute")] |>
                   `rownames<-`(NULL), expected_result_half_life)
    expect_equal(result$intervals[result$intervals$cmax & !is.na(result$intervals$cmax),
                                  c("analyte", "cmax", "impute")] |>
                   `rownames<-`(NULL), expected_result_cmax)
  })

  it("handles target_groups correctly", {
    expected_result_analyte1 <- data.frame(analyte = c("Analyte1", "Analyte1"),
                                           half.life = c(TRUE, TRUE),
                                           cmax = c(TRUE, TRUE),
                                           impute = c("start_predose", NA_character_))
    expected_result_analyte2 <- o_data$intervals[o_data$intervals$analyte == "Analyte2",
                                                 c("analyte", "half.life", "cmax", "impute")]
    result <- interval_remove_impute(o_data, target_impute = "start_conc0",
                                     target_groups = data.frame(analyte = "Analyte1"))
    expect_equal(result$intervals[result$intervals$analyte == "Analyte1",
                                  c("analyte", "half.life", "cmax", "impute")] |>
                   `rownames<-`(NULL), expected_result_analyte1)
    expect_equal(result$intervals[result$intervals$analyte == "Analyte2",
                                  c("analyte", "half.life", "cmax", "impute")],
                 expected_result_analyte2)
  })

  it("handles multiple target_params correctly", {
    expected_result <- data.frame(analyte = c("Analyte1", "Analyte2", "Analyte1"),
                                  half.life = c(TRUE, TRUE, TRUE),
                                  cmax = c(TRUE, TRUE, TRUE),
                                  impute = c("start_predose", "start_predose", NA_character_))
    result <- interval_remove_impute(
      o_data,
      target_impute = "start_conc0",
      target_params = c("half.life", "cmax")
    )
    expect_equal(result$intervals[, c("analyte", "half.life", "cmax", "impute")], expected_result)
  })

  it("makes no changes and warns when no matching intervals found", {
    expect_warning({
      res <- interval_remove_impute(o_data,
                                    target_impute = "start_conc0",
                                    target_groups = data.frame(analyte = "Analyte3"))
      expect_equal(res, o_data)
      paste0("No intervals found with the specified target parameters,",
             " groups and/or impute method. No changes made.")
    })
  })

  it("handles properly impute character method with multiple imputes", {
    o_data_multiple_imputes <- o_data
    o_data_multiple_imputes$intervals$impute <- "start_conc0,start_predose"
    result <- interval_remove_impute(o_data_multiple_imputes, target_impute = "start_conc0")
    expected_result <- data.frame(
      analyte = c("Analyte1", "Analyte2", "Analyte1"),
      half.life = c(TRUE, TRUE, TRUE),
      cmax = c(TRUE, TRUE, TRUE),
      impute = c("start_predose", "start_predose", "start_predose")
    )
    expect_equal(result$intervals[, c("analyte", "half.life", "cmax", "impute")],
                 expected_result)
  })

  it("handles mixed TRUE/FALSE for cmax and half.life correctly", {
    intervals_mixed <- data.frame(
      start = c(0, 0, 0, 0),
      end = c(24, 48, Inf, 72),
      half.life = c(TRUE, FALSE, TRUE, FALSE),
      cmax = c(FALSE, TRUE, FALSE, TRUE),
      impute = c("start_conc0,start_predose", "start_predose", "start_conc0", "start_predose"),
      analyte = c("Analyte1", "Analyte2", "Analyte1", "Analyte2"),
      ID = c(1, 2, 1, 2)
    )

    o_data_mixed <- PKNCA::PKNCAdata(o_conc, o_dose, intervals = intervals_mixed)

    expected_result <- data.frame(analyte = c("Analyte1", "Analyte2", "Analyte1", "Analyte2"),
                                  half.life = c(TRUE, FALSE, TRUE, FALSE),
                                  cmax = c(FALSE, TRUE, FALSE, TRUE),
                                  impute = c("start_predose", "start_predose",
                                             NA_character_, "start_predose"))
    result <- interval_remove_impute(
      o_data_mixed,
      target_impute = "start_conc0",
      target_params = c("half.life", "cmax")
    )
    expect_equal(result$intervals[, c("analyte", "half.life", "cmax", "impute")], expected_result)
  })

  it("removes all target_impute even if is several times", {
    o_data_multiple_imputes <- o_data
    o_data_multiple_imputes$intervals$impute <- "start_conc0,start_predose,start_conc0"
    result <- interval_remove_impute(o_data_multiple_imputes, target_impute = "start_conc0")
    expected_result <- data.frame(
      analyte = c("Analyte1", "Analyte2", "Analyte1"),
      half.life = c(TRUE, TRUE, TRUE),
      cmax = c(TRUE, TRUE, TRUE),
      impute = c("start_predose", "start_predose", "start_predose")
    )
    expect_equal(result$intervals[, c("analyte", "half.life", "cmax", "impute")],
                 expected_result)
  })

  it("includes new rows right after the original ones", {
    result <- interval_remove_impute(o_data, target_impute = "start_conc0", target_param = "cmax")
    expected_result <- data.frame(
      analyte = c("Analyte1", "Analyte1", "Analyte2", "Analyte1", "Analyte1"),
      half.life = c(TRUE, NA, TRUE, TRUE, NA),
      cmax = c(NA, TRUE, TRUE, NA, TRUE),
      impute = c("start_conc0,start_predose",
                 "start_predose",
                 "start_predose",
                 "start_conc0",
                 NA_character_)
    )
    expect_equal(result$intervals[, c("analyte", "half.life", "cmax", "impute")],
                 expected_result)
  })
})

describe("interval_add_impute and interval_remove_impute", {
  it("are inverses of each other", {
    result_add <- interval_add_impute(o_data, target_impute = "new_impute")
    result_remove <- interval_remove_impute(result_add, target_impute = "new_impute")
    expect_equal(result_remove, o_data)
  })
})

describe("add_impute_method", {
  it("does not crush when impute_vals is empty, returns the empty vector", {
    expect_equal(add_impute_method(character(), "new_impute"), character())
  })
})

describe("remove_impute_method", {
  it("does not crush when impute_vals is empty, returns the empty vector", {
    expect_equal(remove_impute_method(character(), "new_impute"), character())
  })
})
