library(dplyr)
source("generate.data.R")

test_that("prepare_*", {
  tmp_conc <- generate.conc(nsub=5, ntreat=2, time.points=0:24)
  tmp_dose <- generate.dose(tmp_conc)
  o_conc <- PKNCAconc(tmp_conc, formula=conc~time|treatment+ID)
  o_dose <- PKNCAdose(tmp_dose, formula=dose~time|treatment+ID)

  expect_equal(
    prepare_PKNCAconc(o_conc),
    tidyr::nest(
      dplyr::mutate(
        tmp_conc[, c("treatment", "ID", "conc", "time")],
        volume=NA_real_,
        duration=0
      ),
      data_conc=!c("treatment", "ID")
    )
  )
  expect_equal(
    prepare_PKNCAdose(o_dose),
    tidyr::nest(
      dplyr::mutate(
        tmp_dose,
        duration=0,
        route="extravascular"
      ),
      data_dose=!c("treatment", "ID")
    )
  )
  # No groups
  expect_equal(
    prepare_PKNCAintervals(.dat=PKNCA.options("single.dose.aucs")),
    tibble(
      data_intervals=list(as_tibble(PKNCA.options("single.dose.aucs")))
    )
  )
  # With groups
  tmp_intervals <- PKNCA.options("single.dose.aucs")
  tmp_intervals$g <- "A"
  expect_equal(
    prepare_PKNCAintervals(.dat=tmp_intervals, vars="g"),
    tibble::tibble(
      g="A",
      data_intervals=list(tibble::as_tibble(PKNCA.options("single.dose.aucs")))
    )
  )
})

test_that("full_join for PKNCAconc, PKNCAdose, and PKNCAdata", {
  tmp_conc <- generate.conc(nsub=5, ntreat=2, time.points=0:24)
  tmp_dose <- generate.dose(tmp_conc)
  o_conc <- PKNCAconc(tmp_conc, formula=conc~time|treatment+ID)
  o_dose <- PKNCAdose(tmp_dose, formula=dose~time|treatment+ID)
  o_data <- PKNCAdata(o_conc, o_dose)

  expect_equal(
    full_join_PKNCAconc_PKNCAdose(o_conc, o_dose),
    dplyr::full_join(
      prepare_PKNCAconc(o_conc),
      prepare_PKNCAdose(o_dose),
      by=c("treatment", "ID")
    )
  )
  expect_equal(
    full_join_PKNCAdata(o_data),
    tidyr::crossing(
      dplyr::full_join(
        prepare_PKNCAconc(o_conc),
        prepare_PKNCAdose(o_dose),
        by=c("treatment", "ID")
      ),
      data_intervals=list(as_tibble(PKNCA.options("single.dose.aucs")))
    )
  )
  # When intervals have no groups
  o_data_manual_interval <-
    PKNCAdata(
      o_conc,
      o_dose,
      intervals=PKNCA.options("single.dose.aucs")[1,]
    )
  expect_equal(
    full_join_PKNCAdata(o_data_manual_interval),
    tidyr::crossing(
      dplyr::full_join(
        prepare_PKNCAconc(o_conc),
        prepare_PKNCAdose(o_dose),
        by=c("treatment", "ID")
      ),
      data_intervals=list(as_tibble(PKNCA.options("single.dose.aucs")[1,]))
    )
  )
  # When dosing is not provided
  o_data_no_dose <- PKNCAdata(o_conc, intervals=PKNCA.options("single.dose.aucs")[1,])
  expect_equal(
    full_join_PKNCAdata(o_data_no_dose),
    tidyr::crossing(
      prepare_PKNCAconc(o_conc),
      tibble(data_dose=list(NA)),
      data_intervals=list(as_tibble(PKNCA.options("single.dose.aucs")[1,]))
    )
  )
})

test_that("check_reserved_column_names", {
  expect_null(
    check_reserved_column_names(data.frame())
  )
  expect_error(
    check_reserved_column_names(data.frame(data_dose=1)),
    regexp="The column 'data_dose' is reserved for internal use in PKNCA.  Change the name and retry.",
    fixed=TRUE
  )
  expect_error(
    check_reserved_column_names(data.frame(data_dose=1, data_conc=1, data_intervals=1)),
    regexp="The columns 'data_conc', 'data_dose', 'data_intervals' are reserved for internal use in PKNCA.  Change the names and retry.",
    fixed=TRUE
  )
})

test_that("standardize_column_names", {
  # One column works
  expect_equal(
    standardize_column_names(data.frame(a=1), cols=list(b="a")),
    data.frame(b=1)
  )
  # Two columns work
  expect_equal(
    standardize_column_names(data.frame(a=1, b=2), cols=list(c="a", d="b")),
    data.frame(c=1, d=2)
  )
  # group_cols overlap with cols values fails
  expect_error(
    standardize_column_names(data.frame(a=1, b=2), cols=list(c="a", d="b"), group_cols="b"),
    regexp="group_cols must not overlap with other column names"
  )
  # group_cols overlap with cols names fails
  expect_error(
    standardize_column_names(data.frame(a=1, b=2), cols=list(c="a", d="b"), group_cols="c"),
    regexp="group_cols must not overlap with standardized column names"
  )
  # group_cols works
  expect_equal(
    standardize_column_names(data.frame(a=1, b=2), cols=list(d="b"), group_cols="a"),
    data.frame(group1=1, d=2)
  )
  # Missing values are inserted correctly
  expect_equal(
    standardize_column_names(
      data.frame(a=1, b=2),
      cols=list(d="b"),
      group_cols="a",
      insert_if_missing=list(dose=1)
    ),
    data.frame(group1=1, d=2, dose=1)
  )
  # Missing values are not inserted when something is already in the data for
  # that column
  expect_equal(
    standardize_column_names(
      data.frame(a=1, b=2, c=3),
      cols=list(d="b", dose="c"),
      group_cols="a",
      insert_if_missing=list(dose=2)
    ),
    data.frame(group1=1, d=2, dose=3)
  )
})

test_that("restore_group_col_names", {
  d <- data.frame(a=1, group1=2)
  # zero group cols
  expect_equal(
    restore_group_col_names(d),
    d
  )
  # One group col
  expect_equal(
    restore_group_col_names(d, group_cols="b"),
    data.frame(a=1, b=2)
  )
  # More than one group col
  expect_equal(
    restore_group_col_names(data.frame(a=1, group1=2, group2=3), c("b", "c")),
    data.frame(a=1, b=2, c=3)
  )
  expect_error(
    restore_group_col_names(d, group_cols=c("b", "c")),
    regexp="missing intermediate group_cols names"
  )
  expect_error(
    restore_group_col_names(data.frame(a=1, group2=3, group1=2), c("b", "c")),
    regexp="Intermediate group_cols are out of order"
  )
})
