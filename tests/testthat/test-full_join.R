library(dplyr)
source("generate.data.R")

test_that("nest", {
  tmp_conc <- generate.conc(nsub=5, ntreat=2, time.points=0:24)
  tmp_dose <- generate.dose(tmp_conc)
  o_conc <- PKNCAconc(tmp_conc, formula=conc~time|treatment+ID)
  o_dose <- PKNCAdose(tmp_dose, formula=dose~time|treatment+ID)

  expect_equal(
    nest(o_conc),
    tidyr::nest(
      dplyr::mutate(
        tmp_conc,
        exclude=NA_character_,
        volume=NA_real_,
        duration=0
      ),
      data_conc=!c("treatment", "ID")
    )
  )
  expect_equal(
    nest(o_dose),
    tidyr::nest(
      dplyr::mutate(
        tmp_dose,
        exclude=NA_character_,
        route="extravascular",
        duration=0
      ),
      data_dose=!c("treatment", "ID")
    )
  )
  # No groups
  expect_equal(
    nest_PKNCAintervals(.dat=PKNCA.options("single.dose.aucs")),
    tibble(
      data_intervals=list(as_tibble(PKNCA.options("single.dose.aucs")))
    )
  )
  # With groups
  tmp_intervals <- PKNCA.options("single.dose.aucs")
  tmp_intervals$g <- "A"
  expect_equal(
    nest_PKNCAintervals(.dat=tmp_intervals, vars="g"),
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
    dplyr::full_join(nest(o_conc), nest(o_dose), by=c("treatment", "ID"))
  )
  expect_equal(
    full_join_PKNCAdata(o_data),
    tidyr::crossing(
      dplyr::full_join(
        nest(o_conc), nest(o_dose),
        by=c("treatment", "ID")
      ),
      data_intervals=list(as_tibble(PKNCA.options("single.dose.aucs")))
    )
  )
  # When intervals have no groups
  o_data_manual_interval <- PKNCAdata(o_conc, o_dose, intervals=PKNCA.options("single.dose.aucs")[1,])
  expect_equal(
    full_join_PKNCAdata(o_data_manual_interval),
    tidyr::crossing(
      dplyr::full_join(
        nest(o_conc), nest(o_dose),
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
      nest(o_conc),
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
