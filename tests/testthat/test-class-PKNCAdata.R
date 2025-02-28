test_that("PKNCAdata", {
  tmp.conc <- generate.conc(nsub=5, ntreat=2, time.points=0:24)
  tmp.conc.analyte <- generate.conc(nsub=5, ntreat=2, time.points=0:24,
                                    nanalytes=2)
  tmp.conc.study <- generate.conc(nsub=5, ntreat=2, time.points=0:24,
                                  nstudies=2)
  tmp.conc.analyte.study <- generate.conc(nsub=5, ntreat=2, time.points=0:24,
                                          nanalytes=2, nstudies=2)
  tmp.dose <- generate.dose(tmp.conc)
  tmp.dose.analyte <- generate.dose(tmp.conc.analyte)
  tmp.dose.study <- generate.dose(tmp.conc.study)
  tmp.dose.analyte.study <- generate.dose(tmp.conc.analyte.study)
  obj.conc <-
    PKNCAconc(tmp.conc, formula=conc~time|treatment+ID)
  obj.conc.analyte <-
    PKNCAconc(tmp.conc.analyte,
              formula=conc~time|treatment+ID/analyte)
  obj.conc.study <-
    PKNCAconc(tmp.conc.study,
              formula=conc~time|study+treatment+ID)
  obj.conc.analyte.study <-
    PKNCAconc(tmp.conc.analyte.study,
              formula=conc~time|study+treatment+ID/analyte)

  obj.dose <- PKNCAdose(tmp.dose, formula=dose~time|treatment+ID)
  obj.dose.analyte <- PKNCAdose(tmp.dose.analyte, formula=dose~time|treatment+ID)
  obj.dose.study <- PKNCAdose(tmp.dose.study, formula=dose~time|study+treatment+ID)
  obj.dose.analyte.study <- PKNCAdose(tmp.dose.analyte.study, formula=dose~time|study+treatment+ID)

  expect_equal(PKNCAdata(obj.conc, obj.dose),
               PKNCAdata(obj.dose, obj.conc),
               info="Input arguments are reversible")
  expect_equal(PKNCAdata(obj.conc.analyte, obj.dose),
               PKNCAdata(obj.dose, obj.conc.analyte),
               info="Combination of dose and analyte works")
  expect_equal(PKNCAdata(data.conc=tmp.conc, formula.conc=conc~time|treatment+ID,
                         data.dose=tmp.dose, formula.dose=dose~time|treatment+ID),
               PKNCAdata(obj.conc, obj.dose),
               info="Concentration and dose data can be created on the fly")

  # Input checking
  expect_error(PKNCAdata(obj.conc, obj.dose, options="a"),
               regexp="options must be a list.",
               info="Option class")
  expect_error(PKNCAdata(obj.conc, obj.dose, options=list(1)),
               regexp="options must have names.",
               info="Option structure")
  expect_error(PKNCAdata(obj.conc, obj.dose, options=list(foo=1)),
               regexp="Invalid setting for PKNCA.*foo",
               info="Option names")

  # Single dose AUCs are appropriately selected
  expect_equal(
    PKNCAdata(obj.conc, obj.dose),
    {
      tmp.intervals <- tibble::as_tibble(merge(PKNCA.options("single.dose.aucs"), tmp.dose))
      tmp.intervals <- tmp.intervals[order(tmp.intervals$treatment, tmp.intervals$ID),]
      tmp.intervals$time <- NULL
      tmp.intervals$dose <- NULL
      tmp <- list(
        conc=obj.conc,
        dose=obj.dose,
        options=list(),
        intervals=tmp.intervals,
        impute=NA_character_
      )
      class(tmp) <- c("PKNCAdata", "list")
      tmp
    },
    ignore_attr=FALSE,
    info="Selection of single dose AUCs"
  )

  tmp.conc <- generate.conc(nsub=5, ntreat=2, time.points=0:24)
  tmp.dose <- generate.dose(tmp.conc)
  tmp.conc <- tmp.conc[!(tmp.conc$ID %in% 1),]
  obj.conc <-
    PKNCAconc(tmp.conc, formula=conc~time|treatment+ID)
  obj.dose <- PKNCAdose(tmp.dose, formula=dose~time|treatment+ID)
  expect_warning(expect_warning(
    PKNCAdata(obj.conc, obj.dose),
    class = "pknca_no_intervals_generated"),
    class = "pknca_no_intervals_generated",
    info="Missing concentration data with dose data gives a warning."
  )

  expect_warning(expect_warning(expect_warning(
    PKNCAdata(obj.conc, obj.dose, formula.conc=a~b),
    class = "pknca_dataconc_formulaconc"),
    class = "pknca_no_intervals_generated"),
    class = "pknca_no_intervals_generated"
  )
  expect_warning(expect_warning(expect_warning(
    PKNCAdata(obj.conc, obj.dose, formula.dose=a~b),
    class = "pknca_dataconc_formuladose"),
    class = "pknca_no_intervals_generated"),
    class = "pknca_no_intervals_generated"
  )
})

test_that("PKNCAdata with no or limited dose information", {
  tmp.conc <- generate.conc(nsub=5, ntreat=2, time.points=0:24)
  obj.conc <-
    PKNCAconc(tmp.conc, formula=conc~time|treatment+ID)

  expect_error(PKNCAdata(obj.conc),
               regexp="If data.dose is not given, intervals must be given",
               info="One of dose and intervals is required (no dose)")
  expect_error(PKNCAdata(obj.conc, data.dose=NA),
               regexp="If data.dose is not given, intervals must be given",
               info="One of dose and intervals is required (NA dose)")
  expect_equal(
    PKNCAdata(obj.conc, intervals=data.frame(start=0, end=24, aucinf.obs=TRUE)),
    {
      tmp <-
        list(
          conc=obj.conc,
          dose=NA,
          options=list(),
          intervals=check.interval.specification(
            data.frame(start=0, end=24, aucinf.obs=TRUE)),
          impute=NA_character_
        )
      class(tmp) <- c("PKNCAdata", "list")
      tmp
    }
  )

  tmp.conc <- generate.conc(nsub=5, ntreat=2, time.points=0:24)
  tmp.dose <- generate.dose(tmp.conc)
  obj.conc <-
    PKNCAconc(tmp.conc, formula=conc~time|treatment+ID)
  obj.dose <- PKNCAdose(tmp.dose, formula=dose~.|treatment+ID)
  expect_error(PKNCAdata(obj.conc, obj.dose),
               regexp="Dose times were not given, so intervals must be manually specified.",
               info="No dose times requires intervals.")
})

test_that("print.PKNCAdata", {
  tmp.conc <- generate.conc(nsub=2, ntreat=2, time.points=0:24)
  tmp.dose <- generate.dose(tmp.conc)
  obj.conc <- PKNCAconc(tmp.conc, formula=conc~time|treatment+ID)
  obj.dose <- PKNCAdose(tmp.dose, formula=dose~time|treatment+ID)
  obj.data.nodose <- PKNCAdata(obj.conc,
                               intervals=data.frame(start=0, end=24, aucinf.obs=TRUE))
  obj.data.nodose.opt <-
    PKNCAdata(obj.conc,
              intervals=data.frame(start=0, end=24, aucinf.obs=TRUE),
              options=list(min.hl.r.squared=0.95))
  obj.data.dose <- PKNCAdata(obj.conc, data.dose=obj.dose)
  obj.data.units <- PKNCAdata(obj.conc,
                               intervals=data.frame(start=0, end=24, aucinf.obs=TRUE))
  obj.data.units$units <- "mg"

  expect_output(print.PKNCAdata(obj.data.nodose),
                regexp="Formula for concentration:
 conc ~ time | treatment + ID
With 2 subjects defined in the 'ID' column.
Nominal time column is not specified.

First 6 rows of concentration data:
 treatment ID time      conc exclude
     Trt 1  1    0 0.0000000    <NA>
     Trt 1  1    1 0.7052248    <NA>
     Trt 1  1    2 0.7144320    <NA>
     Trt 1  1    3 0.8596094    <NA>
     Trt 1  1    4 0.9998126    <NA>
     Trt 1  1    5 0.7651474    <NA>
No dosing information.

With 1 rows of interval specifications.
No options are set differently than default.",
                info="Generic print.PKNCAdata works with no dosing")
  expect_output(print.PKNCAdata(obj.data.dose),
                regexp="Formula for concentration:
 conc ~ time | treatment + ID
With 2 subjects defined in the 'ID' column.
Nominal time column is not specified.

First 6 rows of concentration data:
 treatment ID time      conc exclude
     Trt 1  1    0 0.0000000    <NA>
     Trt 1  1    1 0.7052248    <NA>
     Trt 1  1    2 0.7144320    <NA>
     Trt 1  1    3 0.8596094    <NA>
     Trt 1  1    4 0.9998126    <NA>
     Trt 1  1    5 0.7651474    <NA>
Formula for dosing:
 dose ~ time | treatment + ID
Nominal time column is not specified.

Data for dosing:
 treatment ID dose time exclude         route duration
     Trt 1  1    1    0    <NA> extravascular        0
     Trt 1  2    1    0    <NA> extravascular        0
     Trt 2  1    2    0    <NA> extravascular        0
     Trt 2  2    2    0    <NA> extravascular        0
With 1 rows of interval specifications.
No options are set differently than default.",
                info="Generic print.PKNCAdata works with dosing")

  expect_output(print.PKNCAdata(obj.data.nodose.opt),
                regexp="Formula for concentration:
 conc ~ time | treatment + ID
With 2 subjects defined in the 'ID' column.
Nominal time column is not specified.

First 6 rows of concentration data:
 treatment ID time      conc exclude
     Trt 1  1    0 0.0000000    <NA>
     Trt 1  1    1 0.7052248    <NA>
     Trt 1  1    2 0.7144320    <NA>
     Trt 1  1    3 0.8596094    <NA>
     Trt 1  1    4 0.9998126    <NA>
     Trt 1  1    5 0.7651474    <NA>
No dosing information.

With 1 rows of interval specifications.
Options changed from default are:
$min.hl.r.squared
[1] 0.95",
                info="Generic print.PKNCAdata works with no dosing and with options changed")
  expect_output(print.PKNCAdata(obj.data.units),
                regexp="Formula for concentration:
 conc ~ time | treatment + ID
With 2 subjects defined in the 'ID' column.
Nominal time column is not specified.

First 6 rows of concentration data:
 treatment ID time      conc exclude
     Trt 1  1    0 0.0000000    <NA>
     Trt 1  1    1 0.7052248    <NA>
     Trt 1  1    2 0.7144320    <NA>
     Trt 1  1    3 0.8596094    <NA>
     Trt 1  1    4 0.9998126    <NA>
     Trt 1  1    5 0.7651474    <NA>
No dosing information.
With units
With 1 rows of interval specifications.
No options are set differently than default.",
                info="Generic print.PKNCAdata works with no dosing")
})

test_that("summary.PKNCAdata", {
  tmp.conc <- generate.conc(nsub=2, ntreat=2, time.points=0:24)
  tmp.dose <- generate.dose(tmp.conc)
  obj.conc <- PKNCAconc(tmp.conc, formula=conc~time|treatment+ID)
  obj.dose <- PKNCAdose(tmp.dose, formula=dose~time|treatment+ID)
  obj.data.nodose <- PKNCAdata(obj.conc,
                               intervals=data.frame(start=0, end=24, aucinf.obs=TRUE))

  expect_output(summary(obj.data.nodose),
                regexp="Formula for concentration:
 conc ~ time | treatment + ID
With 2 subjects defined in the 'ID' column.
Nominal time column is not specified.

Group summary:
 Group Name Count
  treatment     2
         ID     4

First 6 rows of concentration data:
 treatment ID time      conc exclude
     Trt 1  1    0 0.0000000    <NA>
     Trt 1  1    1 0.7052248    <NA>
     Trt 1  1    2 0.7144320    <NA>
     Trt 1  1    3 0.8596094    <NA>
     Trt 1  1    4 0.9998126    <NA>
     Trt 1  1    5 0.7651474    <NA>
No dosing information.

With 1 rows of interval specifications.
No options are set differently than default.",
                info="Generic summary.PKNCAdata works.")
})

test_that("no intervals auto-determined (Fix GitHub issue #84)", {
  tmp_conc <-
    data.frame(
      Subject=1,
      Treatment=c(1, rep(2, 6)),
      Time=c(0, 1:6),
      Conc=1
    )
  tmp_dose <-
    data.frame(
      Subject=1,
      Treatment=c(1, 2),
      Time=c(0, 1),
      Dose=1
    )

  interval_1 <- PKNCA.options("single.dose.aucs")[c(1:2, 1:2),]
  interval_1$start <- rep(0:1, each=2)
  interval_1$end <- c(interval_1$end[1:2], interval_1$end[3:4] + 1)
  interval_1 <- cbind(interval_1, data.frame(Treatment=rep(1:2, each=2), Subject=1))
  two_single_dose_treatments <-
    PKNCAdata(
      PKNCAconc(data=tmp_conc, Conc~Time|Treatment+Subject),
      PKNCAdose(data=tmp_dose, Dose~Time|Treatment+Subject)
    )
  expect_equal(
    two_single_dose_treatments$intervals,
    interval_1,
    ignore_attr=TRUE
  )
  interval_2 <-
    check.interval.specification(
      tibble::tibble(
        start=1, end=c(2, Inf),
        auclast=c(TRUE, FALSE),
        cmax=c(TRUE, FALSE),
        tmax=c(TRUE, FALSE),
        half.life=c(FALSE, TRUE),
        Treatment=2,
        Subject=1
      )
    )
  expect_warning(
    two_multiple_dose_treatments <-
      PKNCAdata(
        PKNCAconc(data=tmp_conc, Conc~Time|Treatment+Subject),
        PKNCAdose(data=tmp_dose, Dose~Time|Subject)
      ),
    regexp="No intervals generated"
  )
  expect_equal(
    two_multiple_dose_treatments$intervals,
    interval_2
  )
})

test_that("Ensure that unexpected arguments to PKNCAdata give an error (related to issue #83)", {
  tmp.conc <- generate.conc(nsub=2, ntreat=1, time.points=0:24)
  tmp.dose <- generate.dose(tmp.conc)
  obj.conc <-
    PKNCAconc(tmp.conc, formula=conc~time|treatment+ID)
  obj.dose <-
    PKNCAdose(tmp.dose, formula=dose~time|treatment+ID)
  expect_error(mydata <- PKNCAdata(obj.conc, obj.dose, 1),
               regexp="Unknown argument")
})

test_that("intervals may be a tibble", {
  tmp.conc <- generate.conc(nsub=2, ntreat=1, time.points=0:24)
  tmp.dose <- generate.dose(tmp.conc)
  obj.conc <-
    PKNCAconc(tmp.conc, formula=conc~time|treatment+ID)
  obj.dose <-
    PKNCAdose(tmp.dose, formula=dose~time|treatment+ID)
  intervals <- data.frame(start=0, end=24, aucinf.obs=TRUE)
  mydata_tibble <- PKNCAdata(obj.conc, obj.dose, intervals=dplyr::as_tibble(intervals))
  mydata <- PKNCAdata(obj.conc, obj.dose, intervals=intervals)
  expect_equal(
    as.data.frame(pk.nca(mydata_tibble)),
    as.data.frame(pk.nca(mydata))
  )
})

test_that("PKNCAdata units (#336)", {
  # Typical use
  d_conc <- data.frame(conc = 1, time = 0, concu_x = "A", timeu_x = "B", amountu_x = "C")
  d_dose <- data.frame(dose = 1, time = 0, doseu_x = "D")

  o_conc <- PKNCAconc(data = d_conc, conc~time, concu = "concu_x", timeu = "timeu_x")
  o_dose <- PKNCAdose(data = d_dose, dose~time, doseu = "doseu_x")
  o_data <- PKNCAdata(o_conc, o_dose)
  expect_equal(
    o_data$units,
    pknca_units_table(concu = "A", doseu = "D", timeu = "B")
  )
  suppressWarnings(o_nca <- pk.nca(o_data))
  expect_true("Cmax (A)" %in% names(summary(o_nca)))

  # NA unit values are ignored
  d_conc <- data.frame(conc = 1, time = 0:1, concu_x = c("A", NA), timeu_x = "B", amountu_x = "C")
  d_dose <- data.frame(dose = 1, time = 0, doseu_x = "D")

  o_conc <- PKNCAconc(data = d_conc, conc~time, concu = "concu_x", timeu = "timeu_x")
  o_dose <- PKNCAdose(data = d_dose, dose~time, doseu = "doseu_x")
  o_data <- PKNCAdata(o_conc, o_dose)
  expect_equal(
    o_data$units,
    pknca_units_table(concu = "A", doseu = "D", timeu = "B")
  )
  suppressWarnings(o_nca <- pk.nca(o_data))
  expect_true("Cmax (A)" %in% names(summary(o_nca)))

  # multiple unit values cause an error
  d_conc <- data.frame(conc = 1, time = 0:1, concu_x = c("A", "C"), timeu_x = "B", amountu_x = "C")
  d_dose <- data.frame(dose = 1, time = 0, doseu_x = "B")

  o_conc <- PKNCAconc(data = d_conc, conc~time, concu = "concu_x")
  o_dose <- PKNCAdose(data = d_dose, dose~time, doseu = "doseu_x")
  expect_error(
    PKNCAdata(o_conc, o_dose),
    regexp = "Only one unit may be provided at a time: A, C"
  )
})

test_that("getGroups works", {
  # Check that it works with grouping [contains only the grouping column(s)]
  o_conc_group <- PKNCAconc(as.data.frame(datasets::Theoph), conc~Time|Subject)
  data_group <- as.data.frame(datasets::Theoph)
  expected_group <- data.frame(Subject = data_group$Subject)
  o_data_group <- PKNCAdata(o_conc_group, intervals = data.frame(start = 0, end = 1, cmax = TRUE))
  expect_equal(getGroups(o_data_group), expected_group)

  # Check that it works without groupings as expected [empty]
  o_conc_nongroup <- PKNCAconc(as.data.frame(datasets::Theoph)[datasets::Theoph$Subject == 1,], conc~Time)
  o_data_nogroup <- PKNCAdata(o_conc_nongroup, intervals = data.frame(start = 0, end = 1, cmax = TRUE))

  # It should be an empty data.frame with 11 rows
  expect_equal(getGroups(o_data_nogroup), data.frame(A = 1:11)[, -1])
})

test_that("group_vars.PKNCAdata", {
  o_conc_group <- PKNCAconc(as.data.frame(datasets::Theoph), conc~Time|Subject)
  o_data_group <- PKNCAdata(o_conc_group, intervals = data.frame(start = 0, end = 1, cmax = TRUE))

  expect_equal(dplyr::group_vars(o_data_group), "Subject")

  # Check that it works without groupings as expected [empty]
  o_conc_nongroup <- PKNCAconc(as.data.frame(datasets::Theoph)[datasets::Theoph$Subject == 1,], conc~Time)
  o_data_nogroup <- PKNCAdata(o_conc_nongroup, intervals = data.frame(start = 0, end = 1, cmax = TRUE))

  expect_equal(dplyr::group_vars(o_data_nogroup), character(0))
})
