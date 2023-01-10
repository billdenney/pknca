source("generate.data.R")

test_that("PKNCAconc expected errors", {
  tmp.conc <- generate.conc(nsub=1, ntreat=1, time.points=0:24)
  tmp.conc$foo <- "A"
  expect_error(
    PKNCAconc(conc~time, volume="foo", data=tmp.conc),
    regexp="Volume must be numeric"
  )
  expect_error(
    PKNCAconc(conc~time, duration="foo", data=tmp.conc),
    regexp="duration must be numeric without missing (NA) or infinite values, and all values must be >= 0",
    fixed=TRUE
  )
})

test_that("PKNCAconc", {
  tmp.conc <- generate.conc(nsub=5, ntreat=2, time.points=0:24)
  tmp.conc.analyte <- generate.conc(nsub=5, ntreat=2, time.points=0:24,
                                    nanalytes=2)
  tmp.conc.study <- generate.conc(nsub=5, ntreat=2, time.points=0:24,
                                  nstudies=2)
  tmp.conc.analyte.study <- generate.conc(nsub=5, ntreat=2, time.points=0:24,
                                          nanalytes=2, nstudies=2)
  # Data exists
  expect_error(PKNCAconc(data.frame()),
               regexp="data must have at least one row.",
               info="PKNCAconc requires data")

  # Variables present
  expect_error(
    PKNCAconc(tmp.conc, formula=XXX~time|treatment+ID),
    regexp="All of the variables in the formula must be in the data.  Missing: XXX",
    info="All formula parameters must be in the data (LHS)"
  )
  expect_error(
    PKNCAconc(tmp.conc, formula=conc~XXX|treatment+ID),
    regexp="All of the variables in the formula must be in the data.  Missing: XXX",
    info="All formula parameters must be in the data (RHS)"
  )
  expect_error(
    PKNCAconc(tmp.conc, formula=conc~time|XXX+ID),
    regexp="All of the variables in the formula must be in the data.  Missing: XXX",
    info="All formula parameters must be in the data (groups)"
  )

  # Number of variables
  expect_error(PKNCAconc(tmp.conc, formula=conc+ID~time|treatment+ID),
               regexp="The left hand side of the formula must have exactly one variable",
               info="The right number of parameters in the formula (LHS)")
  expect_error(PKNCAconc(tmp.conc, formula=conc~time+ID|treatment+ID),
               regexp="The right hand side of the formula \\(excluding groups\\) must have exactly one variable",
               info="The right number of parameters in the formula (RHS)")

  # Subject assignment
  expect_equal(PKNCAconc(tmp.conc.analyte, formula=conc~time|treatment+ID/analyte),
               PKNCAconc(tmp.conc.analyte, formula=conc~time|treatment+ID/analyte, subject="ID"))
  expect_error(PKNCAconc(tmp.conc.analyte, formula=conc~time|treatment+ID/analyte, subject=5),
               regexp="subject must be a character string")
  expect_error(PKNCAconc(tmp.conc.analyte, formula=conc~time|treatment+ID/analyte, subject=c("", "foo")),
               regexp="subject must be a scalar")
  expect_error(PKNCAconc(tmp.conc.analyte, formula=conc~time|treatment+ID/analyte, subject="foo"),
               regexp="The subject parameter must map to a name in the data")

  # Keys must be unique
  expect_error(PKNCAconc(tmp.conc.analyte, formula=conc~time|treatment+ID),
               regexp="Rows that are not unique per group and time",
               info="Duplicated key rows")

  expect_equal(
    PKNCAconc(
      tmp.conc.analyte,
      formula=conc~time|treatment+ID/analyte
    ),
    PKNCAconc(
      dplyr::as_tibble(tmp.conc.analyte),
      formula=conc~time|treatment+ID/analyte
    ),
    info="tibble and data.frame classes both work and create identical objects"
  )
})

test_that("PKNCAconc with input other than data.frames", {
  tmp <- structure(list(), class="foo")
  captured_message <- tryCatch(as.data.frame(tmp), error=function(e) e$message)
  expect_error(PKNCAconc(tmp, formula=conc~time|treatment+ID),
               regexp=captured_message,
               fixed=TRUE,
               info="Attempt to coerce into a data.frame.")
})

test_that("model frame and parameter extractions", {
  tmp.conc <- generate.conc(nsub=5, ntreat=2, time.points=0:24)
  tmp_conc_single <- generate.conc(nsub=1, ntreat=1, time.points=0:24)

  expect_equal(model.frame(PKNCAconc(tmp.conc, formula=conc~time|treatment+ID)),
               tmp.conc[,c("conc", "time", "treatment", "ID")],
               info="model.frame.PKNCAconc extracts the correct components")
  expect_equal(getDepVar(PKNCAconc(tmp.conc, formula=conc~time|treatment+ID)),
               tmp.conc$conc,
               info="getDepVar.PKNCAconc extracts the correct component")
  expect_equal(getIndepVar(PKNCAconc(tmp.conc, formula=conc~time|treatment+ID)),
               tmp.conc$time,
               info="getIndepVar.PKNCAconc extracts the correct component")
  expect_equal(getGroups.PKNCAconc(PKNCAconc(tmp.conc, formula=conc~time|treatment+ID)),
               tmp.conc[,c("treatment", "ID")],
               info="getGroups.PKNCAconc extracts the correct components")

  tmp.conc <- generate.conc(nsub=5, ntreat=1, time.points=0:24)
  expect_equal(getGroups.PKNCAconc(PKNCAconc(tmp.conc, formula=conc~time|ID)),
               tmp.conc[,"ID", drop=FALSE],
               info="getGroups.PKNCAconc returns a data.frame even with a single grouping level")
  # Levels
  tmp.conc <- generate.conc(nsub=5, ntreat=2, time.points=0:24)
  expect_equal(getGroups.PKNCAconc(PKNCAconc(tmp.conc, formula=conc~time|treatment+ID), level=1),
               tmp.conc[,"treatment", drop=FALSE],
               info="getGroups.PKNCAconc the correct level (numeric scalar)")
  expect_equal(getGroups.PKNCAconc(PKNCAconc(tmp.conc, formula=conc~time|treatment+ID), level=-1),
               tmp.conc[,"ID", drop=FALSE],
               info="getGroups.PKNCAconc the correct level (negative numeric scalar)")
  expect_equal(getGroups.PKNCAconc(PKNCAconc(tmp.conc, formula=conc~time|treatment+ID), level=1:2),
               tmp.conc[,c("treatment", "ID"), drop=FALSE],
               info="getGroups.PKNCAconc the correct level (numeric vector)")
  expect_equal(getGroups.PKNCAconc(PKNCAconc(tmp.conc, formula=conc~time|treatment+ID), level=2:1),
               tmp.conc[,c("ID", "treatment"), drop=FALSE],
               info="getGroups.PKNCAconc the correct level (numeric vector)")
  expect_equal(getGroups.PKNCAconc(PKNCAconc(tmp.conc, formula=conc~time|treatment+ID), level="ID"),
               tmp.conc[,"ID", drop=FALSE],
               info="getGroups.PKNCAconc the correct level (character string scalar)")
  expect_equal(getGroups.PKNCAconc(PKNCAconc(tmp.conc, formula=conc~time|treatment+ID), level=c("ID", "treatment")),
               tmp.conc[,c("ID", "treatment"), drop=FALSE],
               info="getGroups.PKNCAconc the correct level (character string vector)")
  expect_error(getGroups.PKNCAconc(PKNCAconc(tmp.conc, formula=conc~time|treatment+ID), level="foo"),
               regexp="Not all levels are listed in the group names",
               info="getGroups.PKNCAconc gives an error if a group name is not present")

  expect_equal(
    group_vars.PKNCAconc(PKNCAconc(tmp.conc, formula=conc~time|treatment+ID)),
    c("treatment", "ID")
  )
  expect_equal(
    group_vars.PKNCAconc(PKNCAconc(tmp_conc_single, formula=conc~time)),
    character(0),
    info="Ungrouped data works with group_vars"
  )
})

test_that("print.PKNCAconc", {
  tmp.conc <- generate.conc(nsub=2, ntreat=2, time.points=0:24)
  myconc <- PKNCAconc(tmp.conc, formula=conc~time|treatment+ID)

  expect_output(print.PKNCAconc(myconc),
                regexp="Formula for concentration:
 conc ~ time | treatment + ID
                With 2 subjects defined in the 'ID' column.
                Nominal time column is not specified.

                First 6 rows of concentration data:
                treatment ID time      conc
                Trt 1  1    0 0.0000000
                Trt 1  1    1 0.7052248
                Trt 1  1    2 0.7144320
                Trt 1  1    3 0.8596094
                Trt 1  1    4 0.9998126
                Trt 1  1    5 0.7651474",
                info="Generic print.PKNCAconc works")
  expect_output(print.PKNCAconc(myconc, n=0),
                regexp="Formula for concentration:
 conc ~ time | treatment + ID
With 2 subjects defined in the 'ID' column.
Nominal time column is not specified.",
                info="print.PKNCAconc respects the n argument.")
  expect_output(print.PKNCAconc(myconc, n=-98),
                regexp="Formula for concentration:
 conc ~ time | treatment + ID
                With 2 subjects defined in the 'ID' column.
                Nominal time column is not specified.

                First 2 rows of concentration data:
                treatment ID time      conc
                Trt 1  1    0 0.0000000
                Trt 1  1    1 0.7052248",
                info="print.PKNCAconc accurately uses negative n argument.")

  tmp.conc <- generate.conc(nsub=1, ntreat=1, time.points=0:24)
  myconc <- PKNCAconc(tmp.conc, formula=conc~time)
  expect_output(
    print(myconc),
    regexp="As a single-subject dataset"
  )
  expect_output(
    print(myconc, summarize=TRUE),
    regexp="No groups\\."
  )
  expect_output(
    print(myconc, n=1e6),
    regexp="Data for concentration"
  )

  myconc <- PKNCAconc(tmp.conc, formula=conc~time|treatment+ID, time.nominal="time")
  expect_output(
    print(myconc),
    regexp = "Nominal time column is: time"
  )
})

test_that("summary.PKNCAconc", {
  tmp.conc <- generate.conc(nsub=2, ntreat=2, time.points=0:24)
  myconc <- PKNCAconc(tmp.conc, formula=conc~time|treatment+ID)

  expect_output(summary(myconc),
                regexp="Formula for concentration:
 conc ~ time | treatment + ID
                With 2 subjects defined in the 'ID' column.
                Nominal time column is not specified.

                Group summary:
                Group Name Count
                treatment     2
                ID     4",
                info="Generic summary.PKNCAconc works.")
})

test_that("PKNCAconc with exclusions", {
  tmp.conc <- generate.conc(nsub=2, ntreat=2, time.points=0:24)
  tmp.conc$excl <- NA_character_
  myconc <- PKNCAconc(tmp.conc, formula=conc~time|treatment+ID, exclude="excl")
  expect_equal(
    myconc,
    structure(
      list(data=cbind(tmp.conc,
                      volume=NA_real_,
                      duration=0),
           formula=conc~time|treatment+ID,
           columns=
             list(
               concentration="conc",
               time="time",
               groups=
                 list(
                   group_vars=c("treatment", "ID"),
                   group_analyte=character()
                 ),
               subject="ID",
               exclude="excl",
               volume="volume",
               duration="duration"
             )
      ),
      class=c("PKNCAconc", "list")
    )
  )
})

test_that("PKNCAconc with duration", {
  tmp.conc <- generate.conc(nsub=2, ntreat=2, time.points=0:24)
  tmp.conc$duration_test <- 0.1
  myconc <- PKNCAconc(tmp.conc,
                      formula=conc~time|treatment+ID,
                      duration="duration_test")
  expect_equal(
    myconc,
    structure(
      list(data=cbind(tmp.conc,
                      data.frame(exclude=NA_character_,
                                 volume=NA_real_,
                                 stringsAsFactors=FALSE)),
           formula=conc~time|treatment+ID,
           columns=
             list(
               concentration="conc",
               time="time",
               groups=
                 list(
                   group_vars=c("treatment", "ID"),
                   group_analyte=character()
                 ),
               subject="ID",
               exclude="exclude",
               volume="volume",
               duration="duration_test"
             )
      ),
      class=c("PKNCAconc", "list")
    )
  )
})

test_that("PKNCAconc with nominal time added", {
  tmp.conc <- generate.conc(nsub=2, ntreat=2, time.points=0:24)
  tmp.conc$tnom <- tmp.conc$time
  myconc <- PKNCAconc(tmp.conc, formula=conc~time|treatment+ID, time.nominal="tnom")
  expect_equal(
    myconc,
    structure(
      list(data=cbind(tmp.conc,
                      data.frame(exclude=NA_character_,
                                 volume=NA_real_,
                                 duration=0,
                                 stringsAsFactors=FALSE)),
           formula=conc~time|treatment+ID,
           columns=
             list(
               concentration="conc",
               time="time",
               groups=
                 list(
                   group_vars=c("treatment", "ID"),
                   group_analyte=character()
                 ),
               subject="ID",
               exclude="exclude",
               volume="volume",
               duration="duration",
               time.nominal="tnom")),
      class=c("PKNCAconc", "list")
    )
  )
  expect_equal(
    PKNCAconc(tmp.conc, formula=conc~time|treatment+ID, time.nominal="foo"),
    structure(
      list(data=cbind(tmp.conc,
                      data.frame(exclude=NA_character_,
                                 volume=NA_real_,
                                 duration=0,
                                 foo=NA,
                                 stringsAsFactors=FALSE)),
           formula=conc~time|treatment+ID,
           columns=
             list(
               concentration="conc",
               time="time",
               groups=
                 list(
                   group_vars=c("treatment", "ID"),
                   group_analyte=character()
                 ),
               subject="ID",
               exclude="exclude",
               volume="volume",
               duration="duration",
               time.nominal="foo"
             )
      ),
      class=c("PKNCAconc", "list"))
  )
})

test_that("PKNCAconc with volume added", {
  tmp.conc <- generate.conc(nsub=2, ntreat=2, time.points=0:24)
  tmp.conc$vol <- seq_len(nrow(tmp.conc))
  myconc <- PKNCAconc(tmp.conc, formula=conc~time|treatment+ID, volume="vol")
  expect_equal(
    myconc,
    structure(
      list(data=cbind(tmp.conc,
                      data.frame(exclude=NA_character_,
                                 duration=0,
                                 stringsAsFactors=FALSE)),
           formula=conc~time|treatment+ID,
           columns=list(
             concentration="conc",
             time="time",
             groups=
               list(
                 group_vars=c("treatment", "ID"),
                 group_analyte=character()
               ),
             subject="ID",
             exclude="exclude",
             volume="vol",
             duration="duration"
           )
      ),
      class=c("PKNCAconc", "list"))
  )
  myconc_manual_vol <- PKNCAconc(tmp.conc, formula=conc~time|treatment+ID, volume=2)
  expect_equal(
    myconc_manual_vol,
    structure(
      list(data=cbind(tmp.conc,
                      data.frame(exclude=NA_character_,
                                 volume=2,
                                 duration=0,
                                 stringsAsFactors=FALSE)),
           formula=conc~time|treatment+ID,
           columns=
             list(
               concentration="conc",
               time="time",
               groups=
                 list(
                   group_vars=c("treatment", "ID"),
                   group_analyte=character()
                 ),
               subject="ID",
               exclude="exclude",
               volume="volume",
               duration="duration"
             )
      ),
      class=c("PKNCAconc", "list")
    )
  )
  myconc_manual_vol_vector <- PKNCAconc(tmp.conc, formula=conc~time|treatment+ID, volume=seq_len(nrow(tmp.conc)))
  expect_equal(
    myconc_manual_vol_vector,
    structure(
      list(
        data=
          cbind(
            tmp.conc,
            data.frame(
              exclude=NA_character_,
              volume=seq_len(nrow(tmp.conc)),
              duration=0
            )
          ),
        formula=conc~time|treatment+ID,
        columns=
          list(
            concentration="conc",
            time="time",
            groups=
              list(
                group_vars=c("treatment", "ID"),
                group_analyte=character()
              ),
            subject="ID",
            exclude="exclude",
            volume="volume",
            duration="duration"
          )
      ),
      class=c("PKNCAconc", "list")
    )
  )
})

test_that("as.data.frame.PKNCAconc", {
  tmp_conc <- generate.conc(nsub=1, ntreat=1, time.points=0:24)
  result_conc <- tmp_conc
  result_conc$exclude <- NA_character_
  result_conc$volume <- NA_real_
  result_conc$duration <- 0
  expect_equal(
    as.data.frame(PKNCAconc(conc~time, data=tmp_conc)),
    result_conc
  )
})

test_that("PKNCAconc with sparse data", {
  d_sparse <-
    data.frame(
      id = c(1L, 2L, 3L, 1L, 2L, 3L, 1L, 2L, 3L, 4L, 5L, 6L, 4L, 5L, 6L, 7L, 8L, 9L, 7L, 8L, 9L),
      conc = c(0, 0, 0,  1.75, 2.2, 1.58, 4.63, 2.99, 1.52, 3.03, 1.98, 2.22, 3.34, 1.3, 1.22, 3.54, 2.84, 2.55, 0.3, 0.0421, 0.231),
      time = c(0, 0, 0, 1, 1, 1, 6, 6, 6, 2, 2, 2, 10, 10, 10, 4, 4, 4, 24, 24, 24),
      dose = c(100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100)
    )
  o_conc_sparse <- PKNCAconc(d_sparse, conc~time|id, sparse=TRUE)
  expect_true("data_sparse" %in% names(o_conc_sparse))
  expect_false("data" %in% names(o_conc_sparse))

  d_sparse_aug <- d_sparse
  d_sparse_aug$exclude <- NA_character_
  d_sparse_aug$volume <- NA_real_
  d_sparse_aug$duration <- 0
  expect_equal(
    o_conc_sparse$data_sparse,
    d_sparse_aug
  )
})

test_that("print.PKNCAconc with sparse data", {
  d_sparse <-
    data.frame(
      id = c(1L, 2L, 3L, 1L, 2L, 3L, 1L, 2L, 3L, 4L, 5L, 6L, 4L, 5L, 6L, 7L, 8L, 9L, 7L, 8L, 9L),
      conc = c(0, 0, 0,  1.75, 2.2, 1.58, 4.63, 2.99, 1.52, 3.03, 1.98, 2.22, 3.34, 1.3, 1.22, 3.54, 2.84, 2.55, 0.3, 0.0421, 0.231),
      time = c(0, 0, 0, 1, 1, 1, 6, 6, 6, 2, 2, 2, 10, 10, 10, 4, 4, 4, 24, 24, 24),
      dose = c(100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100)
    )
  o_conc_sparse <- PKNCAconc(d_sparse, conc~time|id, sparse=TRUE)
  expect_output(print.PKNCAconc(o_conc_sparse),
                regexp="Formula for concentration:
 conc ~ time | id
Data are sparse PK.
With 9 subjects defined in the 'id' column.
Nominal time column is not specified.

First 6 rows of concentration data:
 id conc time dose exclude volume duration
  1 0.00    0  100    <NA>     NA        0
  2 0.00    0  100    <NA>     NA        0
  3 0.00    0  100    <NA>     NA        0
  1 1.75    1  100    <NA>     NA        0
  2 2.20    1  100    <NA>     NA        0
  3 1.58    1  100    <NA>     NA        0"
  )
})
