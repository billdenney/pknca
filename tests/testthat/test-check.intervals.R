context("AUC interval checking")

test_that(
  "check.interval.specification", {

    ## Get the current name order of the expected results
    nameorder <- names(check.interval.deps(data.frame(start=0, end=1, cmax=TRUE)))

    d1 <- data.frame(start=0, end=1)
    r1 <- data.frame(start=0,
                     end=1,
                     stringsAsFactors=FALSE)
    r1[,setdiff(nameorder, names(r1))] <- FALSE
    expect_warning(check.interval.specification(as.matrix(d1)),
                   regexp="Interval specification must be a data.frame",
                   info="Interval must be a data.frame or coercable into a data frame")
    expect_warning(d1.check <- check.interval.specification(d1),
                   regexp="Nothing to be calculated in interval specification number\\(s\\): 1",
                   info="Warn if nothing is to be calculated in an interval specification")
    expect_equal(d1.check, r1[,nameorder],
                 info="Expand a minimal data frame for interval specification")

    ## Giving one parameter will fill in everything else as false
    d2 <- data.frame(start=0, end=1, auclast=TRUE)
    r2 <- data.frame(start=0,
                     end=1,
                     auclast=TRUE,
                     stringsAsFactors=FALSE)
    r2[,setdiff(nameorder, names(r2))] <- FALSE
    expect_equal(check.interval.specification(d2),
                 r2[,nameorder],
                 info="Expand a data frame interval specification with only one request given")

    ## start and end must both be specified
    d3 <- data.frame(start=0)
    expect_error(check.interval.specification(d3),
                 regexp="Column\\(s\\) 'end' missing from interval specification",
                 info="Confirm end column is in interval specification")
    d4 <- data.frame(end=1)
    expect_error(check.interval.specification(d4),
                 regexp="Column\\(s\\) 'start' missing from interval specification",
                 info="Confirm start column is in interval specification")
    d5 <- data.frame(blah=5)
    expect_error(check.interval.specification(d5),
                 regexp="Column\\(s\\) 'start', 'end' missing from interval specification",
                 info="Confirm start and end columns are in interval specification")

    ## Ensure that there are data
    d6 <- data.frame()
    expect_error(check.interval.specification(d6),
                 regexp="interval specification has no rows",
                 info="It is an error to have an interval specification with no rows")

    ## Confirm specific column values required
    d7 <- data.frame(start=as.numeric(NA), end=1)
    expect_error(
      check.interval.specification(d7),
      regexp="Interval specification may not have NA for the starting time",
      info="Interval specification may not have NA for the starting time")
    d8 <- data.frame(start=0, end=as.numeric(NA))
    expect_error(
      check.interval.specification(d8),
      regexp="Interval specification may not have NA for the end time",
      info="Interval specification may not have NA for the end time")
    d9 <- data.frame(start=1, end=1)
    expect_error(check.interval.specification(d9),
                 regexp="start must be < end",
                 info="In interval specification, start must be < end (they are equal).")
    d10 <- data.frame(start=1, end=0)
    expect_error(check.interval.specification(d10),
                 regexp="start must be < end",
                 info="In interval specification, start must be < end (end is less).")
    d11 <- data.frame(start=Inf, end=1)
    expect_error(check.interval.specification(d11),
                 regexp="start may not be infinite",
                 info="In interval specification, start may not be infinite (positive infinity).")
    d12 <- data.frame(start=-Inf, end=1)
    expect_error(check.interval.specification(d12),
                 regexp="start may not be infinite",
                 info="In interval specification, start may not be infinite (negative infinity).")

    ## But it is OK to have an infinite end
    d13 <- data.frame(start=0, end=Inf)
    r13 <- data.frame(start=0,
                      end=Inf,
                      stringsAsFactors=FALSE)
    r13[,setdiff(nameorder, names(r13))] <- FALSE
    expect_warning(d13.check <- check.interval.specification(d13))
    expect_equal(d13.check, r13[,nameorder],
                 info="In interval specification, end may be infinite (positive infinity).")
    expect_error(check.interval.specification(data.frame(start=0, end=-Inf)),
                 info="In interval specification, end may not be negative infinity (start is 0).")
    expect_error(check.interval.specification(data.frame(start=-Inf, end=-Inf)),
                 info="In interval specification, end may not be negative infinity (start is -Inf).")
    
    ## When the no-calculation interval specification is not the first,
    ## ensure that is warned correctly
    d14 <- data.frame(start=0, end=24, auclast=c(rep(FALSE, 3), TRUE))
    expect_warning(check.interval.specification(d14),
                   regexp="Nothing to be calculated in interval specification number\\(s\\): 1, 2, 3",
                   info="Warn when nothing is to be calculated in all rows of the specification.")

    d14 <- data.frame(start=0, end=24, auclast=c(rep(TRUE, 3), FALSE))
    expect_warning(check.interval.specification(d14),
                   regexp="Nothing to be calculated in interval specification number\\(s\\): 4",
                   info="Warn when nothing is to be calculated in one but not all rows of the specification.")

    ## Other information is passed through untouched after all the
    ## calculation columns
    d15 <- data.frame(start=0, end=Inf, treatment="foo",
                      stringsAsFactors=FALSE)
    r15 <- data.frame(start=0,
                      end=Inf,
                      treatment="foo",
                      stringsAsFactors=FALSE)
    r15[,setdiff(nameorder, names(r15))] <- FALSE
    expect_warning(v15 <- check.interval.specification(d15))
    expect_equal(v15, r15[,c(nameorder, "treatment")],
                 info="Extra information is maintained in the interval specification.")

    d16 <- data.frame(start=factor(0), end=1)
    expect_error(check.interval.specification(d16),
                 regexp="Interval column 'start' should not be a factor",
                 info="Start must be numeric and not a factor.")
    d17 <- data.frame(start=0, end=factor(1))
    expect_error(check.interval.specification(d17),
                 regexp="Interval column 'end' should not be a factor",
                 info="End must be numeric and not a factor.")
})

test_that("check.interval.deps", {
  ## Get the current name order of the expected results
  nameorder <- names(check.interval.deps(data.frame(start=0, end=1, cmax=TRUE)))

  r1 <- data.frame(start=0,
                   end=24,
                   lambda.z=TRUE,
                   clast.obs=TRUE,
                   aucinf.obs=TRUE)
  r1[,setdiff(nameorder, names(r1))] <- FALSE
  expect_equal(check.interval.deps(data.frame(start=0, end=24, aucinf.obs=TRUE)),
               r1[,nameorder],
               info="Confirm that the interval dependencies are accurately added")
})

test_that("get.parameter.deps", {
  expect_error(
    get.parameter.deps("foo"),
    regexp="`x` must be the name of an NCA parameter listed by the function `get.interval.cols()`",
    info="The argument must be a parameter.",
    fixed=TRUE
  )
  expect_equal(
    get.parameter.deps("kel.obs"),
    c("kel.obs"),
    info="Parameters that have nothing that depend on them return themselves only."
  )
  expect_equal(
    get.parameter.deps("ctrough"),
    c("ctrough", "ctrough.dn", "ptr"),
    info="Parameters with formalsmap-related dependencies return themselves and the formalsmap-related dependencies."
  )
  expect_equal(
    get.parameter.deps("start"),
    character(0),
    info="Special columns that are not actually parameters have no dependencies (including themselves)."
  )
  expect_equal(
    get.parameter.deps("cl.obs"),
    c("cl.obs", "vss.iv.obs", "vss.obs", "vz.obs"),
    info="Parameters with dependencies return them."
  )
})

test_that("check.intervals requires a valid value", {
  expect_error(
    check.interval.specification(data.frame(start=0, end=1, cmax="A")),
    regexp="Invalid value(s) in column cmax:A",
    fixed=TRUE
  )
})

test_that("check.intervals works with tibble input (fix #141)", {
  e.dat <-
    data.frame(
      conc=c(0.5,2,5,9.2,12,2,1.85,1.08,0.5,0.3,2.4,4.5,10.2,15,2.6,1.65,1.1,
             0.5,2,5,9.2,12,2,1.85,1.08,NA,0.3,2.4,4.5,10.2,15,2.6,1.65,1.1),
      time=c(seq(264.2,312.2,3),seq(264,312,3)),
      ARM=rep(c(rep(1,8),rep(2,9)),2),
      SUBJ=c(rep(1,17),rep(2,17)),
      Dose=c(rep(5,17)),rep(5,17)
    )
  
  intervals_manual_first <-
    e.dat %>%
    dplyr::group_by(SUBJ) %>%
    dplyr::summarize(
      start=time[dplyr::between(time, 264, 265)],
      end=time[dplyr::between(time, 288, 289)]
    )
  intervals_manual_second <-
    e.dat %>%
    dplyr::group_by(SUBJ) %>%
    dplyr::summarize(
      start=time[dplyr::between(time, 288, 289)],
      end=time[dplyr::between(time, 312, 313)]
    )
  intervals_manual <-
    dplyr::bind_rows(
      intervals_manual_first,
      intervals_manual_second
    ) %>%
    dplyr::mutate(
      auclast=TRUE,
      aucall=TRUE,
      tlast=TRUE
    )
  # There is some other issue here where intervals are having an issue being a tibble
  expect_equal(
    check.interval.specification(intervals_manual)$start,
    intervals_manual$start
  )
})
