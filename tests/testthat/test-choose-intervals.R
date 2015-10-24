context("Determining tau for AUCs")

test_that("find.tau", {
  ## Regularly spaced intervals give the regular spacing
  expect_equal(
    find.tau(sort(unique(c(seq(0, 168, 12),
                           seq(4, 168, 12)))),
             tau.choices=NA),
    12)
  expect_equal(
    find.tau(sort(unique(c(seq(0, 168, 12),
                           seq(4, 168+12, 12)))),
             tau.choices=NA),
    12)
  expect_equal(find.tau(0:10, tau.choices=NA), 1)
  ## It overrides tau.choices if everything is equally spaced.
  expect_equal(find.tau(0:10, tau.choices=c(24, 168)), 1)
  expect_equal(find.tau(seq(0, 100, by=10),
                        tau.choices=c(24, 168)), 10)
  expect_equal(find.tau(seq(0, 48, by=24),
                        tau.choices=c(24, 168)), 24)
  ## Alternatively spaced intervals give the alternative spacing
  expect_equal(find.tau(c(seq(0, 48, by=24),
                          seq(10, 48, by=24)),
                        tau.choices=c(24, 168)),
               24)
  ## Smaller interval spacing adheres to tau.choices with unequal
  ## spacing.
  expect_equal(find.tau(c(seq(0, 48, by=12),
                          seq(10, 48, by=12)),
                        tau.choices=c(24, 168)),
               24)
  ## It works with more complex spacing and many intervals
  expect_equal(find.tau(
    sort(unique(c(seq(0, 168, by=24),
                  seq(4, 168, by=24),
                  seq(10, 168, by=24)))),
    tau.choices=c(24, 168)),
    24)
  expect_equal(find.tau(c(0, 5, 19, 30) + rep((0:5)*168, each=4),
                        tau.choices=c(24, 168)),
               168)
  ## If there is only one dosing time, return 0-- regardless of the
  ## option.
  expect_equal(find.tau(rep(5, 5), tau.choices=c(24, 168)), 0)
  expect_equal(find.tau(rep(0, 5), tau.choices=c(24, 168)), 0)
  ## If everything is NA, return NA
  expect_equal(find.tau(NA,
                        tau.choices=c(24, 168)),
               NA)
  expect_equal(find.tau(rep(NA, 10),
                        tau.choices=c(24, 168)),
               NA)
  ## If there is no sequence, return NA
  expect_equal(find.tau(c(0, 1, 3, 5, 9),
                        tau.choices=c(24, 168)),
               NA)
  expect_equal(find.tau(c(0, 1, 3, 5, 9, 24),
                        tau.choices=NA),
               NA)
})

test_that("choose.auc.intervals", {
  tmp.single.dose.auc <-
    check.interval.specification(
      data.frame(start=0,
                 end=c(24, Inf),
                 auclast=c(TRUE, FALSE),
                 aucinf=c(FALSE, TRUE),
                 half.life=c(FALSE, TRUE),
                 stringsAsFactors=FALSE))

  ## Check the inputs
  expect_error(choose.auc.intervals(NA, 1, tmp.single.dose.auc),
               regexp="time.conc may not have any NA values")
  expect_error(choose.auc.intervals(1, NA, tmp.single.dose.auc),
               regexp="time.dosing may not have any NA values")
  ## The below test is a bit of a non-sequeter-- essentially, it just
  ## needs to return a 0-row data frame.
  expect_error(choose.auc.intervals(1, 1,
                                    single.dose.aucs=data.frame()),
               regexp="interval specification has no rows")
  ## It adjusts single dose AUCs by the starting time
  expect_equal(choose.auc.intervals(1, 1,
                                    single.dose.aucs=tmp.single.dose.auc),
               check.interval.specification(
                 data.frame(start=1,
                            end=c(25, Inf),
                            auclast=c(TRUE, FALSE),
                            aucinf=c(FALSE, TRUE),
                            half.life=c(FALSE, TRUE))))

  ## Find intervals for two doses with PK at both points and one in
  ## between.
  expect_equal(choose.auc.intervals(c(1, 2, 3), c(1, 3),
                                    single.dose.aucs=tmp.single.dose.auc),
               check.interval.specification(
                 data.frame(start=1,
                            end=3,
                            cmax=TRUE,
                            tmax=TRUE,
                            auclast=TRUE)))
  ## Find intervals for two doses with PK at both points, one in
  ## between, and one after asking for AUClast after the second dose
  ## but no half-life.
  expect_equal(choose.auc.intervals(1:5, c(1, 3),
                                    single.dose.aucs=tmp.single.dose.auc),
               check.interval.specification(
                 data.frame(start=c(1, 3),
                            end=c(3, 5),
                            cmax=TRUE,
                            tmax=TRUE,
                            auclast=TRUE)))
  ## Find intervals for two doses with PK at both points, one in
  ## between, and one after asking for AUClast after the second dose
  ## with half-life.
  expect_equal(choose.auc.intervals(1:6, c(1, 3),
                                    single.dose.aucs=tmp.single.dose.auc),
               check.interval.specification(
                 data.frame(start=c(1, 3, 3),
                            end=c(3, 5, Inf),
                            auclast=c(TRUE, TRUE, FALSE),
                            cmax=c(TRUE, TRUE, FALSE),
                            tmax=c(TRUE, TRUE, FALSE),
                            half.life=c(FALSE, FALSE, TRUE))))
  ## Some doses have PK betwen them, some not.
  expect_equal(choose.auc.intervals(1:6, c(1, 3, 5, 7, 9),
                                    single.dose.aucs=tmp.single.dose.auc),
               check.interval.specification(
                 data.frame(start=c(1, 3),
                            end=c(3, 5),
                            auclast=TRUE,
                            cmax=TRUE,
                            tmax=TRUE)))
  ## Find intervals when some doses do not have AUCs between them
  ## (pairs of doses with trough but no PK between)  
  expect_equal(choose.auc.intervals(c(1, 2, 3, 5, 6, 7),
                                    c(1, 3, 5, 7, 9),
                                    single.dose.aucs=tmp.single.dose.auc),
               check.interval.specification(
                 data.frame(start=c(1, 5),
                            end=c(3, 7),
                            cmax=TRUE,
                            tmax=TRUE,
                            auclast=TRUE)))
  ## Find intervals for two doses with PK at both points, one in
  ## between, and one after asking for AUClast after the second dose
  ## with half-life.  Since tau is not detectable, no half-life at the
  ## end.
  expect_equal(choose.auc.intervals(1:6, c(1, 3, 5, 9),
                                    single.dose.aucs=tmp.single.dose.auc),
               check.interval.specification(
                 data.frame(start=c(1, 3),
                            end=c(3, 5),
                            auclast=TRUE,
                            cmax=TRUE,
                            tmax=TRUE)))
})
