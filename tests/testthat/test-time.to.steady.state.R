context("Time to steady-state")

test_that("pk.tss.data.prep", {
  conc.test <- 1:5
  time.test <- 0:4
  subject.test <- letters[1:5]
  treatment.test <- LETTERS[1:5]
  time.dosing.test <- 0
  ## Confirm that any NAs in time.dosing are an error
  expect_error(pk.tss.data.prep(conc=conc.test,
                                time=time.test,
                                subject=subject.test,
                                treatment=treatment.test,
                                time.dosing=NA),
               regexp="time.dosing may not contain any NA values")

  expect_error(pk.tss.data.prep(conc=conc.test,
                                time=time.test,
                                subject=subject.test,
                                treatment=treatment.test,
                                time.dosing=c(0, NA)),
               regexp="time.dosing may not contain any NA values")

  ## Confirm that conc and subject must be the same length
  expect_error(pk.tss.data.prep(conc=conc.test,
                                time=time.test,
                                subject=subject.test[-1],
                                treatment=treatment.test,
                                time.dosing=time.dosing.test),
               regexp="arguments imply differing number of rows: 5, 4")

  ## Confirm that conc and treatment must be the same length
  expect_error(pk.tss.data.prep(conc=conc.test,
                                time=time.test,
                                subject=subject.test,
                                treatment=treatment.test[-1],
                                time.dosing=time.dosing.test),
               regexp="arguments imply differing number of rows: 5, 4")

  ## If removed down to one treatment, treatment is not a column of
  ## the output.
  expect_equal(pk.tss.data.prep(conc=conc.test,
                                time=time.test,
                                subject=subject.test,
                                treatment=treatment.test,
                                time.dosing=time.dosing.test),
               data.frame(conc=1, time=0))

  ## If no treatment is given, it still works.
  expect_equal(pk.tss.data.prep(conc=conc.test,
                                time=time.test,
                                subject=subject.test,
                                time.dosing=time.dosing.test),
               data.frame(conc=1, time=0))

  ## If no subject is given, it still works
  expect_equal(pk.tss.data.prep(conc=conc.test,
                                time=time.test,
                                treatment=treatment.test,
                                time.dosing=time.dosing.test),
               data.frame(conc=1, time=0))
  
  ## What do we actually expect to get out?
  ## Check a single row output dropping treatment
  expect_equal(pk.tss.data.prep(conc=conc.test,
                                time=time.test,
                                subject=subject.test,
                                treatment=treatment.test,
                                time.dosing=time.dosing.test,
                                subject.dosing=subject.test),
               data.frame(time=0, conc=1))

  ## Check a single row output with no treatment given
  expect_equal(pk.tss.data.prep(conc=conc.test,
                                time=time.test,
                                subject=subject.test,
                                time.dosing=time.dosing.test),
               data.frame(conc=1, time=0))

  ## Check a multi-row output with treatments kept
  conc.test <- 1:10
  time.test <- rep(0:4, 2)
  subject.test <- letters[1:10]
  treatment.test <- LETTERS[1:10]
  time.dosing.test <- 0
  expect_equal(pk.tss.data.prep(conc=conc.test,
                                time=time.test,
                                subject=subject.test,
                                treatment=treatment.test,
                                time.dosing=time.dosing.test),
               data.frame(conc=c(1, 6), time=c(0, 0),
                          subject=c("a", "f"),
                          treatment=c("A", "F")),
               check.attributes=FALSE)

  ## Check a multi-row output with treatments dropped
  conc.test <- 1:10
  time.test <- rep(0:4, 2)
  subject.test <- letters[1:10]
  treatment.test <- rep(LETTERS[1:5], 2)
  time.dosing.test <- 0
  expect_equal(pk.tss.data.prep(conc=conc.test,
                                time=time.test,
                                subject=subject.test,
                                treatment=treatment.test,
                                time.dosing=time.dosing.test),
               data.frame(conc=c(1, 6), time=c(0, 0),
                          subject=c("a", "f")),
               check.attributes=FALSE)
})

## This data will be used multiple times in testing, and it is
## nontrivial to create.
generate.data <- function() {
  set.seed(5)
  tmpdata <-
    merge(data.frame(subject=factor(1:10),
                     css.re=rnorm(10, sd=0.2),
                     tss.re=rnorm(10, sd=0.2),
                     treatment=rep(c("A", "B"), each=5)),
          data.frame(treatment=c("A", "B"),
                     css.mean=c(5, 10),
                     tss.mean=5))
  tmpdata <- merge(tmpdata, data.frame(time=0:14))
  tmpdata$conc.resid <- rnorm(nrow(tmpdata), sd=0.05)
  tmpdata$conc <- with(tmpdata,
                       css.mean*exp(css.re+conc.resid)*
                       (1-exp(log(1-0.9)*time/(tss.mean*exp(tss.re)))))
  tmpdata
}
## Note that this graphically represents the test
## library(latticeExtra)
## (xyplot(conc~time|treatment,
##         groups=subject,
##         data=generate.data(),
##         type="l") +
##  layer(panel.abline(h=5, col="gray", lty=2), packets=1) +
##  layer(panel.abline(h=10, col="gray", lty=2), packets=2))

test_that("pk.tss.stepwise.linear", {
  tmpdata <- generate.data()
  expect_equal(
    pk.tss.stepwise.linear(conc=tmpdata$conc,
                           time=tmpdata$time,
                           subject=tmpdata$subject,
                           treatment=tmpdata$treatment,
                           time.dosing=0:14,
                           verbose=FALSE),
    data.frame(tss.stepwise.linear=7),
    info="pk.tss.stepwise.linear 1")

  expect_warning(
    pk.tss.stepwise.linear(conc=tmpdata$conc,
                           time=tmpdata$time,
                           subject=tmpdata$subject,
                           treatment=tmpdata$treatment,
                           min.points=c(3, 4),
                           time.dosing=0:14,
                           verbose=FALSE),
    regex="Only first value of min.points is used",
    info="pk.tss.stepwise.linear 2")

  ## Check the level input
  expect_error(
    pk.tss.stepwise.linear(conc=tmpdata$conc,
                           time=tmpdata$time,
                           subject=tmpdata$subject,
                           treatment=tmpdata$treatment,
                           time.dosing=0:14,
                           level="A",
                           verbose=FALSE),
    regex="level must be a number",
    info="pk.tss.stepwise.linear 3")

  expect_error(
    pk.tss.stepwise.linear(conc=tmpdata$conc,
                           time=tmpdata$time,
                           subject=tmpdata$subject,
                           treatment=tmpdata$treatment,
                           time.dosing=0:14,
                           level=2,
                           verbose=FALSE),
    regex="level must be between 0 and 1, exclusive",
    info="pk.tss.stepwise.linear 4")
  
  expect_error(
    pk.tss.stepwise.linear(conc=tmpdata$conc,
                           time=tmpdata$time,
                           subject=tmpdata$subject,
                           treatment=tmpdata$treatment,
                           time.dosing=0:14,
                           level=-1,
                           verbose=FALSE),
    regex="level must be between 0 and 1, exclusive",
    info="pk.tss.stepwise.linear 5")
  expect_error(
    pk.tss.stepwise.linear(conc=tmpdata$conc,
                           time=tmpdata$time,
                           subject=tmpdata$subject,
                           treatment=tmpdata$treatment,
                           time.dosing=0:14,
                           level=0,
                           verbose=FALSE),
    regex="level must be between 0 and 1, exclusive",
    info="pk.tss.stepwise.linear 6")
  expect_error(
    pk.tss.stepwise.linear(conc=tmpdata$conc,
                           time=tmpdata$time,
                           subject=tmpdata$subject,
                           treatment=tmpdata$treatment,
                           time.dosing=0:14,
                           level=1,
                           verbose=FALSE),
    regex="level must be between 0 and 1, exclusive",
    info="pk.tss.stepwise.linear 7")

  expect_warning(
    pk.tss.stepwise.linear(conc=tmpdata$conc,
                           time=tmpdata$time,
                           subject=tmpdata$subject,
                           treatment=tmpdata$treatment,
                           time.dosing=0:14,
                           level=c(0.95, 0.99),
                           verbose=FALSE),
    regex="Only first value of level is being used",
    info="pk.tss.stepwise.linear 8")

  ## This is mainly to test verbosity
  expect_warning(
    pk.tss.stepwise.linear(conc=tmpdata$conc,
                           time=tmpdata$time,
                           subject=tmpdata$subject,
                           treatment=tmpdata$treatment,
                           time.dosing=0:14,
                           level=c(0.95, 0.99),
                           verbose=FALSE),
    regex="Only first value of level is being used",
    info="pk.tss.stepwise.linear 9")

  ## Ensure that the first value really is used
  expect_warning(v1 <-
    pk.tss.stepwise.linear(conc=tmpdata$conc,
                           time=tmpdata$time,
                           subject=tmpdata$subject,
                           treatment=tmpdata$treatment,
                           time.dosing=0:14,
                           level=c(0.8, 0.99),
                           verbose=FALSE))
  expect_equal(
    v1,
    pk.tss.stepwise.linear(conc=tmpdata$conc,
                           time=tmpdata$time,
                           subject=tmpdata$subject,
                           treatment=tmpdata$treatment,
                           time.dosing=0:14,
                           level=0.8,
                           verbose=FALSE),
    info="pk.tss.stepwise.linear 10")

  ## Confirm testing for minimum number of data points
  expect_warning(
    v1 <- pk.tss.stepwise.linear(conc=tmpdata$conc,
                                 time=tmpdata$time,
                                 subject=tmpdata$subject,
                                 treatment=tmpdata$treatment,
                                 time.dosing=0:1,
                                 min.points=3,
                                 level=0.99,
                                 verbose=FALSE),
    regexp="After removing non-dosing time points, insufficient data remains for tss calculation")
  expect_equal(v1, NA)
  expect_warning(
    pk.tss.stepwise.linear(conc=tmpdata$conc,
                           time=tmpdata$time,
                           subject=tmpdata$subject,
                           treatment=tmpdata$treatment,
                           time.dosing=0:1,
                           min.points=3,
                           level=0.99,
                           verbose=FALSE),
    regexp="After removing non-dosing time points, insufficient data remains for tss calculation")

  ## Confirm the glm model works when subject-level data is not given.
  expect_equal(
    pk.tss.stepwise.linear(conc=tmpdata$conc,
                           time=tmpdata$time,
                           treatment=tmpdata$treatment,
                           time.dosing=0:14,
                           verbose=FALSE),
    data.frame(tss.stepwise.linear=5),
    info="pk.tss.stepwise.linear no subject")
})

test_that("pk.tss.monoexponential", {
  tmpdata <- generate.data()
  expect_equal(
    pk.tss.monoexponential(conc=tmpdata$conc,
                           time=tmpdata$time,
                           subject=tmpdata$subject,
                           treatment=tmpdata$treatment,
                           time.dosing=0:14,
                           verbose=FALSE),
    data.frame(subject=factor(c(1, 10, 2:9)),
               tss.monoexponential.population=4.57618156812974,
               tss.monoexponential.popind=c(
                 5.14156352865421, 4.64862524830397, 4.45956707917941,
                 4.41492203844343, 4.6782583033301, 4.0823047621517,
                 4.96242115751172, 4.52424147509819, 3.70338406668837,
                 5.1465280219363),
               treatment=c("A", "B", "A", "A", "A",
                 "A", "B", "B", "B", "B"),
               tss.monoexponential.individual=c(
                 5.87784329336254, 4.71066285661623, 4.51882509145954,
                 3.91269286106442, 4.74475071729459, 3.99341726779716,
                 5.08737230904342, 4.50068650719192, 3.4876172020751,
                 5.35051537086801),
               tss.monoexponential.single=c(
                 4.49364716304018, 4.5987103152447, 4.49364716304018,
                 4.49364716304018, 4.49364716304018, 4.49364716304018,
                 4.5987103152447, 4.5987103152447, 4.5987103152447,
                 4.5987103152447)),
    tolerance=1e-4,
    check.attributes=FALSE,
    info="pk.tss.monoexponential 1")

  ## Warnings and errors
  expect_warning(
    pk.tss.monoexponential(conc=tmpdata$conc,
                           time=tmpdata$time,
                           subject=tmpdata$subject,
                           treatment=tmpdata$treatment,
                           time.dosing=0:14,
                           tss.fraction=c(0.5, 0.8)),
    regexp="Only first value of tss.fraction is being used")
  expect_error(
    pk.tss.monoexponential(conc=tmpdata$conc,
                           time=tmpdata$time,
                           subject=tmpdata$subject,
                           treatment=tmpdata$treatment,
                           time.dosing=0:14,
                           tss.fraction=0),
    regexp="tss.fraction must be between 0 and 1, exclusive")
  expect_error(
    pk.tss.monoexponential(conc=tmpdata$conc,
                           time=tmpdata$time,
                           subject=tmpdata$subject,
                           treatment=tmpdata$treatment,
                           time.dosing=0:14,
                           tss.fraction=1),
    regexp="tss.fraction must be between 0 and 1, exclusive")
  expect_error(
    pk.tss.monoexponential(conc=tmpdata$conc,
                           time=tmpdata$time,
                           subject=tmpdata$subject,
                           treatment=tmpdata$treatment,
                           time.dosing=0:14,
                           tss.fraction=-1),
    regexp="tss.fraction must be between 0 and 1, exclusive")
  expect_error(
    pk.tss.monoexponential(conc=tmpdata$conc,
                           time=tmpdata$time,
                           subject=tmpdata$subject,
                           treatment=tmpdata$treatment,
                           time.dosing=0:14,
                           tss.fraction=2),
    regexp="tss.fraction must be between 0 and 1, exclusive")
  expect_warning(
    pk.tss.monoexponential(conc=tmpdata$conc,
                           time=tmpdata$time,
                           subject=tmpdata$subject,
                           treatment=tmpdata$treatment,
                           time.dosing=0:14,
                           tss.fraction=0.5),
    regexp="tss.fraction is usually >= 0.8")
})

test_that("pk.tss", {
  ## Ensure that pk.tss will go to the correct type of model
  tmpdata <- generate.data()
  expect_equal(
    pk.tss(conc=tmpdata$conc,
           time=tmpdata$time,
           subject=tmpdata$subject,
           treatment=tmpdata$treatment,
           time.dosing=0:14,
           verbose=FALSE,
           type="monoexponential"),
    pk.tss.monoexponential(conc=tmpdata$conc,
                           time=tmpdata$time,
                           subject=tmpdata$subject,
                           treatment=tmpdata$treatment,
                           time.dosing=0:14,
                           verbose=FALSE))

  expect_equal(
    pk.tss(conc=tmpdata$conc,
           time=tmpdata$time,
           subject=tmpdata$subject,
           treatment=tmpdata$treatment,
           time.dosing=0:14,
           verbose=FALSE,
           type="stepwise.linear"),
    pk.tss.stepwise.linear(conc=tmpdata$conc,
                           time=tmpdata$time,
                           subject=tmpdata$subject,
                           treatment=tmpdata$treatment,
                           time.dosing=0:14,
                           verbose=FALSE))

  ## pk.tss will calculate both if requested
  expect_equal(
    pk.tss(conc=tmpdata$conc,
           time=tmpdata$time,
           subject=tmpdata$subject,
           treatment=tmpdata$treatment,
           time.dosing=0:14,
           verbose=FALSE,
           type=c("monoexponential", "stepwise.linear")),
    merge(pk.tss.monoexponential(conc=tmpdata$conc,
                           time=tmpdata$time,
                           subject=tmpdata$subject,
                           treatment=tmpdata$treatment,
                           time.dosing=0:14,
                           verbose=FALSE),
          pk.tss.stepwise.linear(conc=tmpdata$conc,
                                 time=tmpdata$time,
                                 subject=tmpdata$subject,
                                 treatment=tmpdata$treatment,
                                 time.dosing=0:14,
                                 verbose=FALSE),
          all=TRUE))
})
