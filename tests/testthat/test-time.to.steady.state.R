test_that("pk.tss.data.prep", {
  conc.test <- 1:5
  time.test <- 0:4
  subject.test <- letters[1:5]
  treatment.test <- LETTERS[1:5]
  time.dosing.test <- 0
  # Confirm that any NAs in time.dosing are an error
  expect_error(
    pk.tss.data.prep(
      conc=conc.test,
      time=time.test,
      subject=subject.test,
      treatment=treatment.test,
      time.dosing=NA
    ),
    regexp="time.dosing may not contain any NA values"
  )
  expect_error(
    pk.tss.data.prep(
      conc=conc.test,
      time=time.test,
      subject.dosing=subject.test,
      treatment=treatment.test,
      time.dosing=NA
    ),
    regexp="Cannot give subject.dosing without subject"
  )
  expect_equal(
    pk.tss.data.prep(
      conc=conc.test,
      time=time.test,
      time.dosing=0
    ),
    data.frame(conc=1, time=0)
  )
  
  expect_error(pk.tss.data.prep(conc=conc.test,
                                time=time.test,
                                subject=subject.test,
                                treatment=treatment.test,
                                time.dosing=c(0, NA)),
               regexp="time.dosing may not contain any NA values")

  # Confirm that conc and subject must be the same length
  expect_error(pk.tss.data.prep(conc=conc.test,
                                time=time.test,
                                subject=subject.test[-1],
                                treatment=treatment.test,
                                time.dosing=time.dosing.test),
               regexp="arguments imply differing number of rows: 5, 4")

  # Confirm that conc and treatment must be the same length
  expect_error(pk.tss.data.prep(conc=conc.test,
                                time=time.test,
                                subject=subject.test,
                                treatment=treatment.test[-1],
                                time.dosing=time.dosing.test),
               regexp="arguments imply differing number of rows: 5, 4")

  # If removed down to one treatment, treatment is not a column of
  # the output.
  expect_equal(pk.tss.data.prep(conc=conc.test,
                                time=time.test,
                                subject=subject.test,
                                treatment=treatment.test,
                                time.dosing=time.dosing.test),
               data.frame(conc=1, time=0))

  # If no treatment is given, it still works.
  expect_equal(pk.tss.data.prep(conc=conc.test,
                                time=time.test,
                                subject=subject.test,
                                time.dosing=time.dosing.test),
               data.frame(conc=1, time=0))

  # If no subject is given, it still works
  expect_equal(pk.tss.data.prep(conc=conc.test,
                                time=time.test,
                                treatment=treatment.test,
                                time.dosing=time.dosing.test),
               data.frame(conc=1, time=0))
  
  # What do we actually expect to get out?
  # Check a single row output dropping treatment
  expect_equal(pk.tss.data.prep(conc=conc.test,
                                time=time.test,
                                subject=subject.test,
                                treatment=treatment.test,
                                time.dosing=time.dosing.test,
                                subject.dosing=subject.test),
               data.frame(time=0, conc=1))

  # Check a single row output with no treatment given
  expect_equal(pk.tss.data.prep(conc=conc.test,
                                time=time.test,
                                subject=subject.test,
                                time.dosing=time.dosing.test),
               data.frame(conc=1, time=0))

  # Check a multi-row output with treatments kept
  conc.test <- 1:10
  time.test <- rep(0:4, 2)
  subject.test <- letters[1:10]
  treatment.test <- LETTERS[1:10]
  time.dosing.test <- 0
  expect_equal(
    pk.tss.data.prep(
      conc=conc.test,
      time=time.test,
      subject=subject.test,
      treatment=treatment.test,
      time.dosing=time.dosing.test
    ),
    data.frame(
      conc=c(1, 6), time=c(0, 0),
      subject=factor(c("a", "f")),
      treatment=factor(c("A", "F")),
      stringsAsFactors=FALSE
    ),
    ignore_attr = TRUE
  )

  # Check a multi-row output with treatments dropped
  conc.test <- 1:10
  time.test <- rep(0:4, 2)
  subject.test <- letters[1:10]
  treatment.test <- rep(LETTERS[1:5], 2)
  time.dosing.test <- 0
  expect_equal(
    pk.tss.data.prep(
      conc=conc.test,
      time=time.test,
      subject=subject.test,
      treatment=treatment.test,
      time.dosing=time.dosing.test
    ),
    data.frame(
      conc=c(1, 6),
      time=c(0, 0),
      subject=factor(c("a", "f")),
      stringsAsFactors=FALSE
    ),
    ignore_attr = TRUE
  )
})

# This data will be used multiple times in testing, and it is
# nontrivial to create.
generate.data <- function() {
  set.seed(5)
  tmpdata <-
    merge(
      data.frame(
        subject=factor(1:10),
        css.re=rnorm(10, sd=0.2),
        tss.re=rnorm(10, sd=0.2),
        treatment=rep(c("A", "B"), each=5),
        stringsAsFactors=FALSE
      ),
      data.frame(
        treatment=c("A", "B"),
        css.mean=c(5, 10),
        tss.mean=5,
        stringsAsFactors=FALSE
      )
    )
  tmpdata <- merge(tmpdata, data.frame(time=0:14))
  tmpdata$conc.resid <- rnorm(nrow(tmpdata), sd=0.05)
  tmpdata$conc <- with(tmpdata,
                       css.mean*exp(css.re+conc.resid)*
                       (1-exp(log(1-0.9)*time/(tss.mean*exp(tss.re)))))
  tmpdata
}

# Note that this graphically represents the test
# library(latticeExtra)
# (xyplot(conc~time|treatment,
#         groups=subject,
#         data=generate.data(),
#         type="l") +
#  layer(panel.abline(h=5, col="gray", lty=2), packets=1) +
#  layer(panel.abline(h=10, col="gray", lty=2), packets=2))

test_that("pk.tss.stepwise.linear", {
  tmpdata <- generate.data()
  expect_equal(
    pk.tss.stepwise.linear(conc=tmpdata$conc,
                           time=tmpdata$time,
                           subject=tmpdata$subject,
                           treatment=tmpdata$treatment,
                           time.dosing=0:14,
                           verbose=FALSE),
    data.frame(tss.stepwise.linear=7)
  )

  expect_warning(
    pk.tss.stepwise.linear(conc=tmpdata$conc,
                           time=tmpdata$time,
                           subject=tmpdata$subject,
                           treatment=tmpdata$treatment,
                           min.points=c(3, 4),
                           time.dosing=0:14,
                           verbose=FALSE),
    regexp="Only first value of min.points is used"
  )

  expect_error(
    pk.tss.stepwise.linear(
      conc=tmpdata$conc,
      time=tmpdata$time,
      subject=tmpdata$subject,
      treatment=tmpdata$treatment,
      time.dosing=0:14,
      min.points="A",
      level="A",
      verbose=FALSE
    ),
    regexp="min.points must be a number"
  )
  expect_error(
    pk.tss.stepwise.linear(
      conc=tmpdata$conc,
      time=tmpdata$time,
      subject=tmpdata$subject,
      treatment=tmpdata$treatment,
      time.dosing=0:14,
      min.points=1,
      level="A",
      verbose=FALSE
    ),
    regexp="min.points must be at least 3"
  )

  expect_error(
    pk.tss.stepwise.linear(conc=tmpdata$conc,
                           time=tmpdata$time,
                           subject=tmpdata$subject,
                           treatment=tmpdata$treatment,
                           time.dosing=0:14,
                           level="A",
                           verbose=FALSE),
    regexp="level must be a number"
  )

  expect_error(
    pk.tss.stepwise.linear(conc=tmpdata$conc,
                           time=tmpdata$time,
                           subject=tmpdata$subject,
                           treatment=tmpdata$treatment,
                           time.dosing=0:14,
                           level=2,
                           verbose=FALSE),
    regexp="level must be between 0 and 1, exclusive"
  )
  
  expect_error(
    pk.tss.stepwise.linear(conc=tmpdata$conc,
                           time=tmpdata$time,
                           subject=tmpdata$subject,
                           treatment=tmpdata$treatment,
                           time.dosing=0:14,
                           level=-1,
                           verbose=FALSE),
    regexp="level must be between 0 and 1, exclusive"
  )
  expect_error(
    pk.tss.stepwise.linear(conc=tmpdata$conc,
                           time=tmpdata$time,
                           subject=tmpdata$subject,
                           treatment=tmpdata$treatment,
                           time.dosing=0:14,
                           level=0,
                           verbose=FALSE),
    regexp="level must be between 0 and 1, exclusive"
  )
  expect_error(
    pk.tss.stepwise.linear(conc=tmpdata$conc,
                           time=tmpdata$time,
                           subject=tmpdata$subject,
                           treatment=tmpdata$treatment,
                           time.dosing=0:14,
                           level=1,
                           verbose=FALSE),
    regexp="level must be between 0 and 1, exclusive"
  )

  expect_warning(
    pk.tss.stepwise.linear(conc=tmpdata$conc,
                           time=tmpdata$time,
                           subject=tmpdata$subject,
                           treatment=tmpdata$treatment,
                           time.dosing=0:14,
                           level=c(0.95, 0.99),
                           verbose=FALSE),
    regexp="Only first value of level is being used"
  )

  # This is mainly to test verbosity
  expect_warning(
    pk.tss.stepwise.linear(conc=tmpdata$conc,
                           time=tmpdata$time,
                           subject=tmpdata$subject,
                           treatment=tmpdata$treatment,
                           time.dosing=0:14,
                           level=c(0.95, 0.99),
                           verbose=FALSE),
    regexp="Only first value of level is being used"
  )

  # Check outputs
  suppressMessages(
    withr::with_options(
      list(try.outFile=nullfile()),
      expect_message(
        pk.tss.stepwise.linear(
          conc=tmpdata$conc,
          time=tmpdata$time,
          subject=tmpdata$subject,
          treatment=tmpdata$treatment,
          time.dosing=0:14,
          level=0.8,
          verbose=TRUE
        ),
        regexp="Trying 0"
      )
    )
  )
  suppressMessages(
    withr::with_options(
      list(try.outFile=nullfile()),
      expect_message(
        pk.tss.stepwise.linear(
          conc=tmpdata$conc,
          time=tmpdata$time,
          subject=tmpdata$subject,
          treatment=tmpdata$treatment,
          time.dosing=0:14,
          level=0.8,
          verbose=TRUE
        ),
        regexp="Current interval"
      )
    )
  )
  
  # Ensure that the first value really is used
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
    pk.tss.stepwise.linear(
      conc=tmpdata$conc,
      time=tmpdata$time,
      subject=tmpdata$subject,
      treatment=tmpdata$treatment,
      time.dosing=0:14,
      level=0.8,
      verbose=FALSE
    )
  )

  # Confirm testing for minimum number of data points
  expect_warning(
    v1 <- pk.tss.stepwise.linear(
      conc=tmpdata$conc,
      time=tmpdata$time,
      subject=tmpdata$subject,
      treatment=tmpdata$treatment,
      time.dosing=0:1,
      min.points=3,
      level=0.99,
      verbose=FALSE
    ),
    regexp="After removing non-dosing time points, insufficient data remains for tss calculation"
  )
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

  # Confirm the glm model works when subject-level data is not given.
  suppressMessages(
    expect_equal(
      pk.tss.stepwise.linear(conc=tmpdata$conc,
                             time=tmpdata$time,
                             treatment=tmpdata$treatment,
                             time.dosing=0:14,
                             verbose=FALSE),
      data.frame(tss.stepwise.linear=5),
      info="pk.tss.stepwise.linear no subject"
    )
  )
})

test_that("pk.tss.monoexponential", {
  tmpdata <- generate.data()
  expect_warning(
    expect_equal(
      pk.tss.monoexponential(
        conc=tmpdata$conc,
        time=tmpdata$time,
        subject=tmpdata$subject,
        treatment=tmpdata$treatment,
        time.dosing=0:14,
        verbose=FALSE
      ),
      data.frame(
        subject=factor(as.character(c(1, 10, 2:9))),
        tss.monoexponential.population=4.57618156812974,
        tss.monoexponential.popind=c(
          5.14156352865421, 4.64862524830397, 4.45956707917941,
          4.41492203844343, 4.6782583033301, 4.0823047621517,
          4.96242115751172, 4.52424147509819, 3.70338406668837,
          5.1465280219363),
        treatment=
          factor(c("A", "B", "A", "A", "A", "A", "B", "B", "B", "B")),
        tss.monoexponential.individual=c(
          5.87784329336254, 4.71066285661623, 4.51882509145954,
          3.91269286106442, 4.74475071729459, 3.99341726779716,
          5.08737230904342, 4.50068650719192, 3.4876172020751,
          5.35051537086801),
        tss.monoexponential.single=4.56067603534,
        stringsAsFactors=FALSE
      ),
      tolerance=1e-4
    )
  )
  expect_output(
    expect_equal(
      pk.tss.monoexponential(
        conc=c(0, 1000),
        time=0:1,
        subject=c(1, 1),
        treatment=c("A", "A"),
        time.dosing=0:1,
        tss.fraction=0.9,
        output="single"
      ),
      data.frame(tss.monoexponential.single=NA_real_),
      info="Single-subject data fitting works when it does not converge."
    ),
    regexp = "approximate covariance matrix for parameter estimates not of full rank"
  )
})

test_that("pk.tss.monoexponential without treatment", {
  tmpdata <- generate.data()
    expect_equal(
      pk.tss.monoexponential(
        conc=tmpdata$conc,
        time=tmpdata$time,
        subject=tmpdata$subject,
        time.dosing=0:14,
        verbose=FALSE
      ),
      data.frame(
        subject=factor(as.character(c(1, 10, 2:9))),
        tss.monoexponential.population=4.561580,
        tss.monoexponential.popind=c(
          5.028357, 4.672555, 4.536178,
          4.295758, 4.660561, 4.158700,
          4.940674, 4.507900, 3.676531,
          5.138583),
        tss.monoexponential.individual=c(
          5.87784329336254, 4.71066285661623, 4.51882509145954,
          3.91269286106442, 4.74475071729459, 3.99341726779716,
          5.08737230904342, 4.50068650719192, 3.4876172020751,
          5.35051537086801),
        tss.monoexponential.single=4.56067603534,
        stringsAsFactors=FALSE
      ),
      tolerance=1e-4
    )
})

test_that("pk.tss.monoexponential corner case tests", {
  tmpdata <- generate.data()
  # population output, only
  expect_warning(
    expect_equal(
      pk.tss.monoexponential(
        conc=tmpdata$conc,
        time=tmpdata$time,
        subject=tmpdata$subject,
        treatment=tmpdata$treatment,
        time.dosing=0:14,
        output = "population",
        verbose=FALSE
      ),
      data.frame(
        subject=factor(as.character(seq_len(10))),
        tss.monoexponential.population=4.57618156812974,
        stringsAsFactors=FALSE
      ),
      tolerance=1e-4
    )
  )
  # (Pseudo) single treatment, only
  expect_equal(
    pk.tss.monoexponential(
      conc=tmpdata$conc,
      time=tmpdata$time,
      subject=tmpdata$subject,
      time.dosing=0:14,
      output = "population",
      verbose=FALSE
    ),
    data.frame(
      subject=factor(as.character(seq_len(10))),
      tss.monoexponential.population=4.56157960341961,
      stringsAsFactors=FALSE
    ),
    tolerance=1e-4
  )
})

test_that("pk.tss.monoexponential expected warnings and errors", {
  tmpdata <- generate.data()
  expect_error(
    pk.tss.monoexponential(conc=tmpdata$conc,
                           time=tmpdata$time,
                           subject=tmpdata$subject,
                           treatment=tmpdata$treatment,
                           time.dosing=0:14,
                           tss.fraction=factor(1)),
    regexp="tss.fraction must be a number"
  )
  suppressWarnings(
    expect_warning(
      pk.tss.monoexponential(
        conc=tmpdata$conc,
        time=tmpdata$time,
        subject=tmpdata$subject,
        treatment=tmpdata$treatment,
        time.dosing=0:14,
        tss.fraction=c(0.5, 0.8)
      ),
      regexp="Only first value of tss.fraction is being used"
    )
  )
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
  suppressWarnings(
    expect_warning(
      pk.tss.monoexponential(conc=tmpdata$conc,
                             time=tmpdata$time,
                             subject=tmpdata$subject,
                             treatment=tmpdata$treatment,
                             time.dosing=0:14,
                             tss.fraction=0.5),
      regexp="tss.fraction is usually >= 0.8"
    )
  )
  suppressWarnings(
    expect_warning(
      pk.tss.monoexponential(conc=tmpdata$conc,
                             time=tmpdata$time,
                             subject=tmpdata$subject,
                             treatment=tmpdata$treatment,
                             time.dosing=0:14,
                             tss.fraction=c(0.5, 0.8)),
      regexp="Only first value of tss.fraction is being used"
    )
  )
  tmpdata_single_subject <- tmpdata[tmpdata$subject == 1 & tmpdata$treatment == "A", ]
  expect_warning(
    tss_single_subject <-
      pk.tss.monoexponential(
        conc=tmpdata_single_subject$conc,
        time=tmpdata_single_subject$time,
        subject=tmpdata_single_subject$subject,
        treatment=tmpdata_single_subject$treatment,
        time.dosing=0:14,
        tss.fraction=0.8
      ),
    regexp="Cannot give 'population', 'popind', or 'individual' output without multiple subjects of data",
    fixed=TRUE
  )
  expect_equal(
    tss_single_subject,
    data.frame(tss.monoexponential.single=4.108541),
    tolerance=1e-5
  )
})

test_that("pk.tss", {
  # Ensure that pk.tss will go to the correct type of model
  tmpdata <- generate.data()
  suppressWarnings(
    expect_equal(
      pk.tss(
        conc=tmpdata$conc,
        time=tmpdata$time,
        subject=tmpdata$subject,
        treatment=tmpdata$treatment,
        time.dosing=0:14,
        verbose=FALSE,
        type="monoexponential"
      ),
      pk.tss.monoexponential(
        conc=tmpdata$conc,
        time=tmpdata$time,
        subject=tmpdata$subject,
        treatment=tmpdata$treatment,
        time.dosing=0:14,
        verbose=FALSE
      )
    )
  )

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

  # pk.tss will calculate both if requested
  suppressWarnings(
    expect_equal(
      pk.tss(
        conc=tmpdata$conc,
        time=tmpdata$time,
        subject=tmpdata$subject,
        treatment=tmpdata$treatment,
        time.dosing=0:14,
        verbose=FALSE,
        type=c("monoexponential", "stepwise.linear")
      ),
      merge(
        pk.tss.monoexponential(
          conc=tmpdata$conc,
          time=tmpdata$time,
          subject=tmpdata$subject,
          treatment=tmpdata$treatment,
          time.dosing=0:14,
          verbose=FALSE
        ),
        pk.tss.stepwise.linear(
          conc=tmpdata$conc,
          time=tmpdata$time,
          subject=tmpdata$subject,
          treatment=tmpdata$treatment,
          time.dosing=0:14,
          verbose=FALSE
        ),
        all=TRUE
      )
    )
  )
})

test_that("pk.tss.monoexponential with single-subject data", {
  d_prep <- datasets::Theoph[datasets::Theoph$Subject %in% 2, ]
  dose_times <- seq(0, 96-1, by=6)
  d_multidose <-
    superposition(
      conc=d_prep$conc,
      time=d_prep$Time,
      tau=96, # 48 hours
      n.tau=1, # One tau interval (0 to 48 hours)
      dose.times=dose_times
    )
  expect_equal(
    pk.tss.monoexponential(
      conc=d_multidose$conc, time=d_multidose$time, subject=rep(1, nrow(d_multidose)),
      time.dosing=dose_times, subject.dosing=rep(1, length(dose_times)),
      output="single"
    ),
    data.frame(tss.monoexponential.single=22.53),
    tolerance=0.001
  )
})

test_that("verbose being TRUE", {
  d_prep <- datasets::Theoph[datasets::Theoph$Subject %in% 2, ]
  dose_times <- seq(0, 96-1, by=6)
  d_multidose <-
    superposition(
      conc=d_prep$conc,
      time=d_prep$Time,
      tau=96, # 48 hours
      n.tau=1, # One tau interval (0 to 48 hours)
      dose.times=dose_times
    )
  expect_equal(
    pk.tss.monoexponential(
      conc=d_multidose$conc, time=d_multidose$time, subject=rep(1, nrow(d_multidose)),
      time.dosing=dose_times, subject.dosing=rep(1, length(dose_times)),
      output="single", verbose = TRUE
    ),
    data.frame(tss.monoexponential.single=22.53),
    tolerance=0.001
  )
})

test_that("tss Monoexponential no models converged", {
  bad_data <- data.frame(subject=rep(1:2, each = 4), time=0, conc=0, tss.constant=0.9)
  expect_warning(pk.tss.monoexponential.population(data = bad_data, output = "population"),
  regexp = "No population model for monoexponential Tss converged, no results given")
})

test_that("tss Monoexponential no models converged with verbose = TRUE", {
  bad_data <- data.frame(subject=rep(1:2, each = 4), time=0, conc=0, tss.constant=0.9)
  expect_warning(pk.tss.monoexponential.population(data = bad_data, output = "population", verbose = TRUE),
                 regexp = "No population model for monoexponential Tss converged, no results given")
})
