context("PKNCA option setting")

test_that("PKNCA.options", {
  ## Missing/incorrect option names give an error indicating all the
  ## options that are missing.
  expect_error(PKNCA.options("foo"),
               regexp="PKNCA.options does not have value\\(s\\) for foo.")
  expect_error(PKNCA.options("foo", "bar"),
               regexp="PKNCA.options does not have value\\(s\\) for foo, bar.")
  ## A mix of mixxing and extant options only give the missing ones
  expect_error(PKNCA.options("foo", "bar", "tau.choices"),
               regexp="PKNCA.options does not have value\\(s\\) for foo, bar.")
  
  ## Single extant options give their default value
  expect_equal(PKNCA.options("min.hl.points"), 3)

  ## Multiple extant options give a list of their default values
  expect_equal(PKNCA.options("min.hl.points", "min.hl.r.squared"),
               list(min.hl.points=3,
                    min.hl.r.squared=0.9))
  ## The returned options are in the order they were requested
  expect_equal(PKNCA.options("min.hl.points", "min.hl.r.squared"),
               list(min.hl.points=3,
                    min.hl.r.squared=0.9))

  ## Asking for an option using the name argument works the same as
  ## asking for one with a string.
  expect_equal(PKNCA.options(name="single.dose.aucs"),
               PKNCA.options("single.dose.aucs"))
  ## You can request more than one option by using string inputs
  ## combined with a single "name" argument.
  expect_equal(PKNCA.options("first.tmax", name="single.dose.aucs"),
               PKNCA.options("first.tmax", name="single.dose.aucs"))
  ## You cannot both set an option and give it with a name.
  expect_error(PKNCA.options(foo="var", name="foo", value="bar"),
               regexp="Cannot give an option name both with the name argument and as a named argument.")
  ## You cannot both set an option (with a value) and request an
  ## option
  expect_error(PKNCA.options("first.tmax", name="min.span.ratio", value=2),
               regexp="Invalid setting for PKNCA")
  ## You cannot give a value without a name.
  expect_error(PKNCA.options(value=5),
               regexp="Cannot have a value without a name")

  ## Cannot both check and default at the same time
  expect_error(PKNCA.options("adj.r.squared.factor", default=TRUE, check=TRUE),
               regexp="Cannot request both default and check")
  
  ## Confirm that the default state is as expected (setting it first
  ## in case the tests are run in a non-default state)
  PKNCA.options(default=TRUE)
  expect_equal(PKNCA.options(),
               list(adj.r.squared.factor=0.0001,
                    max.missing=0.5,
                    auc.method="lin up/log down",
                    conc.na="drop",
                    conc.blq=list(
                      first="keep",
                      middle="drop",
                      last="keep"),
                    first.tmax=TRUE,
                    allow.tmax.in.half.life=FALSE,
                    min.hl.points=3,
                    min.span.ratio=2,
                    max.aucinf.pext=20,
                    min.hl.r.squared=0.9,
                    tau.choices=NA,
                    single.dose.aucs=check.interval.specification(
                      data.frame(
                        start=0,
                        end=c(24, Inf),
                        auclast=c(TRUE, FALSE),
                        aucinf=c(FALSE, TRUE),
                        half.life=c(FALSE, TRUE),
                        tmax=c(FALSE, TRUE),
                        cmax=c(FALSE, TRUE)))))

  ## Check all the checks on options

  ## adj.r.squared.factor
  expect_error(PKNCA.options(adj.r.squared.factor=c(0.1, 0.9), check=TRUE),
               regexp="adj.r.squared.factor must be a scalar")
  expect_error(PKNCA.options(adj.r.squared.factor=1, check=TRUE),
               regexp="adj.r.squared.factor must be between 0 and 1, exclusive")
  expect_error(PKNCA.options(adj.r.squared.factor=0, check=TRUE),
               regexp="adj.r.squared.factor must be between 0 and 1, exclusive")
  expect_error(PKNCA.options(adj.r.squared.factor="A", check=TRUE),
               regexp="adj.r.squared.factor must be numeric \\(and not a factor\\)")
  expect_warning(v1 <- PKNCA.options(adj.r.squared.factor=0.9, check=TRUE))
  expect_equal(v1, 0.9)
  expect_warning(PKNCA.options(adj.r.squared.factor=0.9, check=TRUE),
                 regexp="adj.r.squared.factor is usually <0.01")

  ## max.missing
  expect_error(PKNCA.options(max.missing=c(1, 2), check=TRUE),
               regexp="max.missing must be a scalar")
  expect_error(PKNCA.options(max.missing="A", check=TRUE),
               regexp="max.missing must be numeric \\(and not a factor\\)")
  expect_error(PKNCA.options(max.missing=-1, check=TRUE),
               regexp="max.missing must be between 0 and 1")
  expect_error(PKNCA.options(max.missing=1, check=TRUE),
               regexp="max.missing must be between 0 and 1")
  expect_equal(PKNCA.options(max.missing=0, check=TRUE),
               0)
  expect_equal(PKNCA.options(max.missing=0.2, check=TRUE),
               0.2)
  expect_warning(PKNCA.options(max.missing=0.6, check=TRUE),
                 regexp="max.missing is usually <= 0.5")

  ## auc.method
  ## All possible methods
  expect_equal(PKNCA.options(auc.method="linear", check=TRUE),
               "linear",
               info="auc.method selection works for linear")
  expect_equal(PKNCA.options(auc.method="lin up/log down", check=TRUE),
               "lin up/log down",
               info="auc.method selection works for lin up/log down")
  expect_error(PKNCA.options(auc.method="foo", check=TRUE),
               regex="should be one of",
               info="auc.method is a valid method")

  ## conc.na
  expect_equal(PKNCA.options(conc.na="drop", check=TRUE),
               "drop")
  expect_warning(v1 <- PKNCA.options(conc.na=factor("drop"), check=TRUE),
                 regexp="conc.na may not be a factor; attempting conversion")
  expect_equal(v1, "drop")
  expect_warning(PKNCA.options(conc.na=factor("drop"), check=TRUE),
                 regexp="conc.na may not be a factor; attempting conversion")
  expect_equal(PKNCA.options(conc.na=1, check=TRUE),
               1)
  expect_warning(v1 <- PKNCA.options(conc.na=-1, check=TRUE),
                 regexp="conc.na is usually not < 0")
  expect_equal(v1, -1)
  expect_warning(PKNCA.options(conc.na=-1, check=TRUE),
                 regexp="conc.na is usually not < 0")
  expect_error(PKNCA.options(conc.na=Inf, check=TRUE),
               regexp="When a number, conc.na must be finite")
  expect_error(PKNCA.options(conc.na="foo", check=TRUE),
               regexp="conc.na must either be a finite number or the text 'drop'")

  ## conc.blq
  ## Confirm all types of single-style inputs
  expect_equal(PKNCA.options(conc.blq="drop", check=TRUE),
               "drop")
  expect_equal(PKNCA.options(conc.blq="keep", check=TRUE),
               "keep")
  expect_equal(PKNCA.options(conc.blq=0, check=TRUE),
               0)
  expect_equal(PKNCA.options(conc.blq=1, check=TRUE),
               1)
  expect_warning(v1 <- PKNCA.options(conc.blq=factor("drop"), check=TRUE),
                 "conc.blq may not be a factor; attempting conversion")
  expect_equal(v1, "drop")
  expect_error(PKNCA.options(conc.blq="foo", check=TRUE),
               regexp="conc.blq must either be a finite number or the text 'drop' or 'keep'")
  expect_error(PKNCA.options(conc.blq=c(1, 2), check=TRUE),
               regexp="conc.blq must be a scalar")
  expect_error(PKNCA.options(conc.blq=NA, check=TRUE),
               regexp="conc.blq must not be NA")

  ## Confirm that list-style input also works
  expect_equal(PKNCA.options(conc.blq=list(first="drop", middle=5, last="keep"),
                             check=TRUE),
               list(first="drop", middle=5, last="keep"))
  expect_error(PKNCA.options(conc.blq=list(first="drop", middle=5, last="keep",
                               foo=5),
                             check=TRUE),
               regexp="When given as a list, conc.blq must only have elements named 'first', 'middle', and 'last'.")
  expect_error(PKNCA.options(conc.blq=list(first="drop", middle=5),
                             check=TRUE),
               regexp="When given as a list, conc.blq must include elements named 'first', 'middle', and 'last'.")

  ## first.tmax
  expect_equal(PKNCA.options(first.tmax=FALSE, check=TRUE),
               FALSE)
  expect_error(PKNCA.options(first.tmax=c(FALSE, TRUE), check=TRUE),
               regexp="first.tmax must be a scalar")
  ## Conversion works
  expect_warning(v1 <- PKNCA.options(first.tmax="T", check=TRUE),
                 regexp="Converting first.tmax to a logical value: TRUE")
  expect_equal(v1, TRUE)
  expect_warning(v1 <- PKNCA.options(first.tmax=1, check=TRUE),
                 regexp="Converting first.tmax to a logical value: TRUE")
  expect_equal(v1, TRUE)
  expect_error(PKNCA.options(first.tmax=NA, check=TRUE),
               regexp="first.tmax may not be NA")
  expect_error(PKNCA.options(first.tmax="x", check=TRUE),
               regexp="Could not convert first.tmax to a logical value")

  ## min.hl.points
  expect_equal(PKNCA.options(min.hl.points=3, check=TRUE),
               3)
  expect_error(PKNCA.options(min.hl.points=c(3, 4), check=TRUE),
               regexp="min.hl.points must be a scalar")
  expect_error(PKNCA.options(min.hl.points=factor(3), check=TRUE),
               regexp="min.hl.points cannot be a factor")
  expect_error(PKNCA.options(min.hl.points="a", check=TRUE),
               regexp="min.hl.points must be a number")
  expect_error(PKNCA.options(min.hl.points=1.5, check=TRUE),
               regexp="min.hl.points must be >=2")
  expect_warning(v1 <- PKNCA.options(min.hl.points=2.5, check=TRUE),
                 regexp="Non-integer given for min.hl.points; rounding to nearest integer")
  ## Note that R uses the engineer's rule of rounding
  expect_equal(v1, 2)

  ## min.span.ratio
  expect_equal(PKNCA.options(min.span.ratio=2, check=TRUE),
               2)
  expect_error(PKNCA.options(min.span.ratio=0, check=TRUE),
               regexp="min.span.ratio must be > 0")
  expect_error(PKNCA.options(min.span.ratio=c(2, 1), check=TRUE),
               regexp="min.span.ratio must be a scalar")
  expect_error(PKNCA.options(min.span.ratio=factor(1), check=TRUE),
               regexp="min.span.ratio cannot be a factor")
  expect_error(PKNCA.options(min.span.ratio="a", check=TRUE),
               regexp="min.span.ratio must be a number")
  expect_warning(PKNCA.options(min.span.ratio=1, check=TRUE),
                 regexp="min.span.ratio is usually >= 2")

  ## max.aucinf.pext
  expect_equal(PKNCA.options(max.aucinf.pext=20, check=TRUE),
               20)
  expect_error(PKNCA.options(max.aucinf.pext=0, check=TRUE),
               regexp="max.aucinf.pext must be > 0")
  expect_error(PKNCA.options(max.aucinf.pext=c(2, 1), check=TRUE),
               regexp="max.aucinf.pext must be a scalar")
  expect_error(PKNCA.options(max.aucinf.pext=factor(1), check=TRUE),
               regexp="max.aucinf.pext cannot be a factor")
  expect_error(PKNCA.options(max.aucinf.pext="a", check=TRUE),
               regexp="max.aucinf.pext must be a number")
  expect_warning(PKNCA.options(max.aucinf.pext=25.1, check=TRUE),
                 regexp="max.aucinf.pext is usually <=25")
  expect_warning(PKNCA.options(max.aucinf.pext=0.1, check=TRUE),
                 regexp="max.aucinf.pext is on the percent not ratio scale, value given is <1%")

  ## min.hl.r.squared
  expect_equal(PKNCA.options(min.hl.r.squared=0.9, check=TRUE),
               0.9)
  expect_error(PKNCA.options(min.hl.r.squared=0, check=TRUE),
               regexp="min.hl.r.squared must be between 0 and 1, exclusive")
  expect_error(PKNCA.options(min.hl.r.squared=c(2, 1), check=TRUE),
               regexp="min.hl.r.squared must be a scalar")
  expect_error(PKNCA.options(min.hl.r.squared=factor(1), check=TRUE),
               regexp="min.hl.r.squared cannot be a factor")
  expect_error(PKNCA.options(min.hl.r.squared="a", check=TRUE),
               regexp="min.hl.r.squared must be a number")
  expect_warning(PKNCA.options(min.hl.r.squared=0.89, check=TRUE),
                 regexp="min.hl.r.squared is usually >= 0.9")

  ## tau.choices
  expect_equal(PKNCA.options(tau.choices=NA, check=TRUE),
               NA)
  expect_equal(PKNCA.options(tau.choices=c(1, 2), check=TRUE),
               c(1, 2))
  expect_error(PKNCA.options(tau.choices=c(NA, 1), check=TRUE),
               regexp="tau.choices may not include NA and be a vector")
  expect_error(PKNCA.options(tau.choices="x", check=TRUE),
               regexp="tau.choices must be a number")

  ## Reset all options to their default to ensure that any subsequent
  ## tests work correctly.
  PKNCA.options(default=TRUE)
})

test_that("PKNCA.choose.option", {
  current.options <- PKNCA.options()
  ## If nothing is given for the non-default options, the default
  ## option is returned.
  expect_equal(PKNCA.choose.option("conc.na"),
               current.options[["conc.na"]])
  ## If an invalid option is requested, it gives an error
  expect_error(PKNCA.choose.option("foo"),
               regexp="PKNCA.options does not have value\\(s\\) for foo.")
  ## It gives an error even if there is a value in the option list.
  ## Note that the error is different because it checks the passed-in
  ## options while it just extracts the default options.
  expect_error(PKNCA.choose.option("foo", options=list(foo="bar")),
               regexp="Invalid setting for PKNCA: foo")
  ## When given in the options list, it will choose that instead of
  ## the default value.
  expect_equal(PKNCA.choose.option("max.aucinf.pext",
                                   options=list(max.aucinf.pext=10)),
               10)
  ## When multiple values are given in the options, it chooses the
  ## right one and ignores all the others (so invalid options can be
  ## listed as long as they are not used).
  expect_equal(PKNCA.choose.option("max.aucinf.pext",
                                   options=list(
                                     foo="bar",
                                     max.aucinf.pext=10)),
               10)
})

context("PKNCA summary setting")

test_that(
  "PKNCA.set.summary input checking",
  {
    ## Get the current state to reset it at the end
    initial.summary.set <- PKNCA.set.summary()
    PKNCA.set.summary(reset=TRUE)
    ## Confirm that reset actually resets the summary settings
    expect_equal(PKNCA.set.summary(), list())
    
    ## name must already be defined
    expect_error(PKNCA.set.summary("blah"),
                 regexp="You must first define the parameter name with add.interval.col")
    ## point must be a function
    expect_error(PKNCA.set.summary("auclast", point="a"),
                 regexp="point must be a function")
    ## spread must be a function
    expect_error(PKNCA.set.summary("auclast", point=mean, spread="a"),
                 regexp="spread must be a function")
    ## Rounding must either be a function or a list
    expect_error(PKNCA.set.summary("auclast", point=mean, spread=sd,
                                   rounding="a"),
                 regexp="rounding must be either a list or a function")
    expect_error(PKNCA.set.summary("auclast", point=mean, spread=sd,
                                   rounding=list(foo=3, bar=4)),
                 regexp="rounding must have a single value in the list")
    expect_error(PKNCA.set.summary("auclast", point=mean, spread=sd,
                                   rounding=list(foo=3)),
                 regexp="When a list, rounding must have a name of either 'signif' or 'round'")
    ## An initial setting works
    expect_equal(PKNCA.set.summary("auclast", point=mean, spread=sd,
                                   rounding=round),
                 list(auclast=list(point=mean, spread=sd, rounding=round)))
    ## Changing a setting works
    expect_equal(PKNCA.set.summary("auclast", point=mean, spread=sd,
                                   rounding=list(round=2)),
                 list(auclast=list(point=mean, spread=sd,
                        rounding=list(round=2))))

    ## Reset all the values to the defaults
    PKNCA.set.summary(reset=TRUE)
    for (n in names(initial.summary.set)) {
      tmp <- initial.summary.set[[n]]
      tmp$name <- n
      do.call(PKNCA.set.summary, tmp)
    }
  })
