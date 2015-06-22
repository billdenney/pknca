library(testthat)

test_that(
  "check.interval.specification", {

    ## Expand a minimal data frame
    d1 <- data.frame(start=0, end=1)
    r1 <- data.frame(start=0,
                     end=1,
                     auc.type=as.character(NA),
                     half.life=FALSE,
                     tfirst=FALSE,
                     tmax=FALSE,
                     tlast=FALSE,
                     cmin=FALSE,
                     cmax=FALSE,
                     clast.obs=FALSE,
                     clast.pred=FALSE,
                     thalf.eff=FALSE,
                     aucpext=FALSE,
                     cl=FALSE,
                     mrt=FALSE,
                     vz=FALSE,
                     vss=FALSE,
                     stringsAsFactors=FALSE)
    expect_equal(check.interval.specification(d1),
                 r1)
    expect_warning(check.interval.specification(d1),
                   regexp="Nothing to be calculated in interval specification number\\(s\\): 1")

    ## AUC specifications that go in as factors come out as characters
    d2 <- data.frame(start=0, end=1,
                     auc.type=factor("auclast"))
    r2 <- data.frame(start=0,
                     end=1,
                     auc.type="auclast",
                     half.life=FALSE,
                     tfirst=FALSE,
                     tmax=FALSE,
                     tlast=FALSE,
                     cmin=FALSE,
                     cmax=FALSE,
                     clast.obs=FALSE,
                     clast.pred=FALSE,
                     thalf.eff=FALSE,
                     aucpext=FALSE,
                     cl=FALSE,
                     mrt=FALSE,
                     vz=FALSE,
                     vss=FALSE,
                     stringsAsFactors=FALSE)
    expect_equal(check.interval.specification(d2),
                 r2)

    ## start and end must both be specified
    d3 <- data.frame(start=0)
    expect_error(check.interval.specification(d3),
                 regexp="Column\\(s\\) 'end' missing from interval specification")
    d4 <- data.frame(end=1)
    expect_error(check.interval.specification(d4),
                 regexp="Column\\(s\\) 'start' missing from interval specification")
    d5 <- data.frame(blah=5)
    expect_error(check.interval.specification(d5),
                 regexp="Column\\(s\\) 'start', 'end' missing from interval specification")

    ## Ensure that there are data
    d6 <- data.frame()
    expect_error(check.interval.specification(d6),
                 regexp="interval specification has no rows")

    ## Confirm specific column values required
    d7 <- data.frame(start=as.numeric(NA), end=1)
    expect_error(
      check.interval.specification(d7),
      regexp="AUC specification may not have NA for the starting time")
    d8 <- data.frame(start=0, end=as.numeric(NA))
    expect_error(
      check.interval.specification(d8),
      regexp="AUC specification may not have NA for the end time")
    d9 <- data.frame(start=1, end=1)
    expect_error(check.interval.specification(d9),
                 regexp="start must be < end")
    d10 <- data.frame(start=1, end=0)
    expect_error(check.interval.specification(d10),
                 regexp="start must be < end")
    d11 <- data.frame(start=Inf, end=1)
    expect_error(check.interval.specification(d11),
                 regexp="start may not be infinite")
    d12 <- data.frame(start=-Inf, end=1)
    expect_error(check.interval.specification(d12),
                 regexp="start may not be infinite")
    ## But it is OK to have an infinite end
    d13 <- data.frame(start=0, end=Inf)
    r13 <- data.frame(start=0,
                     end=Inf,
                     auc.type=as.character(NA),
                     half.life=FALSE,
                     tfirst=FALSE,
                     tmax=FALSE,
                     tlast=FALSE,
                     cmin=FALSE,
                     cmax=FALSE,
                     clast.obs=FALSE,
                     clast.pred=FALSE,
                     thalf.eff=FALSE,
                     aucpext=FALSE,
                     cl=FALSE,
                     mrt=FALSE,
                     vz=FALSE,
                     vss=FALSE,
                     stringsAsFactors=FALSE)
    expect_equal(check.interval.specification(d13), r13)
    
    ## Test that all valid combinations pass through without error or
    ## change
    d13 <- data.frame(start=0, end=Inf)
    r13 <- data.frame(start=0,
                      end=Inf,
                      auc.type=as.character(NA),
                      half.life=FALSE,
                      tfirst=FALSE,
                      tmax=FALSE,
                      tlast=FALSE,
                      cmin=FALSE,
                      cmax=FALSE,
                      clast.obs=FALSE,
                      clast.pred=FALSE,
                      thalf.eff=FALSE,
                      aucpext=FALSE,
                      cl=FALSE,
                      mrt=FALSE,
                      vz=FALSE,
                      vss=FALSE,
                      stringsAsFactors=FALSE)
    expect_equal(check.interval.specification(d13), r13)

    ## All valid auc.type values succeed
    d14 <- data.frame(start=0, end=Inf,
                      auc.type=c("AUClast", "AUCinf", "AUCall", NA))
    r14 <- data.frame(start=0,
                      end=Inf,
                      auc.type=c("AUClast", "AUCinf", "AUCall", NA),
                      half.life=FALSE,
                      tfirst=FALSE,
                      tmax=FALSE,
                      tlast=FALSE,
                      cmin=FALSE,
                      cmax=FALSE,
                      clast.obs=FALSE,
                      clast.pred=FALSE,
                      thalf.eff=FALSE,
                      aucpext=FALSE,
                      cl=FALSE,
                      mrt=FALSE,
                      vz=FALSE,
                      vss=FALSE,
                      stringsAsFactors=FALSE)
    expect_equal(check.interval.specification(d14), r14)

    ## When the no-calculation interval specification is not the first,
    ## ensure that is warned correctly
    expect_warning(check.interval.specification(d14),
                   regexp="Nothing to be calculated in interval specification number\\(s\\): 4")

    ## An invalid auc.type is an error
    d15 <- data.frame(start=0, end=Inf, auc.type="foo")
    expect_error(check.interval.specification(d15),
                 regexp="auc.type must be one of 'aucinf', 'auclast', 'aucall', or NA")

    ## A matrix is converted into a data frame and used
    d16 <- data.frame(start=0, end=1, half.life=c(FALSE, TRUE))
    d16.m <- as.matrix(d16)
    r16 <- data.frame(start=0,
                      end=1,
                      auc.type=as.character(NA),
                      half.life=c(FALSE, TRUE),
                      tfirst=FALSE,
                      tmax=FALSE,
                      tlast=FALSE,
                      cmin=FALSE,
                      cmax=FALSE,
                      clast.obs=FALSE,
                      clast.pred=FALSE,
                      thalf.eff=FALSE,
                      aucpext=FALSE,
                      cl=FALSE,
                      mrt=FALSE,
                      vz=FALSE,
                      vss=FALSE,
                      stringsAsFactors=FALSE)
    expect_equal(check.interval.specification(d16), r16)
    expect_equal(check.interval.specification(d16.m), r16)

    ## Each column that should be logical can also be numeric or
    ## yes/no
    for (n in c("half.life", "tfirst", "tmax", "tlast",
                "cmin", "cmax", "clast.obs", "clast.pred",
                "thalf.eff", "aucpext", "cl", "mrt", "vz", "vss")) {
      d17 <- data.frame(start=0, end=1)
      r17 <- data.frame(
        start=0,
        end=1,
        auc.type=as.character(NA),
        half.life=FALSE,
        tfirst=FALSE,
        tmax=FALSE,
        tlast=FALSE,
        cmin=FALSE,
        cmax=FALSE,
        clast.obs=FALSE,
        clast.pred=FALSE,
        thalf.eff=FALSE,
        aucpext=FALSE,
        cl=FALSE,
        mrt=FALSE,
        vz=FALSE,
        vss=FALSE,
        stringsAsFactors=FALSE)
      d17[,n] <- TRUE
      r17[,n] <- TRUE
      expect_equal(check.interval.specification(d17), r17,
                   label=paste(n, "set to TRUE"))
      d17[,n] <- 1
      expect_equal(check.interval.specification(d17), r17,
                   label=paste(n, "set to 1"))
      d17[,n] <- -1
      expect_equal(check.interval.specification(d17), r17,
                   label=paste(n, "set to -1"))
      d17[,n] <- "yes"
      expect_equal(check.interval.specification(d17), r17,
                   label=paste(n, "set to yes"))
    }

    ## cl can be set to force
    d18 <- data.frame(start=0, end=1, cl='force')
    r18 <- data.frame(
      start=0,
      end=1,
      auc.type=as.character(NA),
      half.life=FALSE,
      tfirst=FALSE,
      tmax=FALSE,
      tlast=FALSE,
      cmin=FALSE,
      cmax=FALSE,
      clast.obs=FALSE,
      clast.pred=FALSE,
      thalf.eff=FALSE,
      aucpext=FALSE,
      cl='force',
      mrt=FALSE,
      vz=FALSE,
      vss=FALSE,
      stringsAsFactors=FALSE)
    expect_equal(check.interval.specification(d18), r18)

    d19 <- data.frame(start=factor(0), end=1)
    expect_error(check.interval.specification(d19),
                 regexp="Must be numeric and not a factor")
})

test_that("make.logical", {
  ## Simple identities
  expect_equal(make.logical(TRUE), TRUE)
  expect_equal(make.logical("TRUE"), TRUE)
  expect_equal(make.logical(factor("TRUE")), TRUE)
  expect_equal(make.logical("T"), TRUE)
  expect_equal(make.logical("YES"), TRUE)
  expect_equal(make.logical("Y"), TRUE)
  expect_equal(make.logical(1), TRUE)
  expect_equal(make.logical(Inf), TRUE)
  expect_equal(make.logical(FALSE), FALSE)
  expect_equal(make.logical("FALSE"), FALSE)
  expect_equal(make.logical(factor("FALSE")), FALSE)
  expect_equal(make.logical("F"), FALSE)
  expect_equal(make.logical("NO"), FALSE)
  expect_equal(make.logical("N"), FALSE)
  expect_equal(make.logical(0), FALSE)
  ## NA conversion works in all types
  expect_equal(make.logical(c(0, NA, 1), na.value=NA),
               c(FALSE, NA, TRUE))
  expect_equal(make.logical(c(0, NA, 1), na.value=FALSE),
               c(FALSE, FALSE, TRUE))
  expect_equal(make.logical(c(0, NA, 1), na.value=TRUE),
               c(FALSE, TRUE, TRUE))
  ## Default na.value is FALSE
  expect_equal(make.logical(c(0, NA, 1)),
               c(FALSE, FALSE, TRUE))
  ## Give a class that cannot be used
  expect_error(make.logical(list()),
               regexp="Cannot handle class: list")
})
