context("AUC interval checking")

test_that(
  "check.interval.specification", {

    ## Expand a minimal data frame
    d1 <- data.frame(start=0, end=1)
    r1 <- data.frame(start=0,
                     end=1,
                     auclast=FALSE,
                     aucall=FALSE,
                     aumclast=FALSE,
                     aumcall=FALSE,
                     cmax=FALSE,
                     cmin=FALSE,
                     tmax=FALSE,
                     tlast=FALSE,
                     tfirst=FALSE,
                     clast.obs=FALSE,
                     f=FALSE,
                     cav=FALSE,
                     ctrough=FALSE,
                     ptr=FALSE,
                     tlag=FALSE,
                     half.life=FALSE,
                     r.squared=FALSE,
                     adj.r.squared=FALSE,
                     lambda.z=FALSE,
                     lambda.z.time.first=FALSE,
                     lambda.z.n.points=FALSE,
                     clast.pred=FALSE,
                     span.ratio=FALSE,
                     aucinf=FALSE,
                     aumcinf=FALSE,
                     aucpext=FALSE,
                     cl=FALSE,
                     mrt=FALSE,
                     vss=FALSE,
                     vd=FALSE,
                     thalf.eff=FALSE,
                     kel=FALSE,
                     vz=FALSE,
                     stringsAsFactors=FALSE)
    expect_warning(d1.check <- check.interval.specification(d1),
                   regexp="Nothing to be calculated in interval specification number\\(s\\): 1")
    expect_equal(d1.check, r1)

    ## Giving one parameter will fill in everything else as false
    d2 <- data.frame(start=0, end=1, auclast=TRUE)
    r2 <- data.frame(start=0,
                     end=1,
                     auclast=TRUE,
                     aucall=FALSE,
                     aumclast=FALSE,
                     aumcall=FALSE,
                     cmax=FALSE,
                     cmin=FALSE,
                     tmax=FALSE,
                     tlast=FALSE,
                     tfirst=FALSE,
                     clast.obs=FALSE,
                     f=FALSE,
                     cav=FALSE,
                     ctrough=FALSE,
                     ptr=FALSE,
                     tlag=FALSE,
                     half.life=FALSE,
                     r.squared=FALSE,
                     adj.r.squared=FALSE,
                     lambda.z=FALSE,
                     lambda.z.time.first=FALSE,
                     lambda.z.n.points=FALSE,
                     clast.pred=FALSE,
                     span.ratio=FALSE,
                     aucinf=FALSE,
                     aumcinf=FALSE,
                     aucpext=FALSE,
                     cl=FALSE,
                     mrt=FALSE,
                     vss=FALSE,
                     vd=FALSE,
                     thalf.eff=FALSE,
                     kel=FALSE,
                     vz=FALSE,
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
                      auclast=FALSE,
                      aucall=FALSE,
                      aumclast=FALSE,
                      aumcall=FALSE,
                      cmax=FALSE,
                      cmin=FALSE,
                      tmax=FALSE,
                      tlast=FALSE,
                      tfirst=FALSE,
                      clast.obs=FALSE,
                      f=FALSE,
                      cav=FALSE,
                      ctrough=FALSE,
                      ptr=FALSE,
                      tlag=FALSE,
                      half.life=FALSE,
                      r.squared=FALSE,
                      adj.r.squared=FALSE,
                      lambda.z=FALSE,
                      lambda.z.time.first=FALSE,
                      lambda.z.n.points=FALSE,
                      clast.pred=FALSE,
                      span.ratio=FALSE,
                      aucinf=FALSE,
                      aumcinf=FALSE,
                      aucpext=FALSE,
                      cl=FALSE,
                      mrt=FALSE,
                      vss=FALSE,
                      vd=FALSE,
                      thalf.eff=FALSE,
                      kel=FALSE,
                      vz=FALSE,
                      stringsAsFactors=FALSE)
    expect_warning(d13.check <- check.interval.specification(d13))
    expect_equal(d13.check, r13)
    
    ## When the no-calculation interval specification is not the first,
    ## ensure that is warned correctly
    d14 <- data.frame(start=0, end=24, auclast=c(rep(FALSE, 3), TRUE))
    expect_warning(check.interval.specification(d14),
                   regexp="Nothing to be calculated in interval specification number\\(s\\): 1, 2, 3")

    d14 <- data.frame(start=0, end=24, auclast=c(rep(TRUE, 3), FALSE))
    expect_warning(check.interval.specification(d14),
                   regexp="Nothing to be calculated in interval specification number\\(s\\): 4")

    ## Other information is passed through untouched after all the
    ## calculation columns
    d15 <- data.frame(start=0, end=Inf, treatment="foo",
                      stringsAsFactors=FALSE)
    r15 <- data.frame(start=0,
                      end=Inf,
                      auclast=FALSE,
                      aucall=FALSE,
                      aumclast=FALSE,
                      aumcall=FALSE,
                      cmax=FALSE,
                      cmin=FALSE,
                      tmax=FALSE,
                      tlast=FALSE,
                      tfirst=FALSE,
                      clast.obs=FALSE,
                      f=FALSE,
                      cav=FALSE,
                      ctrough=FALSE,
                      ptr=FALSE,
                      tlag=FALSE,
                      half.life=FALSE,
                      r.squared=FALSE,
                      adj.r.squared=FALSE,
                      lambda.z=FALSE,
                      lambda.z.time.first=FALSE,
                      lambda.z.n.points=FALSE,
                      clast.pred=FALSE,
                      span.ratio=FALSE,
                      aucinf=FALSE,
                      aumcinf=FALSE,
                      aucpext=FALSE,
                      cl=FALSE,
                      mrt=FALSE,
                      vss=FALSE,
                      vd=FALSE,
                      thalf.eff=FALSE,
                      kel=FALSE,
                      vz=FALSE,
                      treatment="foo",
                      stringsAsFactors=FALSE)
    expect_warning(v15 <- check.interval.specification(d15))
    expect_equal(v15, r15)

    d16 <- data.frame(start=factor(0), end=1)
    expect_error(check.interval.specification(d16),
                 regexp="Interval column 'start' should not be a factor")
})

test_that("check.interval.deps", {
  ## Confirm that the interval dependencies are accurately added
  expect_equal(check.interval.deps(data.frame(start=0, end=24, aucinf=TRUE)),
               data.frame(start=0,
                          end=24,
                          auclast=FALSE,
                          aucall=FALSE,
                          aumclast=FALSE,
                          aumcall=FALSE,
                          cmax=FALSE,
                          cmin=FALSE,
                          tmax=FALSE,
                          tlast=FALSE,
                          tfirst=FALSE,
                          clast.obs=FALSE,
                          f=FALSE,
                          cav=FALSE,
                          ctrough=FALSE,
                          ptr=FALSE,
                          tlag=FALSE,
                          half.life=TRUE,
                          r.squared=FALSE,
                          adj.r.squared=FALSE,
                          lambda.z=FALSE,
                          lambda.z.time.first=FALSE,
                          lambda.z.n.points=FALSE,
                          clast.pred=FALSE,
                          span.ratio=FALSE,
                          aucinf=TRUE,
                          aumcinf=FALSE,
                          aucpext=FALSE,
                          cl=FALSE,
                          mrt=FALSE,
                          vss=FALSE,
                          vd=FALSE,
                          thalf.eff=FALSE,
                          kel=FALSE,
                          vz=FALSE))
            
          })
