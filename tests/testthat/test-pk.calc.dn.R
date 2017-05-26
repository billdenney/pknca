context("Dose-normalized NCA functions")

test_that("pk.calc.dn", {
  ## Ensure correct calculation
  expect_equal(pk.calc.dn(1, 5), 0.2)
  expect_equal(pk.calc.dn(NA, 5), NA_real_)
  expect_equal(pk.calc.dn(1, NA), NA_real_)
})

test_that("pk.calc.cmax", {
  # Ensure that the formalsmap functionality works within pk.nca
  source("generate.data.R")
  tmpconc <- generate.conc(2, 2, 0:24)
  tmpdose <- generate.dose(tmpconc)
  tmpconc <- merge(tmpconc, tmpdose[,c("ID", "treatment", "dose")])
  myconc <- PKNCAconc(tmpconc,
                      conc~time|treatment+dose+ID)
  mydose <- PKNCAdose(tmpdose,
                      dose~time|treatment+ID)
  mydata <- PKNCAdata(myconc, mydose,
                      intervals=data.frame(start=0, end=24, auclast.dn=TRUE))
  myres <- pk.nca(mydata)
  expect_equal(myres$result$PPORRES[myres$result$PPTESTCD %in% "auclast"]/
                 myres$result$PPORRES[myres$result$PPTESTCD %in% "auclast.dn"],
               myres$result$dose[myres$result$PPTESTCD %in% "auclast"],
               info="Dose normalization works when requested as a parameter in pk.nca")
})