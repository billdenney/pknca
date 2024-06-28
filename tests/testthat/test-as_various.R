test_that("as_PKNCA*", {
  tmpconc <- generate.conc(2, 1, 0:24)
  tmpdose <- generate.dose(tmpconc)
  myconc <- PKNCAconc(tmpconc, formula=conc~time|treatment+ID)
  mydose <- PKNCAdose(tmpdose, formula=dose~time|treatment+ID)
  mydata <- PKNCAdata(myconc, mydose)
  myresult <- pk.nca(mydata)

  # as_PKNCAconc
  expect_equal(as_PKNCAconc(myconc), myconc)
  expect_error(as_PKNCAconc(mydose))
  expect_equal(as_PKNCAconc(mydata), myconc)
  expect_equal(as_PKNCAconc(myresult), myconc)

  # as_PKNCAdose
  expect_error(as_PKNCAdose(myconc))
  expect_equal(as_PKNCAdose(mydose), mydose)
  expect_equal(as_PKNCAdose(mydata), mydose)
  expect_equal(as_PKNCAdose(myresult), mydose)

  # as_PKNCAdata
  expect_equal(as_PKNCAdata(mydata), mydata)
  # This is not the same as the input data because options have been applied
  expect_equal(as_PKNCAdata(myresult), myresult$data)

  # as_PKNCAresults
  # Ignore the provenance as that may differ due to recalculation and loaded packages
  expect_equal(as_PKNCAresults(mydata), myresult, ignore_attr = TRUE)
  expect_equal(as_PKNCAresults(myresult), myresult)
})
