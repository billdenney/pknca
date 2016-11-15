context("Provenance")

test_that("provenance", {
  a <- addProvenance("a")
  # It may take a bit of time between setting the value in addProvenance
  # and checking
  expect_true(as.numeric(difftime(Sys.time(), attr(a, "provenance", exact=TRUE)$datetime, units="secs")) < 2,
              info="A correct time is set in provenance")
  expect_equal(attr(a, "provenance", exact=TRUE)$sysInfo, Sys.info(),
               info="Correct system information is set in provenance")
  expect_equal(attr(a, "provenance", exact=TRUE)$sessionInfo, sessionInfo(),
               info="Correct session information is set in provenance")

  expect_true(checkProvenance(a),
              info="Provenance checking works to find correct provenance")
  b <- "b"
  expect_true(is.na(checkProvenance(b)),
              info="Provenance checking works to find missing provenance")
  attr(b, "provenance") <- attr(a, "provenance")
  expect_false(checkProvenance(b),
               info="Provenance checking works to find incorrect provenance")
})
