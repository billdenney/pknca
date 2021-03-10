context("Provenance")

test_that("provenance", {
  a <- addProvenance("a")
  # It may take a bit of time between setting the value in addProvenance
  # and checking
  expect_true(as.numeric(difftime(Sys.time(), attr(a, "provenance", exact=TRUE)$datetime, units="secs")) < 0.1,
              info="A correct time is set in provenance")
  expect_equal(attr(a, "provenance", exact=TRUE)$sysInfo, Sys.info(),
               info="Correct system information is set in provenance")
  expect_equal(attr(a, "provenance", exact=TRUE)$sessionInfo, sessionInfo(),
               info="Correct session information is set in provenance")
  expect_error(addProvenance(a),
               regexp="object already has provenance and the option to replace it was not selected.",
               info="Adding provenance to an object that already has provenance is an error")
  expect_true(
    {
      # Sleep so that there is a difference and it confirms replacement
      Sys.sleep(0.2)
      as.numeric(difftime(Sys.time(),
                          attr(addProvenance(a, replace=TRUE), "provenance", exact=TRUE)$datetime, units="secs")) < 0.1
    },
    info="A correct time is reset in provenance")
  
  expect_true(checkProvenance(a),
              info="Provenance checking works to find correct provenance")
  b <- "b"
  expect_true(is.na(checkProvenance(b)),
              info="Provenance checking works to find missing provenance")
  attr(b, "provenance") <- attr(a, "provenance")
  expect_false(checkProvenance(b),
               info="Provenance checking works to find incorrect provenance")
  
  # Test printing with a fake provenance object
  fakeprov <- list(hash="a",
                   datetime="b",
                   sessionInfo=list(R.version=c(version.string="c")))
  expect_output(print.provenance(fakeprov),
                regexp="Provenance hash a generated on b with c.",
                info="Provenance prints correctly")
  expect_equal(
    print.provenance(fakeprov),
    invisible("Provenance hash a generated on b with c."),
    info="Provenance printing returns an invisible string with the information in it."
  )
})
