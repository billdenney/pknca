context("Interval columns")

# Save the original state
original.state <- get("interval.cols", envir=PKNCA:::.PKNCAEnv)

test_that("add.interval.col", {
  # Invalid inputs fail
  expect_error(add.interval.col(name=1),
               regexp="name must be a character string",
               info="interval column name must be a character string")
  expect_error(add.interval.col(name=c("a", "b")),
               regexp="name must have length",
               info="interval column name must be a scalar character string")

  expect_error(add.interval.col(name="a", FUN=c("a", "b")),
               regexp="FUN must have length == 1",
               info="interval column function must be a scalar character string or NA")
  expect_error(add.interval.col(name="a", FUN=1),
               regexp="FUN must be a character string or NA",
               info="interval column function must be a character string or NA")
  
  expect_error(add.interval.col(name="a", FUN=NA, datatype="individual"),
               regexp="Only the 'interval' datatype is currently supported.",
               info="interval column datatype must be 'interval'")
  
  expect_error(add.interval.col(name="a", FUN=NA, datatype="interval", desc=1:2),
               regexp="desc must have length == 1",
               info="interval column description must be a scalar")
  expect_error(add.interval.col(name="a", FUN=NA, datatype="interval", desc=1),
               regexp="desc must be a character string",
               info="interval column description must be a character scalar")
  expect_error(add.interval.col(name="a", FUN="this function does not exist", datatype="interval", desc="test addition"),
               regexp="The function named '.*' is not defined.  Please define the function before calling add.interval.col.",
               info="interval column function must exist (or be NA)")
  
  expect_equal(
    {
      add.interval.col(name="a", FUN=NA, datatype="interval", desc="test addition")
      get("interval.cols", PKNCA::.PKNCAEnv)[["a"]]
    },
    list(FUN=NA,
         values=c(FALSE, TRUE),
         desc="test addition",
         depends=c(),
         datatype="interval"),
    info="interval column assignment works with FUN=NA")
  expect_equal(
    {
      add.interval.col(name="a", FUN="mean", datatype="interval", desc="test addition")
      get("interval.cols", PKNCA::.PKNCAEnv)[["a"]]
    },
    list(FUN="mean",
         values=c(FALSE, TRUE),
         desc="test addition",
         depends=c(),
         datatype="interval"),
    info="interval column assignment works with FUN=a character string")
})

# Reset the original state
assign("interval.cols", original.state, envir=PKNCA:::.PKNCAEnv)