context("Class generation-general")

test_that("getColumnValueorNot", {
  tmpdata <- data.frame(a=1:2, b=3:4)
  expect_equal(getColumnValueOrNot(tmpdata, "a", "d"),
               list(data=tmpdata, name="a"))
  expect_equal(getColumnValueOrNot(tmpdata, "d", "d"),
               list(data=cbind(tmpdata,
                               data.frame(d="d", stringsAsFactors=FALSE)),
                    name="d"))
  expect_error(getColumnValueOrNot(tmpdata, 1:3, "d"),
               regexp="value was not a column name nor was it a scalar or a vector matching the length of the data")
})
