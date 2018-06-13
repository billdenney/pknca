context("Class generation-general")

getDataName.list <- function(object) "data"

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

test_that("setAttributeColumn", {
  #skip("Issue #226 in testthat prevents these tests from succeeding")
  obj1 <- structure(list(data=data.frame(A=1:3)),
                    class="PKNCAconc")
  # Validation of inputs
  expect_error(setAttributeColumn(object=obj1,
                                  attr_name=c("A", "B")),
               regexp="attr_name must be a character scalar",
               info="attr_name must be a scalar")
  expect_error(setAttributeColumn(object=obj1,
                                  attr_name=1),
               regexp="attr_name must be a character scalar",
               info="attr_name must be a character")
  expect_error(setAttributeColumn(object=obj1,
                                  attr_name="A",
                                  col_name=c("foo", "A"),
                                  default_value=4),
               regexp="col_name must be a character scalar")
  expect_error(setAttributeColumn(object=obj1,
                                  attr_name="A",
                                  col_name=1,
                                  default_value=4),
               regexp="col_name must be a character scalar")
  expect_error(setAttributeColumn(object=obj1,
                                  attr_name="B",
                                  default_value=c(1:2)),
               regexp="default_value must be a scalar or the same length as the rows in the data")
  # Settings
  expect_message(val1 <- setAttributeColumn(object=obj1,
                                            attr_name="A",
                                            default_value=4),
                 regexp="Found column named A, using it for the attribute of the same name")
  expect_equal(val1,
               structure(list(data=data.frame(A=rep(4, 3)),
                              columns=list(A="A")),
                         class="PKNCAconc"),
               info="col_name defaults to attr_name, column values are not automatically replaced")
  
  # Info provided back with *_if_default
  expect_error(setAttributeColumn(object=obj1,
                                  attr_name="foo",
                                  stop_if_default="bar"),
               regexp="bar",
               info="Error is triggered if default value is used and stop_if_default")
  expect_warning(setAttributeColumn(object=obj1,
                                    attr_name="foo",
                                    warn_if_default="bar"),
                 regexp="bar",
                 info="Error is triggered if default value is used and stop_if_default")
  expect_message(setAttributeColumn(object=obj1,
                                    attr_name="foo",
                                    message_if_default="bar"),
                 regexp="bar",
                 info="Message is triggered if default value is used and stop_if_default")

  expect_equal(setAttributeColumn(object=obj1,
                                  attr_name="foo",
                                  default_value=4),
               structure(list(data=data.frame(A=1:3,
                                              foo=4),
                              columns=list(foo="foo")),
                         class="PKNCAconc"),
               info="col_name defaults to attr_name, column values are added if the column doesn't exist")
  expect_equal(setAttributeColumn(object=obj1,
                                  attr_name="foo",
                                  col_name="bar",
                                  default_value=4),
               structure(list(data=data.frame(A=1:3,
                                              bar=4),
                              columns=list(foo="bar")),
                         class="PKNCAconc"),
               info="attr_name is set to col_name, column values are added if the column exists")
  expect_equal(setAttributeColumn(object=obj1,
                                  attr_name="foo",
                                  col_name="A",
                                  default_value=4),
               structure(list(data=data.frame(A=rep(4, 3)),
                              columns=list(foo="A")),
                         class="PKNCAconc"),
               info="attr_name is set to col_name, column values are not added if the column exists")
  obj2 <- setAttributeColumn(object=obj1,
                             attr_name="foo",
                             col_name="A",
                             default_value=4)
  expect_equal(setAttributeColumn(object=obj2,
                                  attr_name="bar",
                                  col_name="B",
                                  default_value=5),
               structure(list(data=data.frame(A=rep(4, 3),
                                              B=5),
                              columns=list(foo="A",
                                           bar="B")),
                         class="PKNCAconc"),
               info="Adding a second attribute works")
  expect_equal(setAttributeColumn(object=obj2,
                                  attr_name="foo",
                                  col_name="B",
                                  default_value=5),
               structure(list(data=data.frame(A=rep(4, 3),
                                              B=5),
                              columns=list(foo="B")),
                         class="PKNCAconc"),
               info="Overwriting an attribute works and is non-destructive to the existing data")

  # col_or_value testing
  expect_equal(setAttributeColumn(object=obj1,
                                  attr_name="foo",
                                  col_name="A"),
               setAttributeColumn(object=obj1,
                                  attr_name="foo",
                                  col_or_value="A"),
               info="col_or_value assigns to col_name when present")
  expect_equal(setAttributeColumn(object=obj1,
                                  attr_name="foo",
                                  default_value=5),
               setAttributeColumn(object=obj1,
                                  attr_name="foo",
                                  col_or_value=5),
               info="col_or_value assigns to default_value when not present as a col_name")
  expect_error(setAttributeColumn(object=obj1,
                                  attr_name="foo",
                                  col_or_value=5,
                                  default_value=1),
               regexp="Cannot provide col_or_value and col_name or default_value")
  expect_error(setAttributeColumn(object=obj1,
                                  attr_name="foo",
                                  col_or_value=5,
                                  col_name="B"),
               regexp="Cannot provide col_or_value and col_name or default_value")
})

test_that("getAttributeColumn", {
  #skip("Issue #226 in testthat prevents these tests from succeeding")
  obj1 <- structure(list(data=data.frame(A=1:3,
                                         B=4:6),
                         columns=list(foo="A")),
                    class="PKNCAconc")
  obj2 <- structure(list(data=data.frame(A=1:3,
                                         B=4:6),
                         columns=list(foo="A",
                                      bar=c("A", "B"))),
                    class="PKNCAconc")
  obj3 <- structure(list(data=data.frame(A=1:3),
                         columns=list(foo="C")),
                    class="PKNCAconc")
  expect_equal(getAttributeColumn(object=obj1, attr_name="foo"),
               data.frame(A=1:3),
               info="A data frame (not a vector) is returned")
  expect_equal(getAttributeColumn(object=obj2, attr_name="bar"),
               data.frame(A=1:3,
                          B=4:6),
               info="A data frame (not a vector) is returned, and that data frame may have multiple columns")
  expect_warning(val1 <- getAttributeColumn(object=obj1, attr_name="bar"),
                 regexp="bar is not set",
                 info="The user is warned for missing attribute")
  expect_true(is.null(val1),
              info="missing attributes return NULL")
  expect_warning(val1 <- getAttributeColumn(object=obj3, attr_name="foo"),
                 regexp="Columns C are not present",
                 info="The user is warned for missing column")
  # when warnings aren't requested, they should not be given
  expect_silent(getAttributeColumn(object=obj1, attr_name="bar", warn_missing="column"))
  expect_silent(getAttributeColumn(object=obj3, attr_name="foo", warn_missing="attr"))
  expect_silent(getAttributeColumn(object=obj3, attr_name="foo", warn_missing=c()))
  expect_silent(getAttributeColumn(object=obj3, attr_name="bar", warn_missing=c()))
  expect_error(getAttributeColumn(object=obj3, attr_name="foo", warn_missing=c("foo")),
               info="warn_missing must have a valid value")
})

test_that("getDataName.default returns NULL", {
  expect_null(getDataName(1:5),
              info="getDataName.default returns NULL (numeric)")
  expect_null(getDataName("a"),
              info="getDataName.default returns NULL (character)")
  expect_null(getDataName(factor("A")),
              info="getDataName.default returns NULL (factor)")
  expect_null(getDataName(TRUE),
              info="getDataName.default returns NULL (logical)")
})
