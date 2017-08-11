context("exclude")

test_that("setExcludeColumn", {
  ## exclude argument not given
  expect_equal(setExcludeColumn(list(data=data.frame(a=1),
                                     exclude="fake")),
               list(data=data.frame(a=1),
                    exclude="fake"),
               info="setExcludeColumn does nothing if the exclude name is already given.")
  expect_equal(setExcludeColumn(list(data=data.frame(a=1))),
               list(data=data.frame(a=1, exclude=NA_character_, stringsAsFactors=FALSE),
                    exclude="exclude"),
               info="setExcludeColumn adds a column named exclude")
  expect_equal(setExcludeColumn(list(data=data.frame(a=1, exclude=2))),
               list(data=data.frame(a=1, exclude=2, exclude.exclude=NA_character_, stringsAsFactors=FALSE),
                    exclude="exclude.exclude"),
               info="setExcludeColumn adds a column named exclude.exclude if 'exclude' is already present")
  expect_equal(setExcludeColumn(list(results=data.frame(a=1)), dataname="results"),
               list(results=data.frame(a=1, exclude=NA_character_, stringsAsFactors=FALSE),
                    exclude="exclude"),
               info="setExcludeColumn works with an alternate dataname")
  
  ## exclude argument given
  expect_equal(setExcludeColumn(list(data=data.frame(a=1, exclude=2),
                                     exclude="exclude"),
                                exclude="exclude"),
               list(data=data.frame(a=1, exclude=2),
                    exclude="exclude"),
               info="setExcludeColumn does nothing if exclude is given and matching")
  expect_error(setExcludeColumn(list(data=data.frame(a=1, exclude=2),
                                     exclude="exclude"),
                                exclude="foo"),
               regexp="exclude is already set for the object.",
               info="setExcludeColumn gives an error if exclude is given and not matching")
  expect_error(setExcludeColumn(list(data=data.frame(a=1)),
                                exclude="exclude"),
               regexp="exclude, if given, must be a column name in the input data.",
               info="setExcludeColumn exclude column must be in the data.")
  expect_equal(setExcludeColumn(list(data=data.frame(a=1, exclude=factor("a"))),
                                exclude="exclude"),
               list(data=data.frame(a=1, exclude="a", stringsAsFactors=FALSE),
                    exclude="exclude"),
               info="setExcludeColumn converts factor column to character")
  expect_equal(setExcludeColumn(list(data=data.frame(a=1, exclude=NA, stringsAsFactors=FALSE)),
                                exclude="exclude"),
               list(data=data.frame(a=1, exclude=NA_character_, stringsAsFactors=FALSE),
                    exclude="exclude"),
               info="setExcludeColumn converts logical NA column to character")
  expect_error(setExcludeColumn(list(data=data.frame(a=1, exclude=FALSE, stringsAsFactors=FALSE)),
                                exclude="exclude"),
               regexp="exclude column must be character vector or something convertable to character without loss of information.",
               info="setExcludeColumn gives error on logical non-NA value")
  expect_error(setExcludeColumn(list(data=data.frame(a=1, exclude=5, stringsAsFactors=FALSE)),
                                exclude="exclude"),
               regexp="exclude column must be character vector or something convertable to character without loss of information.",
               info="setExcludeColumn gives error on non-character value")

  # Zero-row data works
  expect_equal(setExcludeColumn(list(data=data.frame())),
               list(data=data.frame(exclude=NA_character_, stringsAsFactors=FALSE)[-1,,drop=FALSE],
                    exclude="exclude"),
               info="setExcludeColumn works with zero-row data")
  expect_equal(setExcludeColumn(list(data=data.frame()), exclude="foo"),
               list(data=data.frame(foo=NA_character_, stringsAsFactors=FALSE)[-1,,drop=FALSE],
                    exclude="foo"),
               info="setExcludeColumn works with zero-row data")
})

test_that("exclude.helper", {
  ## Check inputs
  obj1 <- list(data=data.frame(a=1:5,
                               exclude=NA_character_,
                               stringsAsFactors=FALSE),
               exclude="exclude")
  expect_error(exclude.helper(obj1,
                              reason="Just because"),
               regexp="Either mask for FUN must be given \\(but not both\\).",
               info="One of mask and FUN must be given")
  expect_error(exclude.helper(obj1,
                              reason="Just because",
                              mask=rep(TRUE, 5),
                              FUN=function(x) rep(TRUE, nrow(x$data))),
               regexp="Either mask for FUN must be given \\(but not both\\).",
               info="Both mask and FUN may not be given")
  obj2 <- list(data=data.frame(a=1:5,
                               exclude=NA_character_,
                               stringsAsFactors=FALSE))
  expect_error(exclude.helper(obj2,
                              reason="Just because",
                              mask=rep(TRUE, 5)),
               regexp="object must have an exclude column specified.",
               info="exclude column is required.")
  obj3 <- list(data=data.frame(a=1:5,
                               stringsAsFactors=FALSE),
               exclude="exclude")
  expect_error(exclude.helper(obj3,
                              reason="Just because",
                              mask=rep(TRUE, 5)),
               regexp="exclude column must exist in object\\[\\['data'\\]\\].",
               info="exclude column must exist in the data")
  expect_error(exclude.helper(obj1,
                              reason="Just because",
                              mask=TRUE),
               regexp="mask or the return value from FUN must match the length of the data.",
               info="mask may not be a scalar")
  expect_error(exclude.helper(obj1,
                              reason="Just because",
                              mask=rep(TRUE, 6)),
               regexp="mask or the return value from FUN must match the length of the data.",
               info="mask must match the length of the data.")
  expect_error(exclude.helper(obj1,
                              reason="Just because",
                              FUN=function(x) TRUE),
               regexp="mask or the return value from FUN must match the length of the data.",
               info="The return from FUN may not be a scalar")
  expect_error(exclude.helper(obj1,
                              reason=1:2,
                              FUN=function(x) TRUE),
               regexp="reason must be a scalar.",
               info="Interpretation of a non-scalar reason is unclear")
  expect_error(exclude.helper(obj1,
                              reason=1,
                              FUN=function(x) TRUE),
               regexp="reason must be a character string.",
               info="Interpretation of a non-character reason is unclear")
  
  ## Check operation
  obj4 <- list(data=data.frame(a=1:5,
                               exclude=c(NA_character_, rep("Just because", 4)),
                               stringsAsFactors=FALSE),
               exclude="exclude")
  
  expect_equal(exclude.helper(obj1,
                              reason="Just because",
                              mask=c(FALSE, rep(TRUE, 4))),
               obj4,
               info="Mask given as a vector works")
  expect_equal(exclude.helper(obj1,
                              reason="Just because",
                              FUN=function(x) c(FALSE, rep(TRUE, nrow(x$data)-1))),
               obj4,
               info="A function returning a vector works")

  obj5 <- list(data=data.frame(a=1:5,
                               exclude=c(NA_character_, "Just because", rep("Just because; really", 3)),
                               stringsAsFactors=FALSE),
               exclude="exclude")
  
  expect_equal(
    exclude.helper(
      exclude.helper(obj1,
                     reason="Just because",
                     FUN=function(x) c(FALSE, rep(TRUE, nrow(x$data)-1))),
      reason="really",
      mask=c(FALSE, FALSE, TRUE, TRUE, TRUE)),
    obj5,
    info="Multiple reasons are tracked.")
})