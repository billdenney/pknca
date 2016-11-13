context("merging lists")

test_that("merge.splitByData", {
  ## Simple example
  d1 <- doBy::splitBy(
    ~a,
    data=data.frame(a=rep(1:3, each=2),
      b=11:16))
  d2 <- doBy::splitBy(
    ~a,
    data.frame(a=1:3, d=4:6))

  expect_equal(merge.splitByData(d1, d2),
               list(list(data.frame(a=1, b=11:12),
                         data.frame(a=1, d=4)),
                    list(data.frame(a=2, b=13:14),
                         data.frame(a=2, d=5)),
                    list(data.frame(a=3, b=15:16),
                         data.frame(a=3, d=6))),
               check.attributes=FALSE)

  ## named arguments
  expect_equal(merge.splitByData(first=d1, second=d2),
               list(list(first=data.frame(a=1, b=11:12),
                         second=data.frame(a=1, d=4)),
                    list(first=data.frame(a=2, b=13:14),
                         second=data.frame(a=2, d=5)),
                    list(first=data.frame(a=3, b=15:16),
                         second=data.frame(a=3, d=6))),
               check.attributes=FALSE)

  ## Missing value in data frame 1
  d3 <- doBy::splitBy(
    ~a,
    data=data.frame(a=rep(c(1, 3), each=2),
      b=c(11:12, 15:16)))
  d4 <- doBy::splitBy(
    ~a,
    data.frame(a=1:3, d=4:6))

  expect_equal(merge.splitByData(d3, d4),
               list(list(data.frame(a=1, b=11:12),
                         data.frame(a=1, d=4)),
                    list(data.frame(a=3, b=15:16),
                         data.frame(a=3, d=6)),
                    list(NULL,
                         data.frame(a=2, d=5))),
               check.attributes=FALSE)

  ## Missing value in data frame 2
  d5 <- doBy::splitBy(
    ~a,
    data=data.frame(a=rep(1:3, each=2),
      b=11:16))
  d6 <- doBy::splitBy(
    ~a,
    data.frame(a=c(1, 3), d=c(4, 6)))

  expect_equal(merge.splitByData(d5, d6),
               list(list(data.frame(a=1, b=11:12),
                         data.frame(a=1, d=4)),
                    list(data.frame(a=2, b=13:14),
                         NULL),
                    list(data.frame(a=3, b=15:16),
                         data.frame(a=3, d=6))),
               check.attributes=FALSE)

  expect_equal(merge.splitByData(d5, d6, missing.value=5),
               list(list(data.frame(a=1, b=11:12),
                         data.frame(a=1, d=4)),
                    list(data.frame(a=2, b=13:14),
                         5),
                    list(data.frame(a=3, b=15:16),
                         data.frame(a=3, d=6))),
               check.attributes=FALSE,
               info="Allow other values for missing.value.")

  expect_error(merge.splitByData(5),
               regexp="Must give at least two lists to merge",
               info="Require > 1 thing to merge")
  ## The inputs must be splitByData (or at least have groupid
  ## attributes
  expect_error(merge.splitByData(d1, data.frame(a=5)),
               regexp="All inputs must have a groupid attribute")

  d7 <- doBy::splitBy(
    ~a,
    data=data.frame(a=rep(1:3, each=2),
      b=11:16))
  d8 <- doBy::splitBy(
    ~a+d,
    data.frame(a=c(1, 1, 3),
               d=4:6,
               e=7:9))
  expect_equivalent(merge.splitByData(d7, d8),
                    list(list(data.frame(a=1, b=11:12),
                              data.frame(a=1, d=4, e=7)),
                         list(data.frame(a=1, b=11:12),
                              data.frame(a=1, d=5, e=8)),
                         list(data.frame(a=2, b=13:14),
                              NULL),
                         list(data.frame(a=3, b=15:16),
                              data.frame(a=3, d=6, e=9))),
                    info="Multiple grouping parameters")
})

test_that("merge.splitByData", {
  d1a <- list(1, 2, 3)
  d1b <- list(4, 5, 6)
  attr(d1a, 'groupid') <- data.frame(A=1:3)
  attr(d1b, 'groupid') <- data.frame(A=1:3)
  r1 <- list(list(1, 4), list(2, 5), list(3, 6))
  attr(r1, 'groupid') <- data.frame(A=1:3)
  attr(r1, 'sourcemap') <- data.frame(1:3, 1:3)
  names(attr(r1, 'sourcemap')) <- NULL
  expect_equal(merge.splitlist(d1a, d1b), r1,
               info="merging splitlists in the simple form works")

  expect_error(merge.splitlist(d1a),
               regexp="At least two lists must be given",
               info="merging splitlists requires >=2 inputs")
  
  expect_error(merge.splitlist(d1a, 1),
               regexp="All arguments must be lists",
               info="merging splitlists requires list inputs")

  d2 <- list(4, 5, 6)
  attr(d2, 'groupid') <- data.frame()
  expect_error(merge.splitlist(d1a, d2),
               regexp="The number of rows in the 'groupid' attribute must match the length of the list.",
               info="merging splitlists requires groupid that matches the list")

  d3a <- list(1, 2, 3)
  d3b <- list(4, 5)
  attr(d3a, 'groupid') <- data.frame(A=1:3)
  attr(d3b, 'groupid') <- data.frame(A=1:2)
  r3 <- list(list(1, 4), list(2, 5), list(3, NULL))
  attr(r3, 'groupid') <- data.frame(A=1:3)
  attr(r3, 'sourcemap') <- data.frame(1:3, c(1, 2, NA))
  names(attr(r3, 'sourcemap')) <- NULL
  
  r4 <- list(list(4, 1), list(5, 2), list(NULL, 3))
  attr(r4, 'groupid') <- data.frame(A=1:3)
  attr(r4, 'sourcemap') <- data.frame(c(1, 2, NA), 1:3)
  names(attr(r4, 'sourcemap')) <- NULL
  expect_equal(merge.splitlist(d3a, d3b), r3,
               info="Missing values give a NULL output (second value)")
  expect_equal(merge.splitlist(d3b, d3a), r4,
               info="Missing values give a NULL output (first value)")

  r5 <- list(list(D=1, E=4), list(D=2, E=5), list(D=3, E=NULL))
  attr(r5, 'groupid') <- data.frame(A=1:3)
  attr(r5, 'sourcemap') <- data.frame(D=1:3, E=c(1, 2, NA))
  expect_equal(merge.splitlist(D=d3a, E=d3b), r5,
               info="Names are preserved")
})