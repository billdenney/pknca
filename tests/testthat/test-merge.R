context("merging lists")

test_that("merge.splitlist", {
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
  expect_error(merge.splitlist(d1a, d2),
               regexp="All arguments must have a 'groupid' attribute",
               info="merging splitlists requires groupid")
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