context("splitBy merging")

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
