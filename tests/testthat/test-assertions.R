test_that("assert_intervaltime_single", {
  expect_equal(
    assert_intervaltime_single(interval = c(1, 2)),
    c(1, 2)
  )
  expect_equal(
    assert_intervaltime_single(start = 1, end = 2),
    c(1, 2)
  )
  expect_equal(
    assert_intervaltime_single(interval = c(1, 2), start = 1, end = 2),
    c(1, 2)
  )
  expect_error(
    assert_intervaltime_single(interval = c(1, 2), start = 1.1, end = 2),
    regexp = "`start` must be the same as the first value in the interval if both are given: 1.1!=1",
    fixed = TRUE
  )
  expect_error(
    assert_intervaltime_single(interval = c(1, 2), start = 1, end = 2.1),
    regexp = "`end` must be the same as the second value in the interval if both are given: 2.1!=2",
    fixed = TRUE
  )
})
