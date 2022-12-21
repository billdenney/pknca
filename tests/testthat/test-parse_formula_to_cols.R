test_that("findOperator", {
  # Check the left, right, and both sides of various operators
  f1 <- a~b
  expect_equal(findOperator(f1, "~", "both"), f1)
  expect_equal(findOperator(f1, "~", "left"), as.name("a"))
  expect_equal(findOperator(f1, "~", "right"), as.name("b"))

  f2 <- a~b+c
  expect_equal(findOperator(f2, "+", "both"),
               call("+", as.name("b"), as.name("c")))
  expect_equal(findOperator(f2, "+", "left"),
               as.name("b"))
  expect_equal(findOperator(f2, "+", "right"),
               as.name("c"))

  f3 <- a+b~c
  expect_equal(findOperator(f3, "+", "both"),
               call("+", as.name("a"), as.name("b")))
  expect_equal(findOperator(f3, "+", "left"),
               as.name("a"))
  expect_equal(findOperator(f3, "+", "right"),
               as.name("b"))

  f4 <- ~b
  expect_equal(findOperator(f4, "~", "both"), f4)
  expect_equal(findOperator(f4, "~", "left"), NA)
  expect_equal(findOperator(f4, "~", "right"), as.name("b"))

  # Parentheses are handled unusually-- check that they work
  f5 <- a+b~(c+d)
  expect_equal(findOperator(f5, "+", "left"), as.name("a"))
  expect_equal(findOperator(f5, "+", "right"), as.name("b"))

  f6 <- a~(c+d)
  expect_equal(findOperator(f6, "~", "both"), f6)
  expect_equal(findOperator(f6, "+", "left"), as.name("c"))
  expect_equal(findOperator(f6, "+", "right"), as.name("d"))
  expect_equal(findOperator(f6, "(", "left"), NA)
  expect_equal(findOperator(f6, "(", "right"),
               call("+", as.name("c"), as.name("d")))

  # Grouping is handled correctly
  f7 <- a+b~c+d|e
  expect_equal(findOperator(f7, "~", "both"), f7)
  expect_equal(findOperator(f7, "|", "both"),
               call("|",
                    call("+", as.name("c"),
                         as.name("d")),
                    as.name("e")))

  # Grouping with parentheses is handled correctly (as used in
  # lme4-type formulas)
  f8 <- a+b~c+(d|e)+(f|g/h)
  expect_equal(findOperator(f8, "~", "both"), f8)
  expect_equal(findOperator(f8, "|", "both"),
               call("|",
                    as.name("d"),
                    as.name("e")))
  # Sub-searching.  Note: that order of operations for formula group
  # from right to left first.
  expect_equal(
    findOperator(
      findOperator(
        findOperator(f8, "~", "right"),
        "+", "right"),
      "|", "both"),
    call("|",
         as.name("f"),
         call("/", as.name("g"), as.name("h"))))

  # When things aren't found, we get a NULL back
  f9 <- a~b
  expect_equal(findOperator(f9, "+", "left"), NULL)
  # When getting the left side of a unary operator, we get an NA back
  f10 <- a~-b
  expect_equal(findOperator(f10, "-", "left"), NA)
  f11 <- a~!b
  expect_equal(findOperator(f11, "!", "left"), NA)
  f12 <- a~+b
  expect_equal(findOperator(f12, "+", "left"), NA)
  f13 <- a~(b)
  expect_equal(findOperator(f13, "(", "left"), NA)

  # The side can only be "left", "right", or "both"
  expect_error(findOperator(f1, "~", "neither"))

  # When given an invalid class, return an error
  expect_error(findOperator(list(), side="both"),
               regexp="Cannot handle class list")

  expect_error(
    findOperator(x=a~b(), op="+", side="left"),
    regexp="call or formula with length 1 found without finding the operator, unknown how to proceed"
  )
})

test_that("parse_formula_to_cols", {
  expect_equal(
    parse_formula_to_cols(form = "~foo"),
    parse_formula_to_cols(form = ~foo)
  )
  expect_equal(
    parse_formula_to_cols(form = ~foo),
    list(
      lhs = character(),
      rhs = "foo",
      groups = character(),
      groups_left_of_slash = character(),
      groups_right_of_slash = character()
    )
  )
  expect_equal(
    parse_formula_to_cols(form = .~foo),
    list(
      lhs = character(),
      rhs = "foo",
      groups = character(),
      groups_left_of_slash = character(),
      groups_right_of_slash = character()
    )
  )
  expect_equal(
    parse_formula_to_cols(form = foo~.),
    list(
      lhs = "foo",
      rhs = character(),
      groups = character(),
      groups_left_of_slash = character(),
      groups_right_of_slash = character()
    )
  )
  expect_equal(
    parse_formula_to_cols(form = foo~bar),
    list(
      lhs = "foo",
      rhs = "bar",
      groups = character(),
      groups_left_of_slash = character(),
      groups_right_of_slash = character()
    )
  )
  expect_equal(
    parse_formula_to_cols(form = foo~.|bar),
    list(
      lhs = "foo",
      rhs = character(),
      groups = "bar",
      groups_left_of_slash = character(),
      groups_right_of_slash = character()
    )
  )
  expect_equal(
    parse_formula_to_cols(form = foo~baz|bar),
    list(
      lhs = "foo",
      rhs = "baz",
      groups = "bar",
      groups_left_of_slash = character(),
      groups_right_of_slash = character()
    )
  )
  expect_equal(
    parse_formula_to_cols(form = foo~baz|bar+quk),
    list(
      lhs = "foo",
      rhs = "baz",
      groups = c("bar", "quk"),
      groups_left_of_slash = character(),
      groups_right_of_slash = character()
    )
  )
  expect_equal(
    parse_formula_to_cols(form = foo~baz|bar+quk/zot),
    list(
      lhs = "foo",
      rhs = "baz",
      groups = character(),
      groups_left_of_slash = c("bar", "quk"),
      groups_right_of_slash = "zot"
    )
  )
  expect_equal(
    parse_formula_to_cols(form = foo~baz|bar+quk/zot+yap),
    list(
      lhs = "foo",
      rhs = "baz",
      groups = character(),
      groups_left_of_slash = c("bar", "quk"),
      groups_right_of_slash = c("zot", "yap")
    )
  )
})

test_that("parse_formula_to_cols expected errors", {
  expect_error(
    parse_formula_to_cols(form = "foo"),
    regexp = "form must be a formula or coercable into one"
  )
})
