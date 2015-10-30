context("Formula parsing")

test_that("parseFormula", {
  tmp.env <- new.env()
  
  ## Parse a rhs-only formula
  f1 <- as.formula("~a", env=tmp.env)
  r1 <- list(model=~a,
             lhs=NA,
             rhs=as.name("a"),
             groups=NA,
             groupFormula=NA,
             env=tmp.env)
  expect_equal(parseFormula(f1),
               r1,
               check.attributes=FALSE)

  ## Parse a lhs and rhs formula
  f2 <- as.formula("a~b", env=tmp.env)
  r2 <- list(model=a~b,
             lhs=as.name("a"),
             rhs=as.name("b"),
             groups=NA,
             groupFormula=NA,
             env=tmp.env)
  expect_equal(parseFormula(f2),
               r2,
               check.attributes=FALSE)

  ## Parse a rhs formula with groups
  f3 <- as.formula("~a|b", env=tmp.env)
  r3 <- list(model=as.formula(~a, env=tmp.env),
             lhs=NA,
             rhs=as.name("a"),
             groups=as.name("b"),
             groupFormula=as.formula(~b, env=tmp.env),
             env=tmp.env)
  expect_equal(parseFormula(f3),
               r3,
               check.attributes=FALSE)

  ## Parse a rhs formula with nested groups
  f4 <- as.formula("~a|b/c", env=tmp.env)
  r4 <- list(model=as.formula(~a, env=tmp.env),
             lhs=NA,
             rhs=as.name("a"),
             groups=call("/", quote(b), quote(c)),
             groupFormula=as.formula(~b/c, env=tmp.env),
             env=tmp.env)
  expect_equal(parseFormula(f4),
               r4,
               check.attributes=FALSE)
  
  ## Parse a rhs formula with crossed groups
  f5 <- as.formula("~a|b+c", env=tmp.env)
  r5 <- list(model=as.formula(~a, env=tmp.env),
             lhs=NA,
             rhs=as.name("a"),
             groups=call("+", quote(b), quote(c)),
             groupFormula=as.formula(~b+c, env=tmp.env),
             env=tmp.env)
  expect_equal(parseFormula(f5),
               r5,
               check.attributes=FALSE)

  ## Parse a rhs formula with crossed and nested groups
  f6 <- as.formula("~a|b/c+d/e", env=tmp.env)
  r6 <- list(model=as.formula(~a, env=tmp.env),
             lhs=NA,
             rhs=as.name("a"),
             groups=call("+",
               call("/", quote(b), quote(c)),
               call("/", quote(d), quote(e))),
             groupFormula=as.formula(~b/c+d/e, env=tmp.env),
             env=tmp.env)
  expect_equal(parseFormula(f6),
               r6,
               check.attributes=FALSE)

  ## Parse a lhs and rhs formula with crossed and nested groups
  f7 <- as.formula("a~b+c|d/e+f/g", env=tmp.env)
  r7 <- list(model=as.formula(a~b+c, env=tmp.env),
             lhs=as.name("a"),
             rhs=call("+", quote(b), quote(c)),
             groups=call("+",
               call("/", quote(d), quote(e)),
               call("/", quote(f), quote(g))),
             groupFormula=as.formula(~d/e+f/g, env=tmp.env),
             env=tmp.env)
  expect_equal(parseFormula(f7),
               r7,
               check.attributes=FALSE)

  ## Ensure that things can be coerced into formulas (ignoring the
  ## environment of the response)
  f8 <- "a~b|c"
  r8 <- list(model=as.formula(a~b),
             lhs=as.name("a"),
             rhs=as.name("b"),
             groups=as.name("c"),
             groupFormula=as.formula(~c))
  expect_equal({t1 <- parseFormula(f8)
                t1$env <- NULL
                t1},
               r8,
               check.attributes=FALSE)

  ## If it can't be made into a formula, get an error.
  expect_error(parseFormula(5),
               regexp="form must be a formula object or coercable into one")

  ## Test the requirements of the formula
  expect_equal(parseFormula(f7, require.groups=TRUE),
               r7,
               check.attributes=FALSE)
  expect_equal(parseFormula(f7, require.two.sided=TRUE),
               r7,
               check.attributes=FALSE)
  expect_equal(parseFormula(f7,
                            require.groups=TRUE,
                            require.two.sided=TRUE),
               r7,
               check.attributes=FALSE)
  expect_error(parseFormula(f2, require.groups=TRUE),
               regex="rhs of formula must be a conditioning expression")
  expect_error(parseFormula(f1, require.two.sided=TRUE),
               regex="formula is one-sided with require.two.sided set to TRUE")
})

test_that("formula.parseFormula", {
  tmp.env <- new.env()
  f1 <- as.formula("~a", env=tmp.env)
  expect_equal(formula(parseFormula(f1)), f1)
  f2 <- as.formula("a~b", env=tmp.env)
  expect_equal(formula(parseFormula(f2)), f2)
  f3 <- as.formula("a+b~c", env=tmp.env)
  expect_equal(formula(parseFormula(f3)), f3)
  f4 <- as.formula("a+b~c+d", env=tmp.env)
  expect_equal(formula(parseFormula(f4)), f4)
  f5 <- as.formula("a+b~c+d|e", env=tmp.env)
  expect_equal(formula(parseFormula(f5)), f5)
  f6 <- as.formula("a+b~c+d|e+f", env=tmp.env)
  expect_equal(formula(parseFormula(f6)), f6)
  f7 <- as.formula("a+b~c+d|e+f/g", env=tmp.env)
  expect_equal(formula(parseFormula(f7)), f7)

  ## Test dropping parts
  r7rhs <- as.formula("~c+d|e+f/g", env=tmp.env)
  expect_equal(formula(parseFormula(f7), drop.lhs=TRUE),
               r7rhs)
  r7g <- as.formula("a+b~c+d", env=tmp.env)
  expect_equal(formula(parseFormula(f7), drop.groups=TRUE),
               r7g)
  r7rhsg <- as.formula("~c+d", env=tmp.env)
  expect_equal(formula(parseFormula(f7),
                       drop.lhs=TRUE,
                       drop.groups=TRUE),
               r7rhsg)

  ## Confirm that the environment is preserved
  expect_false(identical(formula(parseFormula("~a")), f1))
})

test_that("print.parseFormula", {
  expect_output(print(parseFormula("~a")),
                regexp="A one-sided formula without groups[. \n]+~a")
  expect_output(print(parseFormula("~a|b")),
                regexp="A one-sided formula with groups[. \n]+~a \\| b")
  expect_output(print(parseFormula("a~b")),
                regexp="A two-sided formula without groups[. \n]+a ~ b")
  expect_output(print(parseFormula("a~b|c")),
                regexp="A two-sided formula with groups[. \n]+a ~ b \\| c")
})

test_that("findOperator", {
  ## Check the left, right, and both sides of various operators
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

  ## Parentheses are handled unusually-- check that they work
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

  ## Grouping is handled correctly
  f7 <- a+b~c+d|e
  expect_equal(findOperator(f7, "~", "both"), f7)
  expect_equal(findOperator(f7, "|", "both"),
               call("|",
                    call("+", as.name("c"),
                         as.name("d")),
                    as.name("e")))

  ## Grouping with parentheses is handled correctly (as used in
  ## lme4-type formulas)
  f8 <- a+b~c+(d|e)+(f|g/h)
  expect_equal(findOperator(f8, "~", "both"), f8)
  expect_equal(findOperator(f8, "|", "both"),
               call("|",
                    as.name("d"),
                    as.name("e")))
  ## Sub-searching.  Note: that order of operations for formula group
  ## from right to left first.
  expect_equal(
    findOperator(
      findOperator(
        findOperator(f8, "~", "right"),
        "+", "right"),
      "|", "both"),
    call("|",
         as.name("f"),
         call("/", as.name("g"), as.name("h"))))

  ## When things aren't found, we get a NULL back
  f9 <- a~b
  expect_equal(findOperator(f9, "+", "left"), NULL)
  ## When getting the left side of a unary operator, we get an NA back
  f10 <- a~-b
  expect_equal(findOperator(f10, "-", "left"), NA)
  f11 <- a~!b
  expect_equal(findOperator(f11, "!", "left"), NA)
  f12 <- a~+b
  expect_equal(findOperator(f12, "+", "left"), NA)
  f13 <- a~(b)
  expect_equal(findOperator(f13, "(", "left"), NA)

  ## The side can only be "left", "right", or "both"
  expect_error(findOperator(f1, "~", "neither"))

  ## When given an invalid class, return an error
  expect_error(findOperator(list(), side="both"),
               regexp="Cannot handle class list")
})
