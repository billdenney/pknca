test_that("parseFormula", {
  tmp.env <- new.env()
  
  ## Parse a rhs-only formula
  f1 <- as.formula(~a, env=tmp.env)
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
  f2 <- as.formula(a~b, env=tmp.env)
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
  f3 <- as.formula(~a|b, env=tmp.env)
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
  f4 <- as.formula(~a|b/c, env=tmp.env)
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
  f5 <- as.formula(~a|b+c, env=tmp.env)
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
  f6 <- as.formula(~a|b/c+d/e, env=tmp.env)
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
  f7 <- as.formula(a~b+c|d/e+f/g, env=tmp.env)
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

  ## Ensure that things can be coerced into formulas
  f8 <- "a~b|c"
  r8 <- list(model=as.formula(a~b),
             lhs=as.name("a"),
             rhs=as.name("b"),
             groups=as.name("c"),
             groupFormula=as.formula(~c),
             env=tmp.env)
  expect_equal(parseFormula(f8),
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
                regexp="A one-sided formula without groups.\n  ~a")
  expect_output(print(parseFormula("~a|b")),
                regexp="A one-sided formula with groups.\n  ~a \\| b")
  expect_output(print(parseFormula("a~b")),
                regexp="A two-sided formula without groups.\n  a ~ b")
  expect_output(print(parseFormula("a~b|c")),
                regexp="A two-sided formula with groups.\n  a ~ b \\| c")
})
