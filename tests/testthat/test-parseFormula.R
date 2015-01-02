## Note that for all of these tests, we're accepting that the
## environment isn't required.
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
  f4 <- ~a|b/c
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
  f5 <- ~a|b+c
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
  f6 <- ~a|b/c+d/e
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
  f7 <- a~b+c|d/e+f/g
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
})
