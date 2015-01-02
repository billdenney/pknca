#' Parse a formula into its component parts.
#'
#' This function supports parsing 
#'
#' This function extracts the left hand side (\code{lhs}), right hand
#' side (\code{rhs}), groups (\code{groups} and as a formula,
#' \code{grpFormula}), the environment (\code{env}, and the original
#' left/right hand side of the model (\code{model}).
#'
#' This function borrows heavily from the \code{parseGroupFormula}
#' function in the doBy package.
#'
#' @param form the formula to extract into its parts 
#' @param require.groups is it an error not to have groups?
#' @param require.two.sided is it an error to have a one-sided
#' formula?
#' @return A parseFormula class list with elements of
#' \itemize{
#'   \describe{model}{The left~right side of the model (excluding groups)}
#'   \describe{lhs}{The call for the left hand side}
#'   \describe{rhs}{The call for the right hand side (excluding groups)}
#'   \describe{groups}{The call for the groups}
#'   \describe{groupFormula}{A formula form of the groups}
#'   \describe{env}{The original formula's environment}
#' }
#' @examples
#' parseFormula("a~b", require.groups=FALSE)
#' ## parseFormula("a~b", require.groups=TRUE) # This is an error
#' parseFormula("a~b|c")
#' parseFormula("a~b|c")$groups
#' @export
parseFormula <- function (form,
                          require.groups=FALSE,
                          require.two.sided=FALSE) {
  ## If it is not a formula, make it a formula if possible
  if (class(form) != "formula") {
    made.formula <- FALSE
    try({
      form <- as.formula(form)
      made.formula <- TRUE
    }, silent=TRUE)
    if (!made.formula)
      stop("form must be a formula object or coercable into one")
  }
  ## Check how many sides the formula has and extract the left and
  ## right sides
  n.sides <- length(form) - 1
  if (n.sides == 1) {
    if (require.two.sided)
      stop("formula is one-sided with require.two.sided set to TRUE")
    lhs <- NA
    rhs <- form[[2]]
  } else if (n.sides == 2) {
    lhs <- form[[2]]
    rhs <- form[[3]]
  }
  ## Extract the environment
  model <- form
  ## Extract the groups
  if (class(rhs) != "call" ||
      rhs[[1]] != as.symbol("|")) {
    ## If there are no groups
    if (require.groups) {
      stop("rhs of formula must be a conditioning expression")
    } else {
      groups <- NA
      grpFormula <- NA
      model.rhs <- rhs
    }
  } else {
    ## Set the RHS of the model to exclude the groups
    model[[3]] <- rhs[[2]]
    groups <- rhs[[3]]
    model.rhs <- rhs[[2]]
    grpFormula <- as.formula(paste("~", deparse(groups)),
                             env=environment(form))
  }
  if (identical(lhs, NA)) {
    model <- as.formula(call("~", model.rhs), env=environment(form))
  } else {
    model <- as.formula(call("~", lhs, model.rhs), env=environment(form))
  }
  ret <-
    list(model = model,
         lhs = lhs,
         rhs = model.rhs,
         groups = groups,
         groupFormula = grpFormula,
         env=environment(form))
  class(ret) <- c("parseFormula", class(ret))
  ret
}

print.parseFormula <- function(x, ...) {
  if (identical(x$lhs, NA)) {
    cat("A one-sided formula ")
  } else {
    cat("A two-sided formula ")
  }
  if (identical(x$groups, NA)) {
    cat("without groups.\n  ")
  } else {
    cat("with groups.\n  ")
  }
  cat(deparse(formula(x)), "\n")
}

#' Convert the parsed formula back into the original
#'
#' @param x The parsed formula object to revert to the original
#' @param drop.groups logical. Should the returned formula drop the
#' groups?
#' @param drop.lhs logical. Should the returned formula be one-sided
#' dropping the left hand side?
#' @return A formula (optionally with portions removed)
#' @export
formula.parseFormula <- function(x, drop.groups=FALSE, drop.lhs=FALSE, ...) {
  if (identical(x$lhs, NA) | drop.lhs) {
    ret <- as.formula(call("~", x$rhs))
  } else {
    ret <- as.formula(call("~", x$lhs, x$rhs))
  }
  if (!identical(x$groups, NA) & !drop.groups)
    ret <- as.formula(paste0(deparse(ret), "|", deparse(x$groups)))
  environment(ret) <- x$env
  ret
}
