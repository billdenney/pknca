#' Parse a formula into its component parts.
#'
#' This function extracts the left hand side (\code{lhs}), right hand
#' side (\code{rhs}), groups (\code{groups} and as a formula,
#' \code{grpFormula}), the environment (\code{env}, and the original
#' left/right hand side of the model (\code{model}).
#'
#' @param form the formula to extract into its parts 
#' @param require.groups is it an error not to have groups?
#' @param require.two.sided is it an error to have a one-sided
#' formula?
#' @return A parseGroupFormula class list with elements of
#' \itemize{
#'   \describe{model}{The left~right side of the model (excluding groups)}
#'   \describe{lhs}{The call for the left hand side}
#'   \describe{rhs}{The call for the right hand side (excluding groups)}
#'   \describe{groups}{The call for the groups}
#'   \describe{grpFormula}{A formula form of the groups}
#'   \describe{env}{The original formula's environment}
#' }
#' @examples
#' parseGroupFormula("a~b", require.groups=FALSE)
#' ## parseGroupFormula("a~b", require.groups=TRUE) # This is an error
#' parseGroupFormula("a~b|c")
#' parseGroupFormula("a~b|c")$groups
parseGroupFormula <- function (form,
                               require.groups=TRUE,
                               require.two.sided=TRUE) {
  ## If it is not a formula, make it a formula if possible
  if (class(form) != "formula") {
    made.formula <- FALSE
    try({
      form <- as.formula(form)
      made.formula <- TRUE
    })
    if (!made.formula)
      stop("formula must be coercable into a formula object")
  }
  ## Check how many sides the formula has and extract the left and
  ## right sides
  n.sides <- length(form) - 1
  if (n.sides == 1) {
    if (require.two.sided)
      stop("formula is one-sided with require.two.sided as TRUE")
    lhs <- NA
    rhs <- form[[2]]
  } else if (n.sides == 2) {
    lhs <- form[[2]]
    rhs <- form[[3]]
  }
  ## Extract the environment
  model <- form
  ## Extract the groups
  if (class(rhs) != "call" || rhs[[1]] != as.symbol("|")) {
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
    grpFormula <- as.formula(paste("~", deparse(groups)))
  }
  ret <- 
    list(model = model,
         lhs = form[[2]],
         rhs = model.rhs,
         groups = groups,
         groupFormula = grpFormula,
         env=environment(form))
  class(ret) <- c("parseGroupFormula", class(ret))
  ret
}

print.parseGroupFormula <- function(x, ...) {
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

## Convert the parsed formula back into the original
formula.parseGroupFormula <- function(x, drop.groups=FALSE, drop.lhs=FALSE, ...) {
  if (identical(x$lhs, NA) | drop.lhs) {
    ret <- as.formula(paste0("~", deparse(x$rhs)))
  } else {
    ret <- as.formula(paste0(deparse(x$lhs), "~", deparse(x$rhs)))
  }
  if (!identical(x$groups, NA) & !drop.groups)
    ret <- as.formula(paste0(deparse(ret), "|", deparse(x$groups)))
  ret
}
