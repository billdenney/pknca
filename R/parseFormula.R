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
#' \describe{
#'   \item{model}{The left~right side of the model (excluding groups)}
#'   \item{lhs}{The call for the left hand side}
#'   \item{rhs}{The call for the right hand side (excluding groups)}
#'   \item{groups}{The call for the groups}
#'   \item{groupFormula}{A formula form of the groups}
#'   \item{env}{The original formula's environment}
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
  if (!inherits(form, "formula")) {
    made.formula <- FALSE
    try({
      form <- stats::as.formula(form)
      made.formula <- TRUE
    }, silent=TRUE)
    if (!made.formula)
      stop("form must be a formula object or coercable into one")
  }
  ## Check how many sides the formula has and extract the left and
  ## right sides
  lhs <- findOperator(form, "~", "left")
  rhs <- findOperator(form, "~", "right")
  groups <- findOperator(rhs, "|", "right")
  if (identical(groups, NULL)) {
    groups <- NA
    grpFormula <- NA
  } else {
    grpFormula <- stats::as.formula(call("~", groups),
                                    env=environment(form))
    rhs <- findOperator(rhs, "|", "left")
  }
  if (require.two.sided &
      identical(lhs, NA))
    stop("formula is one-sided with require.two.sided set to TRUE")
  if (require.groups &
      identical(groups, NA)) {
    stop("rhs of formula must be a conditioning expression")
  }
  if (identical(lhs, NA)) {
    model <- stats::as.formula(call("~", rhs),
                               env=environment(form))
  } else {
    model <- stats::as.formula(call("~", lhs, rhs),
                               env=environment(form))
  }
  ret <-
    list(model = model,
         lhs = lhs,
         rhs = rhs,
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
  cat(deparse(stats::formula(x)), "\n")
}

#' Convert the parsed formula back into the original
#'
#' @param x The parsed formula object to revert to the original
#' @param drop.groups logical. Should the returned formula drop the
#' groups?
#' @param drop.lhs logical. Should the returned formula be one-sided
#' dropping the left hand side?
#' @param \dots Unused.
#' @return A formula (optionally with portions removed)
#' @export
formula.parseFormula <- function(x, drop.groups=FALSE, drop.lhs=FALSE, ...) {
  if (identical(x$lhs, NA) | drop.lhs) {
    ret <- stats::as.formula(call("~", x$rhs))
  } else {
    ret <- stats::as.formula(call("~", x$lhs, x$rhs))
  }
  if (!identical(x$groups, NA) & !drop.groups)
    ret <- stats::as.formula(paste0(deparse(ret), "|", deparse(x$groups)))
  environment(ret) <- x$env
  ret
}

#' Find the first occurrence of an operator in a formula and return
#' the left, right, or both sides of the operator.
#'
#' @param x The formula to parse
#' @param op The operator to search for (e.g. \code{+}, \code{-},
#' \code{*}, \code{/}, ...)
#' @param side Which side of the operator would you like to see:
#' 'left', 'right', or 'both'.
#' @return The side of the operator requested, NA if requesting the
#' left side of a unary operator, and NULL if the operator is not
#' found.
#' @export
findOperator <- function(x, op, side) {
  side <- match.arg(tolower(side),
                    choices=c("left", "right", "both"))
  if (inherits(x, "name")) {
    ## This is a specific variable, we never found the operator going
    ## down this branch of the tree.
    return(NULL)
  } else if (inherits(x, "call") |
             inherits(x, "formula") |
             inherits(x, "(")) {
    ## This is all or part of a formula
    op <- as.name(op)
    if (identical(x[[1]], op)) {
      ## We found the operator
      if (length(x) == 1) {
        stop("call or formula with length 1 found after finding the operator, unknown how to proceed")
      } else if (length(x) == 2) {
        ## Unary operators have a right hand side only
        if (side == "left") {
          return(NA)
        } else if (side == "right") {
          return(x[[2]])
        } else if (side == "both") {
          return(x)
        }
        stop("Unknown side with a found unary operator")
      } else if (length(x) == 3) {
        ## Binary operator
        if (side == "left") {
          return(x[[2]])
        } else if (side == "right") {
          return(x[[3]])
        } else if (side == "both") {
          return(x)
        }
        stop("Unknown side with a found binary operator")
      }
    } else {
      ## Go down the left then right side of the tree
      if (length(x) == 1)
        stop("call or formula with length 1 found without finding the operator, unknown how to proceed")
      ## First search the left side
      ret <- findOperator(x[[2]], op, side)
      if ((identical(ret, NA) |
           is.null(ret)) &
          length(x) == 3)
        ret <- findOperator(x[[3]], op, side)
    }
  } else {
    ## This should not happen-- find the class that the object is
    stop(sprintf("Cannot handle class %s",
         paste(class(x), sep=", ")))
  }
  ret
}
