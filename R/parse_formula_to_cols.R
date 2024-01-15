#' Find the first occurrence of an operator in a formula and return
#' the left, right, or both sides of the operator.
#'
#' @param x The formula to parse
#' @param op The operator to search for (e.g. `+`, `-`, `*`, `/`, ...)
#' @param side Which side of the operator would you like to see: 'left',
#'   'right', or 'both'.
#' @returns The side of the operator requested, NA if requesting the left side
#'   of a unary operator, and NULL if the operator is not found.
#' @family Formula parsing
findOperator <- function(x, op, side) {
  side <- match.arg(tolower(side),
                    choices=c("left", "right", "both"))
  if (inherits(x, "name")) {
    # This is a specific variable, we never found the operator going
    # down this branch of the tree.
    return(NULL)
  } else if (is.null(x)) {
    return(NULL)
  } else if (inherits(x, "call") |
             inherits(x, "formula") |
             inherits(x, "(")) {
    # This is all or part of a formula
    op <- as.name(op)
    if (identical(x[[1]], op)) {
      # We found the operator
      if (length(x) == 1) {
        stop("call or formula with length 1 found after finding the operator, unknown how to proceed") # nocov
      } else if (length(x) == 2) {
        # Unary operators have a right hand side only
        if (side == "left") {
          return(NA)
        } else if (side == "right") {
          return(x[[2]])
        } else if (side == "both") {
          return(x)
        }
        stop("Unknown side with a found unary operator") # nocov
      } else if (length(x) == 3) {
        # Binary operator
        if (side == "left") {
          return(x[[2]])
        } else if (side == "right") {
          return(x[[3]])
        } else if (side == "both") {
          return(x)
        }
        stop("Unknown side with a found binary operator") # nocov
      }
    } else {
      # Go down the left then right side of the tree
      if (length(x) == 1)
        stop("call or formula with length 1 found without finding the operator, unknown how to proceed")
      # First search the left side
      ret <- findOperator(x[[2]], op, side)
      if ((identical(ret, NA) |
           is.null(ret)) &
          length(x) == 3)
        ret <- findOperator(x[[3]], op, side)
    }
  } else {
    # This should not happen-- find the class that the object is
    stop(sprintf("Cannot handle class %s",
         paste(class(x), sep=", ")))
  }
  ret
}

#' Convert a formula representation to the columns for input data
#'
#' @param form the formula (or something coercible into a formula) to extract
#'   into its parts
#' @returns A list of column names for various formula parts
#' @keywords Internal
#' @family Formula parsing
parse_formula_to_cols <- function(form) {
  if (!inherits(form, "formula")) {
    form <- try({stats::as.formula(form)}, silent = TRUE)
  }
  if (!inherits(form, "formula")) {
    stop("form must be a formula or coercable into one")
  }
  rhs_raw <- findOperator(form, "~", "right")
  groups_raw <- findOperator(rhs_raw, "|", "right")
  if (!is.null(groups_raw)) {
    rhs <- findOperator(rhs_raw, "|", "left")
  } else {
    rhs <- rhs_raw
  }
  # groups_los_raw becomes the last variable to the left of the slash
  groups_raw_c <- all.vars(groups_raw)
  groups_los_raw <- all.vars(findOperator(groups_raw, "/", "left"))
  if (length(groups_los_raw) > 0) {
    groups <- character()
    groups_los <- groups_raw_c[1:which(groups_raw_c == groups_los_raw)]
    groups_ros <- setdiff(groups_raw_c, groups_los)
  } else {
    groups <- groups_raw_c
    groups_los <- character()
    groups_ros <- character()
  }
  ret <-
    list(
      lhs=setdiff(all.vars(findOperator(form, "~", "left")), "."),
      rhs=setdiff(all.vars(rhs), "."),
      groups=groups,
      groups_left_of_slash=groups_los,
      groups_right_of_slash=groups_ros
    )
  ret
}
