#' Get the dependent variable (left hand side of the formula) from a
#' PKNCA object.
#'
#' @param x The object to extract the formula from
#' @param \dots Unused
#' @return The vector of the dependent variable from the object.
#' @export
getDepVar <- function(x, ...)
  UseMethod("getDepVar", x)

#' Get the independent variable (right hand side of the formula) from
#' a PKNCA object.
#'
#' @param x The object to extract the formula from
#' @param \dots Unused
#' @return The vector of the independent variable from the object.
#' @export
getIndepVar <- function(x, ...)
  UseMethod("getIndepVar", x)

#' Get the value from a column in a data frame if the value is a column 
#' there, otherwise, the value should be a scalar or the length of the 
#' data.
#' 
#' @param data A data.frame or similar object
#' @param value A character string giving the name of a column in the 
#'   \code{data}, a scalar, or a vector the same length as the 
#'   \code{data}
#' @param prefix The prefix to use if a column must be added (it will be
#'   used as the full column name if it is not already in the dataset or
#'   it will be prepended to the maximum column name if not.)
#' @return A list with elements named "data", "name" giving the 
#'   \code{data} with a column named "name" with the value in that 
#'   column.
getColumnValueOrNot <- function(data, value, prefix="X") {
  col.name <- setdiff(c(prefix, paste(prefix, max(names(data)), sep=".")), names(data))[1]
  if (is.character(value) && length(value) == 1 && (value %in% names(data))) {
    ## It was a column from the data.frame
    ret <- list(data=data, name=value)
  } else if (length(value) %in% c(1, nrow(data))) {
    data[[col.name]] <- value
    ret <- list(data=data, name=col.name)
  } else {
    stop("value was not a column name nor was it a scalar or a vector matching the length of the data.")
  }
  ret
}
