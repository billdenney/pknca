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

#' Count the number of values that are not NA
#'
#' @param x The object to count non-NA values within.
#' @return A scalar count of the non-NA values.
#' @export
count.non.missing <- function(x)
  sum(!is.na(x))

## Used for setting labels and units
set.name.matching <- function(ret, name, value, data) {
  if (!missing(value)) {
    if (is.null(names(value)))
      stop(paste(name, "must be a named list"))
    if (!(all(names(labels) %in% names(data))))
      stop(paste(name, "names must match data names"))
    ret[[name]] <- value
  }
  ret
}

## Make plotting labels from data, a formula, labels, and units.
make.label <- function(side, data, parsed.formula, labels, units) {
  col.text <- as.character(parsed.formula[[side]])
  label <- col.text
  if (!is.null(labels))
    if (col.text %in% names(labels))
      label <- labels[[col.text]]
  if (!is.null(units))
    if (col.text %in% names(units))
      label <- sprintf("%s (%s)", label, units[[col.text]])
  label
}
