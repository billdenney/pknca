#'  Set Intervals
#'
#'  Takes in two objects, the PKNCAdata object and the proposed intervals.
#'  It will then check that the intervals are valid, given the data object.
#'  If the intervals are valid, it will set them in the object.
#'  It will return the data object with the intervals set.
#'  
#' @param data PKNCAdata object
#' @param intervals Proposed intervals
#' @returns The data object with the intervals set.
#'   
#' @export
set_intervals <- function(data, intervals) {
  valid_intervals <- assert_intervals(intervals, data)
  
  data$intervals <- valid_intervals
  
  data
}

#'  Assert Intervals
#'
#'  Verifies that an interval definition is valid for a PKNCAdata object.
#'  Valid means that intervals are a data.frame (or data.frame-like object), 
#'  that the column names are either the groupings of the PKNCAconc part of 
#'  the PKNCAdata object or that they are one of the NCA parameters allowed 
#'  (i.e. names(get.interval.cols())). 
#'  It will return the intervals argument unchanged, or it will raise an error.
#'  
#' @param intervals Proposed intervals
#' @param data PKNCAdata object
#' @returns The intervals argument unchanged, or it will raise an error.
#'   
#' @export
assert_intervals <- function(intervals, data) {
  if (!is.data.frame(intervals)) {
    stop("The 'intervals' argument must be a data frame or a data frame-like object.")
  }
  
  if (!inherits(data, "PKNCAdata")) {
    stop("The 'data' argument must be a PKNCAdata object.")
  }
  
  allowed_columns <-
    c(
      names(getGroups.PKNCAdata(data)),
      names(get.interval.cols()),
      "conc_above",
      "time_above",
      "impute",
      # If not used, data$options$keep_interval_cols will be NULL
      data$options$keep_interval_cols
    )
  
  invalid_columns <- setdiff(names(intervals), allowed_columns)
  
  if (length(invalid_columns) > 0) {
    stop("The following columns in 'intervals' are not allowed: ", paste(invalid_columns, collapse = ", "))
  }
  
  intervals
}