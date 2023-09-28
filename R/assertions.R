#' Assert that an interval is accurately defined as an interval, and return the interval
#'
#' @param interval Numeric vector of two numbers for the start and end time of
#'   integration
#' @param start The start time of the interval
#' @param end The end time of the interval
#' @return `interval` (or `c(start, end)`)
#' @keywords Internal
assert_intervaltime_single <- function(interval = NULL, start = NULL, end = NULL) {
  if (is.null(interval) & is.null(start) & is.null(end)) {
    stop("One of `interval` or `start` and `end` must be given")
  }
  if (xor(is.null(start), is.null(end))) {
    stop("Both `start` and `end` or neither must be given")
  }
  if (!is.null(interval)) {
    checkmate::assert_numeric(x = interval, sorted = TRUE, unique = TRUE, any.missing = FALSE, len = 2)
    checkmate::assert_number(x = interval[1], na.ok = FALSE, finite = TRUE)
  }

  if (!is.null(start)) {
    # Check start and end
    checkmate::assert_number(start, na.ok = FALSE, finite = TRUE, null.ok = FALSE)
    checkmate::assert_number(end, na.ok = FALSE, finite = FALSE, lower = start, null.ok = FALSE)

    if (is.null(interval)) {
      interval <- c(start, end)
    } else if (start != interval[1]) {
      stop("`start` must be the same as the first value in the interval if both are given: ", start, "!=", interval[1])
    } else if (end != interval[2]) {
      stop("`end` must be the same as the second value in the interval if both are given: ", end, "!=", interval[2])
    }
  }

  interval
}

#' Verify that concentration measurements are valid
#'
#' @param conc Measured concentrations
#' @return `conc` or give an informative error
#' @rdname assert_conc_time
assert_conc <- function(conc) {
  if (length(conc) == 0) {
    rlang::warn(
      message = "No concentration data given",
      class = "pknca_conc_none"
    )
  } else {
    checkmate::assert_numeric(conc)
    if (all(is.na(conc))) {
      rlang::warn(
        message = "All concentration data are missing",
        class = "pknca_conc_all_missing"
      )
    } else if (any(!is.na(conc) & as.numeric(conc) < 0)) {
      # as.numeric(conc) is required for compatibility with units
      warning("Negative concentrations found")
    }
  }
  conc
}

#' Verify that time values are valid
#'
#' @param time Time of the measurement of the concentrations
#' @param sorted_time Must the time be unique and monotonically
#'   increasing?
#' @return `time` or give an informative error
#' @rdname assert_conc_time
assert_time <- function(time, sorted_time = TRUE) {
  if (length(time) == 0) {
    rlang::warn(
      message = "No time data given",
      class = "pknca_time_none"
    )
  } else {
    checkmate::assert_numeric(time, any.missing = FALSE, sorted = sorted_time, unique = sorted_time)
  }
  time
}

#' Verify that the concentrations and times are valid
#'
#' If the concentrations or times are invalid, will provide an error.
#' Reasons for being invalid are
#' \itemize{
#'   \item \code{time} is not a number
#'   \item \code{conc} is not a number
#'   \item Any \code{time} value is NA
#'   \item \code{time} is not monotonically increasing
#'   \item \code{conc} and \code{time} are not the same length
#' }
#'
#' Some cases may generate warnings but allow the data to proceed.
#' \itemize{
#'   \item A negative concentration is often but not always an
#'     error; it will generate a warning.
#' }
#'
#' @return A data.frame with columns named "conc" and "time" or an informative
#'   error
assert_conc_time <- function(conc, time, sorted_time = TRUE) {
  assert_conc(conc)
  assert_time(time, sorted_time = sorted_time)
  checkmate::assert_numeric(conc, len = length(time))
  data.frame(conc = conc, time = time)
}
