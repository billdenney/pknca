#' Assert that an interval is accurately defined as an interval, and return the interval
#'
#' @param interval Numeric vector of two numbers for the start and end time of
#'   integration
#' @param start The start time of the interval
#' @param end The end time of the interval
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
