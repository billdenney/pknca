#' Check that the conversion to a data type does not change the number
#' of NA values
#'
#' @param x the value to convert
#' @param FUN the function to use for conversion
#' @param \dots arguments passed to \code{FUN}
#' @return \code{FUN(x, ...)} or an error if the set of NAs change.
#' @export
check.conversion <- function(x, FUN, ...) {
  ret <- FUN(x, ...)
  new.na <- sum(is.na(x) != is.na(ret))
  if (new.na != 0)
    ## FIXME: It would be nice to have it give the function name as
    ## part of the error
    stop(sprintf("%g new NA value(s) created during conversion",
                 new.na))
  ret
}

#' Verify that the concentration and time are valid
#'
#' If the concentrations or times are invalid, will provide an error.
#' Reasons for being invalid are
#' \itemize{
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
#' @param conc Measured concentrations
#' @param time Time of the measurement of the concentrations
#' @param monotonic.time Must the time be unique and monotonically
#' increasing?
#' @return None
#' @export
check.conc.time <- function(conc, time, monotonic.time=TRUE) {
  if (!missing(conc)) {
    if (length(conc) == 0)
      warning("No concentration data given")
    if (any(!is.na(conc) & conc < 0))
      warning("Negative concentrations found")
    if (all(is.na(conc)))
      warning("All concentration data is missing")
  }
  if (!missing(time)) {
    if (any(is.na(time)))
      stop("Time may not be NA")
    if (monotonic.time) {
      if (!all(time[-1] > time[-length(time)]))
        stop("Time must be monotonically increasing")
      if (!(length(time) == length(unique(time))))
        stop("All time values must be unique")
    }
    if (length(time) == 0)
      warning("No time data given")
  }
  if (!missing(conc) & !missing(time)) {
    if (length(conc) != length(time))
      stop("Conc and time must be the same length")
  }
}

#' Similar to lapplyBy but returning a data frame
#'
#' @param formula See \code{splitBy}
#' @param data See \code{splitBy}
#' @param FUN either a function or a named list of functions
#' @return A data frame with one column for each parameter of
#' \code{formula} and one for each \code{FUN}.  If \code{FUN} is a
#' named list, then the columns will be named the same; otherwise, the
#' column will be named "FUN".
#' @export
sapplyBy <- function(formula, data=parent.frame(), FUN) {
  sb <- doBy::splitBy(formula, data = data)
  gr <- unique(attr(sb, "grps"))
  ret <- attr(sb, "groupid")
  if (is.function(FUN))
    FUN <- list(FUN=FUN)
  for (n in names(FUN)) {
    ddd <- sapply(sb, FUN[[n]])
    ret <- cbind(ret,
                 doBy::renameCol(data.frame(result=ddd[gr]),
                                 "result", n))
  }
  ret
}

#' Round a value to a defined number of digits printing out trailing
#' zeros, if applicable.
#'
#' @param x The number to round
#' @param digits integer indicating the number of decimal places
#' @return A string with the value
#' @seealso \code{\link{round}}, \code{\link{signifString}}
#' @export
roundString <- function(x, digits=0) {
  if (length(digits) == 1) {
    if (digits < 0) {
      formatC(round(x, digits), format='f', digits=0)
    } else {
      formatC(round(x, digits), format='f', digits=digits)
    }
  } else if (length(x) == length(digits)) {
    mapply(roundString, x, digits)
  } else {
    stop("digits must either be a scalar or the same length as x")
  }
}

#' Round a value to a defined number of significant digits printing
#' out trailing zeros, if applicable.
#'
#' @param x The number to round
#' @param digits integer indicating the number of significant digits
#' @return A string with the value
#' @seealso \code{\link{signif}}, \code{\link{roundString}}
#' @export
signifString <- function(x, digits=6) {
  toplog <- bottomlog <- rep(NA, length(x))
  ## When 0 give the digits as the output
  bottomlog[x %in% 0] <- digits
  ## When missing, NaN, or infinite, set digits to 0
  bottomlog[x %in% c(NA, NaN) |
            is.infinite(x)] <- 0
  ## Otherwise set it to digits orders of magnitude lower than the
  ## current value
  toplog <- log10(abs(x))
  ## When the order of magnitude is an exact log 10, move up one so
  ## that the math works for determing the lower log.
  mask.exact.log <- (toplog %% 1) == 0
  toplog[mask.exact.log] <- toplog[mask.exact.log] + 1
  toplog <- ceiling(toplog)
  bottomlog[is.na(bottomlog)] <- digits-toplog[is.na(bottomlog)]
  ## Find times when rounding increases the toplog and shift up the
  ## bottomlog to a corresponding degree. e.g. x=0.9999 and digits=2
  ## should be 1.0 not 1.00.
  newtoplog <- log10(abs(round(x, digits=bottomlog)))
  mask.exact.log <- (newtoplog %% 1) == 0
  newtoplog[mask.exact.log] <- newtoplog[mask.exact.log] + 1
  newtoplog <- ceiling(newtoplog)
  mask.move.up <- toplog < newtoplog
  bottomlog[mask.move.up] <- bottomlog[mask.move.up] - 1
  ## Do the rounding
  roundString(x, digits=bottomlog)
}
