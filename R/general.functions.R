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
#' @param conc Measured concentrations
#' @param time Time of the measurement of the concentrations
#' @param monotonic.time Must the time be unique and monotonically
#' increasing?
#' @return None
#' @export
check.conc.time <- function(conc, time, monotonic.time=TRUE) {
  if (!missing(conc)) {
    if (length(conc) == 0) {
      warning("No concentration data given")
    } else if ((!is.numeric(conc) | is.factor(conc)) &
                   !(is.logical(conc) & all(is.na(conc)))) {
      stop("Concentration data must be numeric and not a factor")
    } else if (all(is.na(conc))) {
      warning("All concentration data is missing")
    } else if (any(!is.na(conc) & conc < 0)) {
      warning("Negative concentrations found")
    }
  }
  if (!missing(time)) {
    if (length(time) == 0) {
      warning("No time data given")
    } else if (any(is.na(time))) {
      stop("Time may not be NA")
    } else if (!is.numeric(time) | is.factor(time)) {
      stop("Time data must be numeric and not a factor")
    }
    if (monotonic.time) {
      if (!all(time[-1] > time[-length(time)]))
        stop("Time must be monotonically increasing")
      if (!(length(time) == length(unique(time))))
        stop("All time values must be unique") # nocov
    }
  }
  if (!missing(conc) & !missing(time)) {
    if (length(conc) != length(time))
      stop("Conc and time must be the same length")
  }
}

#' Round a value to a defined number of digits printing out trailing zeros, if
#' applicable.
#'
#' @param x The number to round
#' @param digits integer indicating the number of decimal places
#' @param sci_range See help for \code{\link{signifString}} (and you likely want
#'   to round with \code{signifString} if you want to use this argument)
#' @param sci_sep The separator to use for scientific notation strings
#'   (typically this will be either "e" or "x10^" for computer- or
#'   human-readable output).
#' @param si_range Deprecated, please use \code{sci_range}
#' @return A string with the value
#' @details Values that are not standard numbers like \code{Inf}, \code{NA}, and
#'   \code{NaN} are returned as \code{"Inf"}, \code{"NA"}, and \code{NaN}.
#' @seealso \code{\link{round}}, \code{\link{signifString}}
#' @export
roundString <- function(x, digits=0, sci_range=Inf, sci_sep="e", si_range) {
  if (!missing(si_range)) {
    .Deprecated(new="roundString with the sci_range argument",
                msg="The si_range argument is deprecated, please use sci_range")
    sci_range <- si_range
  }
  if (length(digits) == 1) {
    mask_na <- is.na(x)
    mask_aschar <- is.nan(x) | is.infinite(x)
    mask_manip <- !(mask_na | mask_aschar)
    ret <- rep(NA, length(x))
    ## Put in the special values
    if (any(mask_na)) {
      ret[mask_na] <- "NA"
    }
    if (any(mask_aschar)) {
      ret[mask_aschar] <- as.character(x[mask_aschar])
    }
    if (any(mask_manip)) {
      xtmp <- round(x[mask_manip], digits)
      mask_sci <-
        xtmp != 0 &
        abs(log10(abs(xtmp))) >= sci_range
      mask_no_sci <- !mask_sci
      if (any(mask_sci)) {
        logval <- floor(log10(abs(xtmp[mask_sci])))
        ret[mask_manip][mask_sci] <-
          paste0(
            formatC(xtmp[mask_sci]/10^logval, format="f", digits=digits + logval),
            sci_sep,
            formatC(logval, format="d"))
      }
      if (any(mask_no_sci)) {
        if (digits < 0) {
          ret[mask_manip][mask_no_sci] <-
            formatC(xtmp[mask_no_sci], format='f', digits=0)
        } else {
          ret[mask_manip][mask_no_sci] <-
            formatC(xtmp[mask_no_sci], format='f', digits=digits)
        }
      }
    }
    ret
  } else if (length(x) == length(digits)) {
    mapply(roundString, x, digits=digits, sci_range=sci_range, sci_sep=sci_sep)
  } else {
    stop("digits must either be a scalar or the same length as x")
  }
}

#' Round a value to a defined number of significant digits printing out trailing
#' zeros, if applicable.
#'
#' @param x The number to round
#' @param digits integer indicating the number of significant digits
#' @param sci_range integer (or \code{Inf}) indicating when to switch to
#'   scientific notation instead of floating point. Zero indicates always use
#'   scientific; \code{Inf} indicates to never use scientific notation;
#'   otherwise, scientific notation is used when \code{abs(log10(x)) > si_range}.
#' @param sci_sep The separator to use for scientific notation strings
#'   (typically this will be either "e" or "x10^" for computer- or
#'   human-readable output).
#' @param si_range Deprecated, please use \code{sci_range}
#' @param ... Arguments passed to methods.
#' @return A string with the value
#' @details Values that are not standard numbers like \code{Inf}, \code{NA}, and
#'   \code{NaN} are returned as \code{"Inf"}, \code{"NA"}, and \code{NaN}.
#' @seealso \code{\link{signif}}, \code{\link{roundString}}
#' @export
signifString <- function(x, ...) 
  UseMethod("signifString")

#' @rdname signifString
#' @export
signifString.data.frame <- function(x, ...) {
  ret <- lapply(x,
                function(y) {
                  if (is.numeric(y) & !is.factor(y)) {
                    signifString(x=y, ...)
                  } else {
                    y
                  }
                })
  ret <- as.data.frame(ret,
                       stringsAsFactors=FALSE)
  rownames(ret) <- rownames(x)
  colnames(ret) <- colnames(x)
  ret
}

#' @rdname signifString
#' @export
signifString.default <- function(x, digits=6, sci_range=6, sci_sep="e", si_range, ...) {
  if (length(list(...))) {
    stop("Additional, unsupported arguments were passed")
  }
  if (!missing(si_range)) {
    .Deprecated(new="roundString with the sci_range argument",
                msg="The si_range argument is deprecated, please use sci_range")
    sci_range <- si_range
  }
  mask_na <- is.na(x)
  mask_aschar <- is.nan(x) | is.infinite(x)
  mask_manip <- !(mask_na | mask_aschar)
  ret <- rep(NA, length(x))
  ## Put in the special values
  if (any(mask_na)) {
    ret[mask_na] <- "NA"
  }
  if (any(mask_aschar)) {
    ret[mask_aschar] <- as.character(x[mask_aschar])
  }
  if (any(mask_manip)) {
    xtmp <- x[mask_manip]
    toplog <- bottomlog <- rep(NA, length(xtmp))
    ## When 0 give the digits as the output
    bottomlog[xtmp %in% 0] <- digits
    ## Otherwise set it to digits orders of magnitude lower than the
    ## current value
    toplog <- log10(abs(xtmp))
    ## When the order of magnitude is an exact log 10, move up one so
    ## that the math works for determing the lower log.
    mask.exact.log <- (toplog %% 1) %in% 0
    toplog[mask.exact.log] <- toplog[mask.exact.log] + 1
    toplog <- ceiling(toplog)
    bottomlog[is.na(bottomlog)] <- digits-toplog[is.na(bottomlog)]
    ## Find times when rounding increases the toplog and shift up the
    ## bottomlog to a corresponding degree. e.g. x=0.9999 and digits=2
    ## should be 1.0 not 1.00.
    newtoplog <- log10(abs(round(xtmp, digits=bottomlog)))
    mask.exact.log <- (newtoplog %% 1) %in% 0
    newtoplog[mask.exact.log] <- newtoplog[mask.exact.log] + 1
    newtoplog <- ceiling(newtoplog)
    mask.move.up <- toplog < newtoplog
    bottomlog[mask.move.up] <- bottomlog[mask.move.up] - 1
    ## Do the rounding
    ret[mask_manip] <- roundString(xtmp, digits=bottomlog,
                                   sci_range=sci_range, sci_sep=sci_sep)
  }
  ret
}

#' Find the summary statistic value with a different value if the input is
#' zero-length.
#'
#' @details If `na.rm` is `TRUE`, then `NA` values will be removed prior to the
#'   check if `length(c(...)) == 0`.
#'
#' @param ... objects to find the summary_statistic for (combined with `c()`
#'   prior to calculation)
#' @param na.rm a logical indicating whether missing values should be removed
#'   (see Details).
#' @param zero_length The value to return if `length(x) == 0`
#' @param FUN the summary statistic function (such as `max` or `min`)
#' @return Either `zero_length` or `FUN(...)`
#' @noRd
#' @importFrom stats na.omit
zero_len_summary <- function(FUN) {
  function(..., na.rm=FALSE, zero_length=NA) { #nocov
    x <- c(...)
    if (na.rm) {
      x <- stats::na.omit(x)
    }
    if (length(x) == 0) {
      zero_length
    } else {
      FUN(x)
    }
  } #nocov
}

#' @describeIn zero_len_summary Find the maximum value with a different value if
#'   the input is zero-length
#' @noRd
max_zero_len <- zero_len_summary(FUN=max)
#' @describeIn zero_len_summary Find the minimum value with a different value if
#'   the input is zero-length
#' @noRd
min_zero_len <- zero_len_summary(FUN=min)
