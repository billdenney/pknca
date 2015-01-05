#' Run any function with a maximum missing fraction of X and 0s
#' possibly counting as missing.
#'
#' Note that all missing values are removed prior to calling the
#' function.  The function is called with the 
#'
#' @param FUN function to run.  The function is called as FUN(x, ...)
#' with missing values removed.
#' @param zero.missing Are zeros counted as missing?  If \code{TRUE}
#' then include them in the missing count.
#' @param options List of changes to the default
#' \code{\link{PKNCA.options}} for calculations.
#' @param max.missing The maximum fraction of values that can be
#' missing to report a value.  Must be between 0 and 1, inclusive.
#' @return A version of FUN that can be called with parameters that
#' are checked for missingness (and zeros) with missing (and zeros)
#' removed before the call.  If \code{PKNCA.options("max.missing")} is
#' exceeded, then NA is returned.
#' @seealso PKNCA.options
#' @export
business <- function(FUN,
                     zero.missing=FALSE,
                     options=list(),
                     max.missing=PKNCA.choose.option("max.missing", options)) {
  PKNCA.options(max.missing=max.missing, check=TRUE)
  function(x, ...) {
    mask.missing <- is.na(x) | (zero.missing & (x %in% 0))
    if (sum(mask.missing)/length(x) > max.missing)
      return(NA)
    FUN(x[!mask.missing], ...)
  }
}

#' Compute the geometric mean
#'
#' @param x A vector to compute the geometric mean of
#' @param na.rm Should missing values be removed?
#' @return The scalar value of the geometric mean, geometric standard
#' deviation, or geometric coefficient of variation.
#' @export
geomean <- function(x, na.rm=FALSE) {
  if (na.rm)
    x <- na.omit(x)
  if (any(is.na(x))) {
    as.numeric(NA)
  } else if (any(x == 0)) {
    0
  } else if (any(x < 0)) {
    ## Protect from overflows by using the logarithm
    prod(sign(x))*exp(sum(log(abs(x)))/length(x))
  } else {
    exp(sum(log(x))/length(x))
  }
}

#' describeIn geomean Compute the geometric standard deviation
geosd <- function(x, na.rm=FALSE)
  exp(sd(log(x), na.rm=na.rm))

#' describeIn geomean Compute the geometric coefficient of variation
geocv <- function(x, na.rm=FALSE)
  sqrt(exp(sd(log(x), na.rm=na.rm)^2)-1)*100

#' Generate functions to do the named function (e.g. mean) applying
#' the business rules.
#'
#' @param x vector to be passed to the various functions
#' @param ... Additional arguments to be passed to the underlying
#' function.
#' @return The value of the various functions or NA if too many values
#' are missing
#' @seealso business
#' @export
business.mean <- business(mean)
business.sd <- business(sd)
business.cv <- business(function(x, ...) {100*sd(x, ...)/mean(x, ...)})
business.geomean <- business(geomean, zero.missing=TRUE)
business.geocv <- business(geocv, zero.missing=TRUE)
business.min <- business(min)
business.max <- business(max)
business.median <- business(median)
