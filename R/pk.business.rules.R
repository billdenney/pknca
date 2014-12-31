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
#' @return A version of FUN that can be called with parameters that
#' are checked for missingness (and zeros) with missing (and zeros)
#' removed before the call.  If \code{PKNCA.options("max.missing")} is
#' exceeded, then NA is returned.
#' @seealso PKNCA.options
#' @export
business <- function(FUN, zero.missing=FALSE)
  function(x, ...) {
    mask.missing <- is.na(x) | (zero.missing & (x %in% 0))
    if (sum(mask.missing)/length(x) > PKNCA.options("max.missing"))
      return(NA)
    FUN(x[!mask.missing], ...)
  }

#' Compute the geometric mean, sd, and CV
#'
#' @param x A vector to compute the geometric mean of
#' @param ... Values passed to the mean or standard deviation function
#' @return The scalar value of the geometric mean, geometric standard
#' deviation, or geometric coefficient of variation.
#' @aliases geosd, geocv
#' @export
geomean <- function(x, ...)
  exp(mean(log(x), ...))

geosd <- function(x, ...)
  exp(sd(log(x), ...))

geocv <- function(x, ...)
  sqrt(exp(sd(log(x), ...)^2)-1)*100

#' Generate functions to do the named function (e.g. mean) applying
#' the business rules.
#'
#' @param x vector to be passed to the various functions
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
