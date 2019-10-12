#' Run any function with a maximum missing fraction of X and 0s possibly
#' counting as missing.  The maximum fraction missing comes from
#' \code{PKNCA.options("max.missing")}.
#'
#' Note that all missing values are removed prior to calling the function.
#' 
#' @param FUN function to run.  The function is called as \code{FUN(x, ...)}
#'   with missing values removed.
#' @param zero.missing Are zeros counted as missing?  If \code{TRUE} then
#'   include them in the missing count.
#' @param max.missing The maximum fraction of the data allowed to be missing (a
#'   number between 0 and 1, inclusive).
#' @return A version of FUN that can be called with parameters that are checked
#'   for missingness (and zeros) with missing (and zeros) removed before the
#'   call.  If \code{max.missing} is exceeded, then NA is returned.
#' @examples
#' my_mean <- pk.business(FUN=mean)
#' mean(c(1:3, NA))
#' # Less than half missing results in the summary statistic of the available
#' # values.
#' my_mean(c(1:3, NA))
#' # More than half missing results in a missing value
#' my_mean(c(1:3, rep(NA, 4)))
#' @export
pk.business <- function(FUN,
                        zero.missing=FALSE,
                        max.missing)
  function(x, ...) {
    ## Allow max.missing to be specified at either function initiation
    ## or for it to use PKNCA.Options("max.missing")
    max.missing <- PKNCA.options("max.missing")
    mask.missing <- is.na(x) | (zero.missing & (x %in% 0))
    if (sum(mask.missing)/length(x) > max.missing)
      return(NA)
    FUN(x[!mask.missing], ...)
  }

#' Compute the geometric mean, sd, and CV
#' 
#' @param x A vector to compute the geometric mean of
#' @param na.rm Should missing values be removed?
#' @return The scalar value of the geometric mean, geometric standard deviation,
#'   or geometric coefficient of variation.
#' @references Kirkwood T. B.L. Geometric means and measures of dispersion.
#'   Biometrics 1979; 35: 908-909
#' @examples
#' geomean(1:3)
#' @export
#' @importFrom stats na.omit
geomean <- function(x, na.rm=FALSE) {
  if (na.rm)
    x <- stats::na.omit(x)
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

#' @describeIn geomean
#' Compute the geometric standard deviation, \code{exp(sd(log(x)))}.
#' @examples
#' geosd(1:3)
#' @export
#' @importFrom stats sd
geosd <- function(x, na.rm=FALSE)
  exp(stats::sd(log(x), na.rm=na.rm))

#' @describeIn geomean
#' Compute the geometric coefficient of variation,
#' \code{sqrt(exp(sd(log(x))^2)-1)*100}.
#' @examples
#' geocv(1:3)
#' @export
#' @importFrom stats sd
geocv <- function(x, na.rm=FALSE)
  sqrt(exp(stats::sd(log(x), na.rm=na.rm)^2)-1)*100

#' Generate functions to do the named function (e.g. mean) applying
#' the business rules.
#'
#' @param x vector to be passed to the various functions
#' @param ... Additional arguments to be passed to the underlying
#' function.
#' @return The value of the various functions or NA if too many values
#' are missing
#' @seealso \code{\link{pk.business}()}
#' @export
business.mean <- pk.business(mean, max.missing=~PKNCA::PKNCA.Options('max.missing'))

#' @describeIn business.mean
#' Compute the standard deviation with business rules.
#' @export
#' @importFrom stats sd
business.sd <- pk.business(stats::sd, max.missing=~PKNCA::PKNCA.Options('max.missing'))

#' @describeIn business.mean
#' Compute the coefficient of variation with business rules.
#' @export
#' @importFrom stats sd
business.cv <-
  pk.business(function(x, ...) {100*stats::sd(x, ...)/mean(x, ...)},
              max.missing=~PKNCA::PKNCA.Options('max.missing'))

#' @describeIn business.mean
#' Compute the geometric mean with business rules.
#' @export
business.geomean <-
  pk.business(geomean, zero.missing=TRUE,
              max.missing=~PKNCA::PKNCA.Options('max.missing'))
#' @describeIn business.mean
#' Compute the geometric coefficient of variation with business rules.
#' @export
business.geocv <-
  pk.business(
    geocv, zero.missing=TRUE,
    max.missing=~PKNCA::PKNCA.Options('max.missing')
  )
#' @describeIn business.mean
#' Compute the minimum with business rules.
#' @export
business.min <-
  pk.business(min, max.missing=~PKNCA::PKNCA.Options('max.missing'))
#' @describeIn business.mean
#' Compute the maximum with business rules.
#' @export
business.max <-
  pk.business(max, max.missing=~PKNCA::PKNCA.Options('max.missing'))
#' @describeIn business.mean
#' Compute the median with business rules.
#' @export
#' @importFrom stats median
business.median <-
  pk.business(stats::median, max.missing=~PKNCA::PKNCA.Options('max.missing'))
#' @describeIn business.mean
#' Compute the range with business rules.
#' @export
business.range <-
  pk.business(range, max.missing=~PKNCA::PKNCA.Options('max.missing'))
