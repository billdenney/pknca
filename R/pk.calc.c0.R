#' Estimate the concentration at dosing time for an IV bolus dose.
#'
#' @inheritParams assert_conc_time
#' @param time.dose The time when dosing occurred
#' @param method The order of methods to test (see details)
#' @param check Check the `conc` and `time` inputs
#' @returns The estimated concentration at time 0.
#'
#' @details Methods available for interpolation are below, and each
#' has its own specific function.
#'
#' \describe{
#'   \item{`c0`}{If the observed `conc` at `time.dose` is nonzero, return that.  This method should usually be used first for single-dose IV bolus data in case nominal time zero is measured.}
#'   \item{`logslope`}{Compute the semilog line between the first two measured times, and use that line to extrapolate backward to `time.dose`}
#'   \item{`c1`}{Use the first point after `time.dose`}
#'   \item{`cmin`}{Set c0 to cmin during the interval.  This method should usually be used for multiple-dose oral data and IV infusion data.}
#'   \item{`set0`}{Set c0 to zero (regardless of any other data).  This method should usually be used first for single-dose oral data.}
#' }
#' @export
pk.calc.c0 <- function(conc, time, time.dose=0,
                       method=c("c0", "logslope", "c1", "cmin", "set0"),
                       check=TRUE) {
  # Check the inputs
  if (check) {
    assert_conc_time(conc = conc, time = time)
  }
  if (length(time.dose) != 1) {
    stop("time.dose must be a scalar")
  } else if (!is.numeric(time.dose) | is.factor(time.dose)) {
    stop("time.dose must be a number")
  }
  if (is.na(time.dose)) {
    warning("time.dose is NA")
    return(structure(NA_real_, exclude = "dose time is missing"))
  } else if (time.dose > max(time)) {
    warning("time.dose is after all available data")
    return(structure(NA_real_, exclude = "dose time is after all available concentration data"))
  }
  method <- match.arg(method, several.ok=TRUE)
  # Find the value
  ret <- NA
  while (is.na(ret) &
         length(method) > 0) {
    current.method <- method[1]
    method <- method[-1]
    ret <- do.call(
      paste("pk.calc.c0.method", current.method, sep="."),
      args=list(
        conc=conc,
        time=time,
        time.dose=time.dose,
        check=FALSE
      )
    )
  }
  ret
}

#' @describeIn pk.calc.c0 Semilog regress the first and second points
#' after time.dose.  This method will return `NA` if the second
#' `conc` after `time.dose` is 0 or greater than the first.
pk.calc.c0.method.logslope <- function(conc, time, time.dose=0,
                                       check=TRUE) {
  if (check) {
    assert_conc_time(conc = conc, time = time)
  }
  mask.positive.time <- (time > time.dose &
                         !(is.na(conc)))
  positive.time <- time[mask.positive.time]
  if (length(positive.time) < 2)
    return(NA)
  # If there is enough data, proceed to calculate
  mask.1 <- time %in% positive.time[1]
  mask.2 <- time %in% positive.time[2]
  c1 <- conc[mask.1]
  c2 <- conc[mask.2]
  t1 <- time[mask.1]
  t2 <- time[mask.2]
  if (c2 < c1 &
      c2 != 0) {
    exp(log(c1) - (log(c2)-log(c1))/(t2-t1)*(t1 - time.dose))
  } else {
    NA
  }
}

#' @describeIn pk.calc.c0 Use `C0` = `conc[time %in% time.dose]` if it is
#'   nonzero.
pk.calc.c0.method.c0 <- function(conc, time, time.dose=0, check=TRUE) {
  if (check) {
    assert_conc_time(conc = conc, time = time)
  }
  # If there is a non-missing and nonzero concentration measurement
  # at time.dose, that's our answer.
  mask.dose <- (time %in% time.dose &
                !(conc %in% c(NA, 0)))
  if (any(mask.dose)) {
    conc[mask.dose]
  } else {
    NA
  }
}

#' @describeIn pk.calc.c0 Use `C0` = `C1`.
pk.calc.c0.method.c1 <- function(conc, time, time.dose=0, check=TRUE) {
  if (check) {
    assert_conc_time(conc = conc, time = time)
  }
  mask.post.dose <- (time > time.dose &
                     !is.na(conc))
  if (any(mask.post.dose)) {
    conc[mask.post.dose][1]
  } else {
    NA
  }
}

#' @describeIn pk.calc.c0 Use `C0` = 0 (typically used for single dose oral and
#'   IV infusion)
pk.calc.c0.method.set0 <- function(conc, time, time.dose=0, check=TRUE) {
  0
}

#' @describeIn pk.calc.c0 Use `C0` = Cmin (typically used for multiple dose oral
#'   and IV infusion but not IV bolus)
pk.calc.c0.method.cmin <- function(conc, time, time.dose=0, check=TRUE) {
  pk.calc.cmin(conc, check=check)
}

# Add the column to the interval specification
add.interval.col("c0",
                 FUN="pk.calc.c0",
                 values=c(FALSE, TRUE),
                 unit_type="conc",
                 pretty_name="C0",
                 desc="Initial concentration after an IV bolus",
                 depends=NULL)
PKNCA.set.summary(
  name="c0",
  description="geometric mean and geometric coefficient of variation",
  point=business.geomean,
  spread=business.geocv
)
