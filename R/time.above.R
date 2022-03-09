#' Determine time at or above a set value
#' 
#' Interpolation is performed aligning with \code{PKNCA.options("auc.method")}.
#' Extrapolation outside of the measured times is not yet implemented.  The
#' \code{method} may be changed by giving a named \code{method} argument, as
#' well.
#'
#' @inheritParams pk.calc.auxc
#' @param conc_above The concentration to be above
#' @param ... Extra arguments.  Currently, the only extra argument that is used
#'   is \code{method} as described in the details section.
#' @return the time above the given concentration
#' @export
pk.calc.time_above <- function(conc, time,
                               conc_above,
                               #method=NULL,
                               ...,
                               options=list(),
                               check=TRUE) {
  arglist <- list(...)
  method <- PKNCA.choose.option(name="auc.method", value=arglist$method, options=options)
  if (missing(conc)) {
    stop("conc must be given")
  }
  if (missing(time)) {
    stop("time must be given")
  }
  stopifnot("conc_above must be a scalar"=length(conc_above) == 1)
  stopifnot("conc_above must not be NA"=!is.na(conc_above))
  stopifnot("conc_above must be numeric"=is.numeric(conc_above))
  if (check) {
    check.conc.time(conc, time)
  }
  # Only keep evaluable rows
  data <- data.frame(conc=conc, time=time)[!is.na(conc), , drop=FALSE]
  # Determine intervals where both values are above, the first value is above,
  # or the second value is above
  data$above <- data$conc >= conc_above
  mask_both <- data$above[-1] & data$above[-nrow(data)]
  mask_first <- !data$above[-1] & data$above[-nrow(data)]
  mask_second <- data$above[-1] & !data$above[-nrow(data)]
  conc1 <- data$conc[-nrow(data)]
  conc2 <- data$conc[-1]
  time1 <- data$time[-nrow(data)]
  time2 <- data$time[-1]
  # Calculations
  if (nrow(data) < 2) {
    ret <- structure(NA_real_, exclude="Too few measured concentrations to assess time_above")
  } else if (method %in% 'lin up/log down') {
    stop("'lin up/log down' is not yet implemented")
  } else if (method %in% 'linear') {
    ret <-
      sum((time2 - time1)[mask_both]) +
      sum(
        ((conc_above - conc1)/(conc2 - conc1)*(time2 - time1))[mask_first]
      ) +
      sum(
        ((conc2 - conc_above)/(conc2 - conc1)*(time2 - time1))[mask_second]
      )
  } else {
    # Should be caught by the method assignment above
    stop("Invalid 'method', please report this as a bug: ", method) # nocov
  }
  ret
}
## Add the column to the interval specification
add.interval.col("time_above",
                 FUN="pk.calc.time_above",
                 values=c(FALSE, TRUE),
                 desc="Time above a given concentration")
PKNCA.set.summary(
  name="time_above",
  description="arithmetic mean and standard deviation",
  point=business.mean,
  spread=business.sd
)
