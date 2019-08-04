#' Times relative to an event (typically dosing)
#' 
#' @details 
#' 
#' @param time_event A vector of times for events
#' @param time_obs A vector of times for observations
#' @param units Passed to `base::as.numeric.difftime()`
#' @return A data.frame with columns for:
#' \itemize{
#'   \item{event_number_before}{The index of `time_event` that is the last one before `time_obs` or `NA` if none are before.}
#'   \item{event_number_after}{The index of `time_event` that is the first one after `time_obs` or `NA` if none are after.}
#'   \item{time_before}{The minimum time that the current `time_obs` is before a `time_event`, 0 if at least one `time_obs == time_event`.}
#'   \item{time_after}{The minimum time that the current `time_obs` is after a `time_event`, 0 if at least one `time_obs == time_event`.}
#' }
#'
#' `time_after` and `time_before` are calculated where 
#' @export
time_calc <- function(time_event, time_obs, units=NULL)
  UseMethod("time_calc")

time_calc.numeric <- function(time_event, time_obs, units=NULL) {
  if (any(order(na.omit(time_event)) != seq_along(na.omit(time_event)))) {
    stop("`time_event` must be sorted.")
  }
  if (!is.numeric(time_obs)) {
    stop("Both `time_event` and `time_obs` must be the same class (numeric).")
  }
  ret <-
    data.frame(
      event_number_before=
        sapply(
          X=time_obs,
          FUN=function(x)
            max_zero_len(
              which(time_event <= x),
              zero_length=NA_integer_
            )
        ),
      event_number_after=
        sapply(
          X=time_obs,
          FUN=function(x)
            min_zero_len(
              which(time_event >= x),
              zero_length=NA_integer_
            )
        )
    )
  ret$time_after_event <-
    time_obs - time_event[ret$event_number_before]
  ret$time_before_event <-
    time_obs - time_event[ret$event_number_after]
  ret
}

time_calc.POSIXt <- function(time_event, time_obs, units=NULL) {
  if (is.null(units)) {
    stop("`units` must be provided.")
  }
  if (!("POSIXt" %in% class(time_obs))) {
    stop("Both `time_event` and `time_obs` must be the same class (POSIXt).")
  }
  first_event <- min(time_event, na.rm=TRUE)
  time_calc(
    time_event=difftime(time_event, first_event, units=units),
    time_obs=difftime(time_obs, first_event, units=units),
    units=units
  )
}

time_calc.difftime <- function(time_event, time_obs, units=NULL) {
  if (is.null(units)) {
    stop("`units` must be provided.")
  }
  if (!("difftime" %in% class(time_obs))) {
    stop("Both `time_event` and `time_obs` must be the same class (difftime).")
  }
  time_calc(
    time_event=as.numeric(time_event, units=units),
    time_obs=as.numeric(time_obs, units=units)
  )
}
