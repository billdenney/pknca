#' Select primary NCA parameters for an interval
#'
#' @inheritParams assert_intervaltime_single
#' @param parameters The NCA parameters to calculate during the interval
#' @param groups The groups to use for parameter calculation
#' @param extrapolation When using an extrapolation method (e.g. for the
#'   parameter "aucint"), what extrapolation method should be used?
#' @param clast_extrap When using an extrapolation method, what clast value
#'   should be used:  "obs" = observed, "pred" = predicted
#' @examples
#' make_intervals(start = 0, end = Inf, parameters = c("cmax", "tmax", "half.life", "aucinf"))
#'
#' @export
make_intervals <- function(start, end,
                           parameters = NULL,
                           groups = data.frame(),
                           route = c("extravascular", "intravascular"),
                           dose_duration = c("bolus", "infusion"),
                           dose_type = c("single", "steady-state"),
                           collection = c("point", "interval"),
                           extrapolation = c("inf", "last", "all"),
                           clast_extrap = c("obs", "pred")) {
  interval <- assert_intervaltime_single(start, end)
  if (!is.data.frame(groups)) {
    stop("`groups` must be a data.frame")
  }
  route <- match.arg(route)
  dose_duration <- match.arg(dose_duration)
  dose_type <- match.arg(dose_type)
  collection <- match.arg(collection)
  extrapolation <- match.arg(extrapolation)
  clast_extrap <- match.arg(clast_extrap)
  if (is.null(parameters)) {
    parameters <-
      auto_parameters(
        start = start, end = end,
        route = route,
        dose_duration = dose_duration,
        dose_type = dose_type,
        collection = collection
      )
    # auto-select parameters
    parameters <- c("cmax", "tmax")
    if (end == "inf") {
      if (dose_type == "single") {
        parameters <- c(parameters, "auclast", "aucinf")
      }
      parameters <- c(parameters, "half.life")
    } else {
      # end is not inf
      parameters <- c(parameters, "aucint")
    }
    if (route == "extravascular") {
      if (dose_duration == "bolus") {
        parameters <- c(parameters, "c0")
      }
    }
  }
}

#' Auto-select NCA parameters based on dosing type
#'
#' @inheritParams assert_intervaltime_single
#' @param route What is the route of administration?
#' @param dose_duration If `route = "intravascular"` how is the dose
#'   administered? Ignored if `route = "extravascular"`.
#' @param dose_type Is the data single-dose or steady-state?
#' @param collection Is the collection at a "point" in time (like plasma, serum,
#'   or blood) or an "interval" of time (like urine or feces)?
#' @returns A character vector of parameters to calculate which can be passed to
#'   `make_intervals()`
#' @export
auto_parameters <- function(start, end,
                            route = c("extravascular", "intravascular"),
                            dose_duration = c("bolus", "infusion"),
                            dose_type = c("single", "steady-state"),
                            collection = c("point", "interval")) {
  if (collection == "point") {
    ret <- c("cmax", "tmax")
    if (is.infinite(end)) {
      ret <- c(ret, "half.life")
      if (dose_type == "single") {
        ret <- c(ret, "aucinf.obs", "auclast")
      }
    } else {
      ret <- c(ret, "aucint.inf.obs")
    }
    if (route == "intravascular") {
      if (dose_duration == "bolus") {
        ret <- c(ret, "c0")
      } else if (dose_duration == "infusion") {
        ret <- c(ret, "ceoi")
      }
    }
  } else if (collection == "interval") {
    ret <- c("ae", "fe")
  }
  ret
}
