#' Compute noncompartmental superposition for repeated dosing
#'
#' @inheritParams assert_conc_time
#' @inheritParams assert_lambdaz
#' @inheritParams PKNCA.choose.option
#' @inheritParams choose_interval_method
#' @param dose.input The dose given to generate the `conc` and `time` inputs. If
#'   missing, output doses will be assumed to be equal to the input dose.
#' @inheritParams assert_dosetau
#' @param dose.times The time of dosing within the dosing interval. The
#'   `min(dose.times)` must be >= 0, and the `max(dose.times)` must be < `tau`.
#'   There may be more than one dose times given as a vector.
#' @param dose.amount The doses given for the output.  Linear proportionality
#'   will be used from the input to output if they are not equal.  The length of
#'   dose.amount must be either 1 or matching the length of `dose.times`.
#' @param n.tau The number of tau dosing intervals to simulate or `Inf` for
#'   steady-state.
#' @param clast.pred To use predicted as opposed to observed Clast, either give
#'   the value for clast.pred here or set it to true (for automatic calculation
#'   from the half-life).
#' @param tlast The time of last observed concentration above the limit of
#'   quantification.  This is calculated if not provided.
#' @param additional.times Times to include in the final outputs in addition to
#'   the standard times (see details).  All `min(additional.times)` must be >=
#'   0, and the `max(additional.times)` must be <= `tau`.
#' @param check.blq Must the first concentration measurement be below the limit
#'   of quantification?
#' @param steady.state.tol The tolerance for assessing if steady-state has been
#'   achieved (between 0 and 1, exclusive).
#' @param ... Additional arguments passed to the `half.life` function if
#'   required to compute `lambda.z`.
#' @returns A data frame with columns named "conc" and "time".
#'
#' @details The returned superposition times will include all of the following
#' times: 0 (zero), `dose.times`, `time modulo tau` (shifting `time` for each
#' dose time as well), `additional.times`, and `tau`.
#'
#' @seealso [interp.extrap.conc()]
#' @export
superposition <- function(conc, ...) {
  UseMethod("superposition", conc)
}

#' @rdname superposition
#' @export
superposition.PKNCAconc <- function(conc, ...) {
  # Split the data by grouping and extract just the concentration and
  # time columns
  nested_data <- prepare_PKNCAconc(conc)
  tmp_results <-
    parallel::mclapply(
      X=seq_len(nrow(nested_data)),
      FUN=function(idx) {
        superposition.numeric(
          conc=nested_data$data_conc[[idx]]$conc,
          time=nested_data$data_conc[[idx]]$time,
          ...
        )
      }
    )
  # Replace the concentration data with the new results
  nested_data$data_conc <- tmp_results
  tidyr::unnest(nested_data, cols="data_conc")
}

#' @rdname superposition
#' @export
superposition.numeric <- function(conc, time, dose.input = NULL,
                                  tau, dose.times=0, dose.amount, n.tau=Inf,
                                  options=list(),
                                  lambda.z, clast.pred=FALSE, tlast,
                                  additional.times=numeric(),
                                  check.blq=TRUE,
                                  method = NULL,
                                  auc.type = "AUCinf",
                                  steady.state.tol=1e-3, ...) {
  # Check the inputs
  # Concentration and time
  assert_conc_time(conc = conc, time = time)
  if (check.blq) {
    if (!(conc[1] %in% 0)) {
      stop("The first concentration must be 0 (and not NA).  To change this set check.blq=FALSE.")
    }
  }
  assert_number_between(dose.input, na.ok = FALSE, null.ok = TRUE, lower = 0)
  assert_dosetau(tau)
  assert_numeric_between(x = dose.times, lower_eq = 0, min.len = 1, upper = tau)
  # dose.amount
  if (!missing(dose.amount)) {
    if (missing(dose.input)) {
      stop("must give dose.input to give dose.amount")
    }
    assert_numeric_between(x = dose.amount, lower = 0, finite = TRUE)
    if (!(length(dose.amount) %in% c(1, length(dose.times))))
      stop("dose.amount must either be a scalar or match the length of dose.times")
  }
  checkmate::assert_number(n.tau, lower = 1)
  if (is.finite(n.tau)) {
    n.tau <- checkmate::assert_integerish(n.tau, lower = 1)
  }
  # lambda.z
  if (!missing(lambda.z)) {
    lambda.z <- assert_number_between(x = lambda.z, lower = 0)
  }
  # clast.pred
  checkmate::assert(
    checkmate::check_number(clast.pred, lower = 0, na.ok = TRUE),
    checkmate::check_logical(clast.pred, any.missing = TRUE, len = 1),
    .var.name = "clast.pred"
  )
  if (is.na(clast.pred)) {
    clast.pred <- FALSE
  } else if (!is.numeric(clast.pred)) {
    checkmate::assert_logical(clast.pred, any.missing = FALSE, len = 1)
  }
  # tlast
  if (!missing(tlast)) {
    checkmate::assert_number(tlast)
  }
  # additional.times
  if (length(additional.times) > 0) {
    if (any(is.na(additional.times))) {
      stop("No additional.times may be NA (to not include any additional.times, enter c() as the function argument)")
    }
    if (!is.numeric(additional.times) | is.factor(additional.times))
      stop("additional.times must be a number")
    if (any(additional.times < 0))
      stop("All additional.times must be nonnegative")
    if (any(additional.times > tau))
      stop("All additional.times must be <= tau")
  }
  # steady.state.tol
  if (length(steady.state.tol) != 1)
    stop("steady.state.tol must be a scalar")
  if (!is.numeric(steady.state.tol) | is.factor(steady.state.tol) | is.na(steady.state.tol))
    stop("steady.state.tol must be a number")
  if (steady.state.tol <= 0 |
      steady.state.tol >= 1)
    stop("steady.state.tol must be between 0 and 1, exclusive.")
  if (steady.state.tol > 0.01)
    warning("steady.state.tol is usually <= 0.01")
  # We get all or none of lambda.z, clast, and tlast
  has.lambda.z <- !missing(lambda.z)
  has.clast.pred <- !is.logical(clast.pred)
  has.tlast <- !missing(tlast)
  if (any(c(has.lambda.z, has.clast.pred, has.tlast)) &
      !all(c(has.lambda.z, has.clast.pred, has.tlast)))
    stop("Either give all or none of the values for these arguments: lambda.z, clast.pred, and tlast")
  # combine dose.input and dose.amount as applicable to scale the
  # outputs.
  if (!missing(dose.amount)) {
    dose.scaling <- dose.amount / dose.input
    if (length(dose.scaling) != length(dose.times)) {
      if (length(dose.scaling) != 1)
        stop("bug in dose.amount, dose.times, and dose.input handling") # nocov
      # it is a scalar and there is more than one dose
      dose.scaling <- rep(dose.scaling, length(dose.times))
    }
  } else {
    dose.scaling <- rep(1, length(dose.times))
  }
  # Determine the output times
  tmp.times <- c(0, tau, dose.times, additional.times)
  # For the output, ensure that we have each of the input times
  # shifted for the dosing times.
  for (d in dose.times)
    tmp.times <- c(tmp.times, (d + time) %% tau)
  ret <- data.frame(conc=0, time=sort(unique(tmp.times)))
  # Check if all the input concentrations are 0, and if so, give that
  # simple answer.
  if (all(conc == 0)) {
    return(ret)
  } else {
    # Check if we will need to extrapolate
    tlast <- pk.calc.tlast(conc, time, check=FALSE)
    if ((tau*n.tau) > tlast) {
      if (missing(lambda.z)) {
        tmp <- pk.calc.half.life(conc, time, options=options,
                                 ..., check=FALSE)
        lambda.z <- tmp$lambda.z
      }
      if (identical(clast.pred, FALSE)) {
        clast <- pk.calc.clast.obs(conc, time, check=FALSE)
      } else if (identical(clast.pred, TRUE)) {
        clast <- tmp$clast.pred
      } else {
        # Use the provided clast.pred
        clast <- clast.pred
      }
    } else {
      # We don't need lambda.z for calculations, but it is simpler to
      # define it.
      lambda.z <- NA
    }
  }
  # cannot continue extrapolating due to missing data (likely due to
  # half-life not calculable)
  if ((n.tau * tau) > tlast & is.na(lambda.z)) {
    ret$conc <- NA
  } else {
    # Do the math! (Finally)
    current.tol <- steady.state.tol + 1
    tau.count <- 0
    # Stop either for reaching steady-state or for reaching the requested number of doses
    while (tau.count < n.tau &
           !is.na(current.tol) &
           current.tol >= steady.state.tol) {
      prev.conc <- ret$conc
      # Perform the dosing for a single dosing interval.
      for (i in seq_along(dose.times)) {
        tmp.time <- ret$time - dose.times[i] + (tau * tau.count)
        # For the first dosing interval, make sure that we do not
        # assign concentrations before the dose is given.
        mask.time <- tmp.time >= 0
        # Update the current concentration (previous concentration +
        # new concentration scaled by the relative dose)
        ret$conc[mask.time] <-
          (ret$conc[mask.time] +
           dose.scaling[i]*
           interp.extrap.conc(
             conc = conc, time = time,
             time.out=tmp.time[mask.time],
             lambda.z=lambda.z, clast=clast,
             options=options, check=FALSE,
             method = method, auc.type = auc.type
           )
          )
      }
      tau.count <- tau.count + 1
      if (any(ret$conc %in% 0)) {
        # prevent division by 0.  Since not all concentrations are 0,
        # all values will eventually be nonzero.
        current.tol <- steady.state.tol + 1
      } else {
        current.tol <- max(1-(prev.conc/ret$conc))
      }
    }
  }
  ret
}
