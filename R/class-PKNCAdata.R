#' Create a PKNCAdata object.
#' 
#' \code{PKNCAdata} combines \code{PKNCAconc} and \code{PKNCAdose} and 
#' adds in the intervals for PK calculations.
#' 
#' @param data.conc Concentration data as a \code{PKNCAconc} object or a
#'   data frame
#' @param data.dose Dosing data as a \code{PKNCAdose} object (see 
#'   details)
#' @param formula.conc Formula for making a \code{PKNCAconc} object with
#'   \code{data.conc}.  This must be given if \code{data.conc} is a 
#'   data.frame, and it must not be given if \code{data.conc} is a 
#'   \code{PKNCAconc} object.
#' @param formula.dose Formula for making a \code{PKNCAdose} object with
#'   \code{data.dose}.  This must be given if \code{data.dose} is a 
#'   data.frame, and it must not be given if \code{data.dose} is a 
#'   \code{PKNCAdose} object.
#' @param intervals A data frame with the AUC interval specifications as
#'   defined in \code{\link{check.interval.specification}}.  If missing,
#'   this will be automatically chosen by 
#'   \code{\link{choose.auc.intervals}}. (see details)
#' @param options List of changes to the default 
#'   \code{\link{PKNCA.options}} for calculations.
#' @param ... arguments passed to \code{PKNCAdata.default}
#' @return A PKNCAdata object with concentration, dose, interval, and 
#'   calculation options stored (note that PKNCAdata objects can also 
#'   have results after a NCA calculations are done to the data).
#' @details If \code{data.dose} is not given or is \code{NA}, then the 
#'   \code{intervals} must be given.  At least one of \code{data.dose}
#'   and \code{intervals} must be given.
#' @family PKNCA objects
#' @seealso \code{\link{choose.auc.intervals}}, \code{\link{pk.nca}}
#' @export
PKNCAdata <- function(data.conc, data.dose, ...)
  UseMethod("PKNCAdata", data.conc)

# Ensure that arguments are reversible
#' @rdname PKNCAdata
#' @export
PKNCAdata.PKNCAconc <- function(data.conc, data.dose, ...)
  PKNCAdata.default(data.conc=data.conc, data.dose=data.dose, ...)

#' @rdname PKNCAdata
#' @export
PKNCAdata.PKNCAdose <- function(data.conc, data.dose, ...) {
  # Swap the arguments
  PKNCAdata.default(data.dose=data.conc, data.conc=data.dose, ...)
}

#' @rdname PKNCAdata
#' @export
PKNCAdata.default <- function(data.conc, data.dose, ...,
                              formula.conc, formula.dose,
                              intervals, options=list()) {
  if (length(list(...))) {
    stop("Unknown argument provided to PKNCAdata.  All arguments other than `data.conc` and `data.dose` must be named.")
  }
  ret <- list()
  # Generate the conc element
  if (inherits(data.conc, "PKNCAconc")) {
    if (!missing(formula.conc)) {
      warning("data.conc was given as a PKNCAconc object.  Ignoring formula.conc")
    }
    ret$conc <- data.conc
  } else {
    ret$conc <- PKNCAconc(data.conc, formula=formula.conc)
  }
  # Generate the dose element
  if (missing(data.dose)) {
    ret$dose <- NA
  } else if (identical(data.dose, NA)) {
    ret$dose <- NA
  } else if (inherits(data.dose, "PKNCAdose")) {
    if (!missing(formula.dose))
      warning("data.dose was given as a PKNCAdose object.  Ignoring formula.dose")
    ret$dose <- data.dose
  } else {
    ret$dose <- PKNCAdose(data.dose, formula.dose)
  }
  # Check the options
  if (!is.list(options)) {
    stop("options must be a list.")
  }
  if (length(options) > 0) {
    if (is.null(names(options)))
      stop("options must have names.")
    for (n in names(options)) {
      tmp.opt <- list(options[[n]], TRUE)
      names(tmp.opt) <- c(n, "check")
      do.call(PKNCA.options, tmp.opt)
    }
  }
  ret$options <- options
  # Check the intervals
  if (missing(intervals) & identical(ret$dose, NA)) {
    stop("If data.dose is not given, intervals must be given")
  } else if (missing(intervals)) {
    # Generate the intervals for each grouping of concentration and
    # dosing.
    if (identical(all.vars(parseFormula(ret$dose)$rhs), ".")) {
      stop("Dose times were not given, so intervals must be manually specified.")
    }
    n_conc_dose <-
      full_join_PKNCAconc_PKNCAdose(
        conc=ret$conc,
        dose=ret$dose
      )
    n_conc_dose$data_intervals <- rep(list(NULL), nrow(n_conc_dose))
    for (idx in seq_len(nrow(n_conc_dose))) {
      current_conc <- n_conc_dose$data_conc[[idx]]
      current_dose <- n_conc_dose$data_dose[[idx]]
      current_group <-
        n_conc_dose[
          idx,
          setdiff(names(n_conc_dose), c("data_conc", "data_dose")),
          drop=FALSE
        ]
      warning_prefix <-
        if (ncol(current_group) > 0) {
          paste0(
            paste(names(current_group), unlist(lapply(current_group, as.character)), sep="=", collapse="; "),
            ": "
          )
        } else {
          ""
        }
      if (!is.null(current_conc)) {
        generated_intervals <-
          choose.auc.intervals(
            current_conc$time,
            current_dose$time,
            options=options
          )
        if (nrow(generated_intervals) > 0) {
          n_conc_dose$data_intervals[[idx]] <- generated_intervals
        } else {
          warning(warning_prefix, "No intervals generated likely due to limited concentration data")
        }
      } else {
        warning(warning_prefix, "No intervals generated due to no concentration data")
      }
    }
    intervals <-
      tidyr::unnest(
        n_conc_dose[, setdiff(names(n_conc_dose), c("data_conc", "data_dose")), drop=FALSE],
        cols="data_intervals"
      )
  }
  ret$intervals <- check.interval.specification(intervals)
  # Verify that either everything or nothing is using units
  units_interval_start <- inherits(ret$intervals$start, "units")
  units_interval_end <- inherits(ret$intervals$end, "units")
  # Assign the class and give it all back to the user.
  class(ret) <- c("PKNCAdata", class(ret))
  ret
}

#' Print a PKNCAdata object
#' @param x The object to print
#' @param ... Arguments passed on to \code{\link{print.PKNCAconc}} and
#' \code{\link{print.PKNCAdose}}
#' @export
print.PKNCAdata <- function(x, ...) {
  print.PKNCAconc(x$conc, ...)
  if (identical(NA, x$dose)) {
    cat("No dosing information.\n")
  } else {
    print.PKNCAdose(x$dose, ...)
  }
  cat(sprintf("\nWith %d rows of AUC specifications.\n",
              nrow(x$intervals)))
  if (length(x$options) == 0) {
    cat("No options are set differently than default.\n")
  } else {
    cat("Options changed from default are:\n")
    print(x$options)
  }
}

#' Extract all the original data from a PKNCAconc or PKNCAdose object
#' @param object R object to extract the data from.
#' @export
getData.PKNCAdata <- function(object) {
  object$data
}

#' @rdname getDataName
getDataName.PKNCAdata <- function(object) {
  "data"
}

#' Summarize a PKNCAdata object showing important details about the
#' concentration, dosing, and interval information.
#' @param object The PKNCAdata object to summarize.
#' @param ... arguments passed on to \code{\link{print.PKNCAdata}}
#' @export
summary.PKNCAdata <- function(object, ...) {
  print.PKNCAdata(object, summarize=TRUE, ...)
}
