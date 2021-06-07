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

## Ensure that arguments are reversible
#' @rdname PKNCAdata
#' @export
PKNCAdata.PKNCAconc <- function(data.conc, data.dose, ...)
  PKNCAdata.default(data.conc=data.conc, data.dose=data.dose, ...)

#' @rdname PKNCAdata
#' @export
PKNCAdata.PKNCAdose <- function(data.conc, data.dose, ...) {
  ## Swap the arguments
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
  ## Generate the conc element
  if (inherits(data.conc, "PKNCAconc")) {
    if (!missing(formula.conc)) {
      warning("data.conc was given as a PKNCAconc object.  Ignoring formula.conc")
    }
    ret$conc <- data.conc
  } else {
    ret$conc <- PKNCAconc(data.conc, formula=formula.conc)
  }
  ## Generate the dose element
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
  ## Check the options
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
  ## Check the AUC intervals
  if (missing(intervals) & identical(ret$dose, NA)) {
    stop("If data.dose is not given, intervals must be given")
  } else if (missing(intervals)) {
    ## Generate the intervals for each grouping of concentration and
    ## dosing.
    tmp.conc.dose <-
      merge.splitlist(
        conc=split(ret$conc),
        dose=split(ret$dose)
      )
    groupid <- attributes(tmp.conc.dose)$groupid
    rownames(groupid) <- NULL
    intervals <- data.frame()
    indep.var.conc <- all.vars(parseFormula(ret$conc)$rhs)
    indep.var.dose <- all.vars(parseFormula(ret$dose)$rhs)
    if (identical(indep.var.dose, ".")) {
      stop("Dose times were not given, so intervals must be manually specified.")
    }
    for (i in seq_len(nrow(groupid))) {
      tmp.group <- groupid[i,,drop=FALSE]
      if (!is.null(tmp.conc.dose[[i]]$conc)) {
        rownames(tmp.group) <- NULL
        generated_intervals <-
          choose.auc.intervals(
            tmp.conc.dose[[i]]$conc$data[,indep.var.conc],
            tmp.conc.dose[[i]]$dose$data[,indep.var.dose],
            options=options
          )
        if (nrow(generated_intervals)) {
          new.intervals <- cbind(tmp.group, generated_intervals)
          intervals <- rbind(intervals, new.intervals)
        } else {
          warning("No intervals generated likely due to limited concentration data for ",
                  paste(names(tmp.group),
                        unlist(lapply(tmp.group, as.character)),
                        sep="=", collapse=", "))
        }
      } else {
        warning("No intervals generated due to no concentration data for ",
                paste(names(tmp.group),
                      unlist(lapply(tmp.group, as.character)),
                      sep="=", collapse=", "))
      }
    }
  }
  ret$intervals <- check.interval.specification(intervals)
  ## Assign the class and give it all back to the user.
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
getData.PKNCAdata <- function(object)
  object$data

#' @rdname getDataName
getDataName.PKNCAdata <- function(object)
  "data"

#' Summarize a PKNCAdata object showing important details about the
#' concentration, dosing, and interval information.
#' @param object The PKNCAdata object to summarize.
#' @param ... arguments passed on to \code{\link{print.PKNCAdata}}
#' @export
summary.PKNCAdata <- function(object, ...)
  print.PKNCAdata(object, summarize=TRUE, ...)

#' @rdname split.PKNCAconc
#' @export
split.PKNCAdata <- function(x, ...) {
  interval.group.cols <- intersect(names(x$intervals),
                                   all.vars(parseFormula(x$conc$formula)$groups))
  if (length(interval.group.cols) > 0) {
    # If the intervals need to be split across the groups
    tmp.interval.split <-
      split.PKNCAconc(list(data=x$intervals),
                      f=x$intervals[, interval.group.cols, drop=FALSE])
    tmp.attr <- attributes(tmp.interval.split)
    tmp.interval.split <- lapply(tmp.interval.split, function(x) x$data)
    attributes(tmp.interval.split) <- tmp.attr
    if (identical(NA, x$dose)) {
      ret <-
        merge.splitlist(conc=split.PKNCAconc(x$conc),
                        intervals=tmp.interval.split)
      ret <- lapply(X=ret,
                    FUN=function(x) {
                      x$dose <- NA
                      x
                    })
    } else {
      ret <-
        merge.splitlist(conc=split.PKNCAconc(x$conc),
                        dose=split.PKNCAdose(x$dose),
                        intervals=tmp.interval.split)
    }
  } else {
    # If the intervals apply to all groups
    if (identical(NA, x$dose)) {
      ret <- lapply(X=split.PKNCAconc(x$conc),
                    FUN=function(x) {
                      list(conc=x,
                           dose=NA)
                    })
    } else {
      ret <-
        merge.splitlist(conc=split.PKNCAconc(x$conc),
                        dose=split.PKNCAdose(x$dose))
    }
    ret <- lapply(X=ret,
                  FUN=function(x, intervals) {
                    x$intervals <- intervals
                    x
                  }, intervals=x$intervals)
  }
  for (n in setdiff(names(x), c("conc", "dose", "intervals"))) {
    # Add any other attributes to all the splits (like options)
    ret <-
      lapply(ret, function(x, name, value) {
        x[[name]] <- value
        x
      },
      name=n, value=x[[n]])
  }
  ret
}
