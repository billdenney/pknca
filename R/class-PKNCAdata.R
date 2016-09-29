#' Create a PKNCAdata object.
#'
#' \code{PKNCAdata} combines \code{PKNCAconc} and \code{PKNCAdose} and
#' adds in the intervals for PK calculations.
#' 
#' @param data.conc Concentration data as a \code{PKNCAconc} object or
#'   a data frame
#' @param data.dose Dosing data as a \code{PKNCAdose} object
#' @param formula.conc Formula for making a \code{PKNCAconc} object
#'   with \code{data.conc}.  This must be given if \code{data.conc} is
#'   a data.frame, and it must not be given if \code{data.conc} is a
#'   \code{PKNCAconc} object.
#' @param formula.dose Formula for making a \code{PKNCAdose} object
#'   with \code{data.dose}.  This must be given if \code{data.dose} is
#'   a data.frame, and it must not be given if \code{data.dose} is a
#'   \code{PKNCAdose} object.
#' @param intervals A data frame with the AUC interval specifications
#'   as defined in \code{\link{check.interval.specification}}.  If
#'   missing, this will be automatically chosen by
#'   \code{\link{choose.auc.intervals}}.
#' @param options List of changes to the default
#'   \code{\link{PKNCA.options}} for calculations.
#' @param ... arguments passed to \code{PKNCAdata.default}
#' @return A PKNCAdata object with concentration, dose, interval, and
#'   calculation options stored (note that PKNCAdata objects can also
#'   have results after a NCA calculations are done to the data).
#' @seealso \code{\link{PKNCAconc}}, \code{\link{PKNCAdose}},
#'   \code{\link{choose.auc.intervals}}, \code{\link{pk.nca}}
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
PKNCAdata.PKNCAdose <- function(data.conc, data.dose, ...)
  ## Swap the arguments
  PKNCAdata.default(data.dose=data.conc, data.conc=data.dose, ...)

#' @rdname PKNCAdata
#' @export
PKNCAdata.default <- function(data.conc, data.dose, ...,
                              formula.conc, formula.dose,
                              intervals, options=list()) {
  ret <- list()
  ## Generate the conc element
  if (inherits(data.conc, "PKNCAconc")) {
    if (!missing(formula.conc))
      warning("data.conc was given as a PKNCAconc object.  Ignoring formula.conc")
    ret$conc <- data.conc
  } else {
    ret$conc <- PKNCAconc(data.conc, formula=formula.conc)
  }
  ## Generate the data element
  if (inherits(data.dose, "PKNCAdose")) {
    if (!missing(formula.dose))
      warning("data.dose was given as a PKNCAdose object.  Ignoring formula.dose")
    ret$dose <- data.dose
  } else {
    ret$dose <- PKNCAdose(data.dose, formula.dose)
  }
  ## Check the options
  if (!is.list(options))
    stop("options must be a list.")
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
  if (missing(intervals)) {
    ## Generate the intervals for each grouping of concentration and
    ## dosing.
    tmp.conc.dose <-
      merge(conc=doBy::splitBy(parseFormula(ret$conc)$groupFormula,
                               ret$conc$data),
            dose=doBy::splitBy(parseFormula(ret$dose)$groupFormula,
                               ret$dose$data))
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
        new.intervals <-
          cbind(
            tmp.group,
            choose.auc.intervals(tmp.conc.dose[[i]]$conc[,indep.var.conc],
                                 tmp.conc.dose[[i]]$dose[,indep.var.dose],
                                 single.dose.aucs=PKNCA.choose.option("single.dose.aucs",
                                                                      options)))
        intervals <-
          rbind(intervals, new.intervals)
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
  print.PKNCAdose(x$dose, ...)
  cat(sprintf("\nWith %d rows of AUC specifications.\n",
              nrow(x$intervals)))
  if (length(x$options) == 0) {
    cat("No options are set differently than default.\n")
  } else {
    cat("Options changed from default are:\n")
    print(x$options)
  }
}

#' Summarize a PKNCAdata object showing important details about the
#' concentration, dosing, and interval information.
#' @param object The PKNCAdata object to summarize.
#' @param ... arguments passed on to \code{\link{print.PKNCAdata}}
#' @export
summary.PKNCAdata <- function(object, ...)
  print.PKNCAdata(object, summarize=TRUE, ...)

