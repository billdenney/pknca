#' Create a PKNCAconc object
#'
#' @param data A data frame with concentration, time, and the groups
#' defined in \code{formula}.
#' @param formula The formula defining the
#' \code{concentration~time|groups}
#' @return A PKNCAconc object that can be used for automated NCA.
#' @export
PKNCAconc <- function(data, formula) {
  ## Verify that all the variables in the formula are columns in the
  ## data.
  if (!all(all.vars(formula) %in% names(data))) {
    stop("All of the variables in the formula must be in the data")
  }
  parsedForm <- parseFormula(formula, require.two.sided=TRUE)
  if (length(all.vars(parsedForm$lhs)) != 1)
    stop("The left hand side of the formula must have exactly one variable")
  if (length(all.vars(parsedForm$rhs)) != 1)
    stop("The right hand side of the formula (excluding groups) must have exactly one variable")
  ret <- list(data=data,
              formula=formula)
  class(ret) <- c("PKNCAconc", class(ret))
  ret
}

#' Create a PKNCAdose object
#'
#' @param data A data frame with time and the groups
#' defined in \code{formula}.
#' @param formula The formula defining the \code{~time|groups} where
#' \code{time} is the time of the dosing.
#' @return A PKNCAconc object that can be used for automated NCA.
#' @export
PKNCAdose <- function(data, formula) {
  ## Verify that all the variables in the formula are columns in the
  ## data.
  if (!all(all.vars(formula) %in% names(data))) {
    stop("All of the variables in the formula must be in the data")
  }
  parsedForm <- parseFormula(formula, require.two.sided=FALSE)
  if (length(all.vars(parsedForm$lhs)) != 0)
    stop("The formula must be one-sided")
  if (length(all.vars(parsedForm$rhs)) != 1)
    stop("The right hand side of the formula (excluding groups) must have exactly one variable")
  ret <- list(data=data,
              formula=formula)
  class(ret) <- c("PKNCAdose", class(ret))
  ret
}

#' Extract the formula from a PKNCAconc object.
#'
#' @param x The object to extract the formula from.
#' @param \dots Unused
#' @return A formula object
#' @export
formula.PKNCAconc <- function(x, ...)
  x$formula

#' @rdname formula.PKNCAconc
#' @export
formula.PKNCAdose <- formula.PKNCAconc

#' Get the dependent variable (left hand side of the formula) from a
#' PKNCA object.
#'
#' @param x The object to extract the formula from
#' @param \dots Unused
#' @return The vector of the dependent variable from the object.
#' @export
getDepVar <- function(x, ...)
  UseMethod("getDepVar", x)

#' @export
getDepVar.PKNCAconc <- function(x, ...) {
  x$data[, all.vars(parseFormula(x)$lhs)]
}

#' @export
getDepVar.PKNCAdose <- getDepVar.PKNCAconc

#' Get the independent variable (right hand side of the formula) from
#' a PKNCA object.
#'
#' @param x The object to extract the formula from
#' @param \dots Unused
#' @return The vector of the independent variable from the object.
#' @export
getIndepVar <- function(x, ...)
  UseMethod("getIndepVar", x)

#' @export
getIndepVar.PKNCAconc <- function(x, ...) {
  x$data[, all.vars(parseFormula(x)$rhs)]
}

#' @export
getIndepVar.PKNCAdose <- getIndepVar.PKNCAconc

## See \code{getGroups} in the nlme package.
getGroups <- function(object, form, level, data, sep)
  UseMethod("getGroups")

#' Get the groups (right hand side after the \code{|} from a PKNCA
#' object.
#'
#' @param object The object to extract the data from
#' @param form The formula to extract the data from (defaults to the
#' formula from \code{object})
#' @param level optional.  If included, this specifies the level(s) of
#' the groups to include.  If a numeric scalar, include the first
#' \code{level} number of groups.  If a numeric vector, include each
#' of the groups specified by the number.  If a character vector,
#' include the named group levels.
#' @param data The data to extract the groups from (defaults to the
#' data from \code{object})
#' @param sep Unused (kept for compatibility with the nlme package)
#' @return A data frame with the (selected) group columns.
#' @export
getGroups.PKNCAconc <- function(object, form=formula(object), level,
                                data=getData(object), sep) {
  grpnames <- all.vars(parseFormula(formula(object))$groups)
  if (!missing(level))
    if (is.factor(level) | is.character(level)) {
      level <- as.character(level)
      if (any(!(level %in% grpnames)))
        stop("Not all levels are listed in the group names.  Missing levels are: ",
             paste(setdiff(level, grpnames), collapse=", "))
      grpnames <- level
    } else if (is.numeric(level)) {
      if (length(level) == 1) {
        grpnames <- grpnames[1:level]
      } else {
        grpnames <- grpnames[level]
      }
    }
  data[, grpnames, drop=FALSE]
}

#' @rdname getGroups.PKNCAconc
#' @export
getGroups.PKNCAdose <- getGroups.PKNCAconc

#' Print and/or summarize a PKNCAconc or PKNCAdose object.
#'
#' @param x The object to print
#' @param n The number of rows of data to show (see
#' \code{\link{head}})
#' @param summarize Summarize the nested number of groups
#' @param \dots Arguments passed to \code{print.formula} and
#' \code{print.data.frame}
#' @export
print.PKNCAconc <- function(x, n=6, summarize=FALSE, ...) {
  if (inherits(x, "PKNCAconc")) {
    data.type <- "concentration"
  } else if (inherits(x, "PKNCAdose")) {
    data.type <- "dosing"
  } else {
    stop("Unknown data type")
  }
  cat(sprintf("Formula for %s:\n ", data.type))
  print(formula(x), ...)
  if (summarize) {
    cat("\n")
    grp <- getGroups(x)
    if (ncol(grp) > 0) {
      tmp.summary <- data.frame(Group.Name=names(grp),
                                Count=0)
      for (i in 1:ncol(grp))
        tmp.summary$Count[i] <- nrow(unique(grp[,1:i,drop=FALSE]))
      cat("Group summary:\n")
      names(tmp.summary) <- gsub("\\.", " ", names(tmp.summary))
      print.data.frame(tmp.summary, row.names=FALSE)
    } else {
      cat("No groups.\n")
    }
  }
  if (n != 0) {
    if (n >= nrow(x$data)) {
      cat(sprintf("\nData for %s:\n", data.type))
    } else if (n < 0) {
      cat(sprintf("\nFirst %d rows of %s data:\n",
                  nrow(x$data)+n, data.type))
    } else {
      cat(sprintf("\nFirst %d rows of %s data:\n",
                  n, data.type))
    }
    print.data.frame(head(x$data, n=n), ..., row.names=FALSE)
  }
}

#' @rdname print.PKNCAconc
#' @export
summary.PKNCAconc <- function(object, n=0, summarize=TRUE, ...)
  print.PKNCAconc(object, n=n, summarize=summarize)

#' @rdname print.PKNCAconc
#' @export
print.PKNCAdose <- print.PKNCAconc

#' @rdname print.PKNCAconc
#' @export
summary.PKNCAdose <- summary.PKNCAconc

#' See \code{getData} in the nlme package.
#'
#' @param object Object to extract data from
#' @return Data from the object
getData <- function(object)
  UseMethod("getData")

#' Extract all the original data from a PKNCAconc or PKNCAdose object
#'
#' @param object The object to extract data from
#' @return The data from the object
#' @export
getData.PKNCAconc <- function(object)
  object$data

#' @rdname getData.PKNCAconc
#' @export
getData.PKNCAdose <- getData.PKNCAconc

#' Extract the columns used in the formula (in order) from a PKNCAconc
#' or PKNCAdose object.
#'
#' @param formula The object to use (parameter name is \code{formula}
#' to use the generic function)
#' @param \dots Unused
#' @return A data frame with the columns from the object in formula
#' order.
#' @export
model.frame.PKNCAconc <- function(formula, ...)
  formula$data[, all.vars(formula$formula), drop=FALSE]

#' @rdname model.frame.PKNCAconc
#' @export
model.frame.PKNCAdose <- model.frame.PKNCAconc

#' Create a PKNCAdata object.
#'
#' \code{PKNCAdata} combines \code{PKNCAconc} and \code{PKNCAdose} and
#' adds in the intervals for PK calculations.
#' 
#' @param data.conc Concentration data as a \code{PKNCAconc} object or
#' a data frame
#' @param formula.conc Formula for making a \code{PKNCAconc} object
#' with \code{data.conc}.  This must be given if \code{data.conc} is a
#' data.frame, and it must not be given if \code{data.conc} is a
#' \code{PKNCAconc} object.
#' @param data.dose Dosing data as a \code{PKNCAdose} object
#' @param formula.dose Formula for making a \code{PKNCAdose} object
#' with \code{data.dose}.  This must be given if \code{data.dose} is a
#' data.frame, and it must not be given if \code{data.dose} is a
#' \code{PKNCAdose} object.
#' @param intervals A data frame with the AUC interval specifications
#' as defined in \code{\link{check.auc.specification}}.  If missing,
#' this will be automatically chosen by
#' \code{\link{choose.auc.intervals}}.
#' @param options List of changes to the default
#' \code{\link{PKNCA.options}} for calculations.
#' @return A PKNCAdata object with concentration, dose, interval, and
#' calculation options stored (note that PKNCAdata objects can also
#' have results after a NCA calculations are done to the data).
#' @seealso \code{\link{PKNCAconc}}, \code{\link{PKNCAdose}},
#' \code{\link{choose.auc.intervals}}
#' @export
PKNCAdata <- function(data.conc, formula.conc,
                      data.dose, formula.dose,
                      intervals, options=list()) {
  ret <- list()
  ## Generate the conc element
  if (inherits(data.conc, "PKNCAconc")) {
    if (!missing(formula.conc))
      warning("data.conc was given as a PKNCAconc object.  Ignoring formula.conc")
    ret$conc <- data.conc
  } else {
    ret$conc <- PKNCAconc(data.conc, formula.conc)
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
      merge(conc=splitBy(parseFormula(ret$conc)$groupFormula, ret$conc$data),
            dose=splitBy(parseFormula(ret$dose)$groupFormula, ret$dose$data))
    groupid <- attributes(tmp.conc.dose)$groupid
    rownames(groupid) <- NULL
    intervals <- data.frame()
    indep.var.conc <- all.vars(parseFormula(ret$conc)$rhs)
    indep.var.dose <- all.vars(parseFormula(ret$dose)$rhs)
    for (i in 1:nrow(groupid)) {
      tmp.group <- groupid[i,,drop=FALSE]
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
    }
  }
  ret$intervals <- check.auc.specification(intervals)
  ## Assign the class and give it all back to the user.
  class(ret) <- c("PKNCAdata", class(ret))
  ret
}

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

#' @export
summary.PKNCAdata <- function(object, ...)
  print.PKNCAdata(object, summarize=TRUE, ...)

#' Generate a PKNCAresults object
#'
#' This function should not be run directly.  The object is created
#' for summarization and plotting.
#'
#' @param result a data frame with NCA calculation results and groups.
#' Each row is one interval and each column is a group name or the
#' name of an NCA parameter.
#' @param formula The formula used for concentration data in the
#' calculations.  The groups are verified to be column names in the
#' \code{result} parameter.
#' @param options Options that are different from the defaults.  All
#' options (default and custom) are stored with the results.
#' @return A PKNCA object with each of the above within.
#' @export
PKNCAresults <- function(result, formula, options=list()) {
  ## Merge the options into the default options.
  tmp.opt <- PKNCA.options()
  tmp.opt[names(options)] <- options
  ## Add all the parts into the object
  ret <- list(result=results,
              formula=formula,
              options=tmp.opt)
  class(ret) <- c("PKNCAresults", class(ret))
  ret
}
