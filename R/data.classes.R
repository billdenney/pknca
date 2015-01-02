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
  x[, parseFormula(formula(x))$lhs]
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
  UseMethod("getDepVar", x)

#' @export
getIndepVar.PKNCAconc <- function(x, ...) {
  x[, parseFormula(formula(x))$rhs]
}

#' @export
getIndepVar.PKNCAdose <- getDepVar.PKNCAconc

#' See \code{getGroups} in the nlme package.
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
  cat("Formula:\n ")
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
      cat("\nData:\n")
    } else if (n < 0) {
      cat(sprintf("\nFirst %d rows of data:\n", nrow(x$data)+n))
    } else {
      cat(sprintf("\nFirst %d rows of data:\n", n))
    }
    print.data.frame(head(x$data, n=n), ..., row.names=FALSE)
  }
}

#' @rdname print.PKNCAconc
#' @export
summary.PKNCAconc <- function(x, n=0, summarize=TRUE)
  print.PKNCAconc(x, n=n, summarize=summarize)

#' @rdname print.PKNCAconc
#' @export
print.PKNCAdose <- print.PKNCAconc

#' @rdname print.PKNCAconc
#' @export
summary.PKNCAdose <- summary.PKNCAconc

#' See \code{getData} in the nlme package.
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
#' @param single.dose.aucs A data frame of AUCs to use for single dose
#' data (see \code{\link{choose.acu.intervals}}).
#' @return A data frame with a row for each interval for each grouping
#' in the concentration data.
#' @seealso \code{\link{PKNCAconc}}, \code{\link{PKNCAdose}},
#' \code{\link{choose.auc.intervals}}
#' @export
PKNCAdata <- function(data.conc, formula.conc,
                      data.dose, formula.dose,
                      intervals,
                      single.dose.aucs=PKNCA.options("single.dose.aucs")) {
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
  ## Check the AUC intervals
  if (missing(intervals)) {
    ret$intervals <-
      choose.auc.intervals(getIndepVar(ret$conc),
                           getIndepVar(ret$dose),
                           single.dose.aucs=single.dose.aucs)
  } else {
    ret$intervals <- check.auc.specification(intervals)
  }
  class(ret) <- c("PKNCAdata", class(ret))
  ret
}

#' @export
print.PKNCAdata <- function(x, ...) {
  print.PKNCAconc(x$conc, ...)
  print.PKNCAdose(x$dose, ...)
  cat(sprintf("\nWith %d rows of AUC specifications.\n",
              nrow(x$intervals)))
}

#' @export
summary.PKNCAdata <- function(x, ...)
  print.PKNCAdata(x, summarize=TRUE, ...)
