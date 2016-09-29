#' Create a PKNCAdose object
#' 
#' @param data A data frame with time and the groups defined in 
#'   \code{formula}.
#' @param formula The formula defining the 
#'   \code{dose.amount~time|groups} where \code{time} is the time of the
#'   dosing and \code{dose.amount} is the amount administered at that 
#'   time (see Details).
#' @param labels (optional) Labels for use when plotting.  They are a 
#'   named list where the names correspond to the names in the data 
#'   frame and the values are used for xlab and/or ylab as appropriate.
#' @param units (optional) Units for use when plotting and calculating 
#'   parameters.  Note that unit conversions and simplifications are not
#'   done; the text is used as-is.
#' @param ... Ignored.
#' @return A PKNCAconc object that can be used for automated NCA.
#' @details The \code{formula} for a \code{PKNCAdose} object can be 
#'   given three ways: one-sided (missing left side), one-sided (missing
#'   right side), or two-sided.  Each of the three ways can be given 
#'   with or without groups.  When given one-sided missing the left 
#'   side, the left side can either be omitted or can be given as a 
#'   period (\code{.}): \code{~time|treatment+subject} and 
#'   \code{.~time|treatment+subject} are identical, and dose-related NCA
#'   parameters will all be reported as not calculable (for example, 
#'   clearance).  When given one-sided missing the right side, the right
#'   side must be specified as a period (\code{.}): 
#'   \code{dose~.|treatment+subject}, and only a single row may be given
#'   per group.  When the right side is missing, PKNCA assumes that the
#'   same dose is given in every interval.  When given as a two-sided formula
#' @seealso \code{\link{PKNCAconc}}, \code{\link{PKNCAdata}}
#' @export
PKNCAdose <- function(data, ...)
  UseMethod("PKNCAdose")

#' @rdname PKNCAdose
#' @export
PKNCAdose.default <- function(data, ...)
  PKNCAdose.data.frame(as.data.frame(data), ...)
#' @rdname PKNCAdose
#' @export
PKNCAdose.tbl_df <- function(data, ...)
  PKNCAdose.data.frame(as.data.frame(data), ...)

#' @rdname PKNCAdose
#' @export
PKNCAdose.data.frame <- function(data, formula, labels, units, ...) {
  ## Verify that all the variables in the formula are columns in the
  ## data.
  parsedForm <- parseFormula(formula, require.two.sided=FALSE)
  ## Check for variable existence and length
  if (!(length(all.vars(parsedForm$lhs)) %in% c(0, 1)))
    stop("The left side of the formula must have zero or one variable")
  if (!(identical(parsedForm$lhs, NA) ||
        all.vars(parsedForm$lhs) %in% c(".", names(data)))) {
    stop("The left side formula must be a variable in the data, empty, or '.'.")
  }
  if (length(all.vars(parsedForm$rhs)) != 1)
    stop("The right side of the formula (excluding groups) must have exactly one variable")
  if (!(all.vars(parsedForm$rhs) %in% c(".", names(data)))) {
    stop("The right side formula must be a variable in the data or '.'.")
  }
  if (!all(all.vars(parsedForm$groups) %in% names(data))) {
    stop("All of the variables in the groups must be in the data")
  }
  ## Values must be unique (one value per measurement)
  key.cols <- c(setdiff(all.vars(parsedForm$rhs), "."),
                all.vars(parsedForm$groupFormula))
  if (any(mask.dup <- duplicated(data[,key.cols])))
    stop("Rows that are not unique per group and time (column names: ",
         paste(key.cols, collapse=", "),
         ") found within dosing data.  Row numbers: ",
         paste(seq_along(mask.dup)[mask.dup], collapse=", "))
  ret <- list(data=data,
              formula=formula)
  mask.indep <- is.na(getIndepVar.PKNCAdose(ret))
  if (any(mask.indep) & !all(mask.indep)) {
    stop("Some but not all values are missing for the independent variable, please see the help for PKNCAdose for how to specify the formula and confirm that your data has dose times for all doses.")
  }
  ## check and add labels and units
  if (!missing(labels))
    ret <- set.name.matching(ret, "labels", labels, data)
  if (!missing(units))
    ret <- set.name.matching(ret, "units", units, data)
  class(ret) <- c("PKNCAdose", class(ret))
  ret
}

#' @rdname formula.PKNCAconc
#' @export
formula.PKNCAdose <-  function(x, ...) {
  x$formula
}

#' @rdname model.frame.PKNCAconc
#' @export
model.frame.PKNCAdose <- function(formula, ...) {
  cbind(getDepVar.PKNCAdose(formula),
        getIndepVar.PKNCAdose(formula),
        getGroups.PKNCAdose(formula))
}

#' @export
getDepVar.PKNCAdose <- function(x, ...) {
  parsedForm <- parseFormula(x$formula, require.two.sided=FALSE)
  if (identical(parsedForm$lhs, NA) ||
      identical(all.vars(parsedForm$lhs), ".")) {
    rep(NA_integer_, nrow(x$data))
  } else {
    x$data[, all.vars(parseFormula(x)$lhs)]
  }
}

#' @export
getIndepVar.PKNCAdose <- function(x, ...) {
  parsedForm <- parseFormula(x$formula, require.two.sided=FALSE)
  if (identical(parsedForm$rhs, NA) ||
      identical(all.vars(parsedForm$rhs), ".")) {
    rep(NA_integer_, nrow(x$data))
  } else {
    x$data[, all.vars(parseFormula(x)$rhs)]
  }
}

#' @rdname getGroups.PKNCAconc
#' @export
getGroups.PKNCAdose <- function(...) {
  getGroups.PKNCAconc(...)
}

#' @rdname getData.PKNCAconc
#' @export
getData.PKNCAdose <-  function(object)
  object$data

#' @rdname print.PKNCAconc
#' @export
print.PKNCAdose <- function(x, n=6, summarize=FALSE, ...) {
  cat("Formula for dosing:\n ")
  print(stats::formula(x), ...)
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
      cat("\nData for dosing:\n")
    } else if (n < 0) {
      cat(sprintf("\nFirst %d rows of dosing data:\n",
                  nrow(x$data)+n))
    } else {
      cat(sprintf("\nFirst %d rows of dosing data:\n",
                  n))
    }
    print.data.frame(utils::head(x$data, n=n), ..., row.names=FALSE)
  }
}

#' @rdname print.PKNCAconc
#' @export
summary.PKNCAdose <- summary.PKNCAconc

#' @rdname plot.PKNCAconc
#' @export
plot.PKNCAdata <- function(x, ...)
  graphics::plot(x$conc, ...)
