#' Create a PKNCAdose object
#'
#' @param data A data frame with time and the groups defined in
#'   \code{formula}.
#' @param formula The formula defining the
#'   \code{dose.amount~time|groups} where \code{time} is the time of
#'   the dosing and \code{dose.amount} is the amount administered at
#'   that time.
#' @param labels (optional) Labels for use when plotting.  They are a
#'   named list where the names correspond to the names in the data
#'   frame and the values are used for xlab and/or ylab as
#'   appropriate.
#' @param units (optional) Units for use when plotting and calculating
#'   parameters.  Note that unit conversions and simplifications are
#'   not done; the text is used as-is.
#' @param ... Ignored.
#' @return A PKNCAconc object that can be used for automated NCA.
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
  if (!all(all.vars(formula) %in% names(data))) {
    stop("All of the variables in the formula must be in the data")
  }
  parsedForm <- parseFormula(formula, require.two.sided=FALSE)
  if (!(length(all.vars(parsedForm$lhs)) %in% c(0, 1)))
    stop("The left hand side of the formula must have zero or one variable")
  if (length(all.vars(parsedForm$rhs)) != 1)
    stop("The right hand side of the formula (excluding groups) must have exactly one variable")
  ## Values must be unique (one value per measurement)
  key.cols <- c(all.vars(parsedForm$rhs),
                all.vars(parsedForm$groupFormula))
  if (any(mask.dup <- duplicated(data[,key.cols])))
    stop("Rows that are not unique per group and time (column names: ",
         paste(key.cols, collapse=", "),
         ") found within dosing data.  Row numbers: ",
         paste(seq_along(mask.dup)[mask.dup], collapse=", "))
  ret <- list(data=data,
              formula=formula)
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
formula.PKNCAdose <-  function(x, ...)
  x$formula

#' @rdname model.frame.PKNCAconc
#' @export
model.frame.PKNCAdose <- function(formula, ...)
  model.frame.PKNCAconc(formula, ...)

#' @export
getDepVar.PKNCAdose <- function(x, ...) {
  x$data[, all.vars(parseFormula(x)$lhs)]
}

#' @export
getIndepVar.PKNCAdose <- function(x, ...) {
  x$data[, all.vars(parseFormula(x)$rhs)]
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
