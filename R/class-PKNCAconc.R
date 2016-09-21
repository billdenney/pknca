#' Create a PKNCAconc object
#'
#' @param data A data frame with concentration, time, and the groups
#' defined in \code{formula}.
#' @param formula The formula defining the
#' \code{concentration~time|groups}
#' @param subject The column indicating the subject number (used for
#' plotting).  If not provided, this defaults to the beginning of the
#' inner groups: For example with
#' \code{concentration~time|Study+Subject/Analyte}, the inner groups
#' start with the first grouping variable before a \code{/},
#' \code{Subject}.  If there is only one grouping variable, it is
#' assumed to be the subject (e.g. \code{concentration~time|Subject}),
#' and if there are multiple grouping variables without a \code{/},
#' subject is assumed to be the last one.  For single-subject data, it
#' is assigned as \code{NULL}.
#' @param labels (optional) Labels for use when plotting.  They are a
#' named list where the names correspond to the names in the data
#' frame and the values are used for xlab and/or ylab as appropriate.
#' @param units (optional) Units for use when plotting and calculating
#' parameters.  Note that unit conversions and simplifications are not
#' done; the text is used as-is.
#' @param ... Ignored.
#' @return A PKNCAconc object that can be used for automated NCA.
#' @seealso \code{\link{PKNCAdata}}, \code{\link{PKNCAdose}}
#' @export
PKNCAconc <- function(data, ...)
  UseMethod("PKNCAconc")

#' @rdname PKNCAconc
#' @export
PKNCAconc.default <- function(data, ...)
  PKNCAconc.data.frame(as.data.frame(data), ...)
#' @rdname PKNCAconc
#' @export
PKNCAconc.tbl_df <- function(data, ...)
  PKNCAconc.data.frame(as.data.frame(data), ...)

#' @rdname PKNCAconc
#' @export
PKNCAconc.data.frame <- function(data, formula, subject, labels, units, ...) {
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
  ## Values must be unique (one value per measurement)
  key.cols <- c(all.vars(parsedForm$rhs),
                all.vars(parsedForm$groupFormula))
  if (any(mask.dup <- duplicated(data[,key.cols])))
    stop("Rows that are not unique per group and time (column names: ",
         paste(key.cols, collapse=", "),
         ") found within concentration data.  Row numbers: ",
         paste(seq_along(mask.dup)[mask.dup], collapse=", "))
  ## Assign the subject
  if (missing(subject)) {
    tmp.groups <- all.vars(parsedForm$groupFormula)
    if (length(tmp.groups) == 1) {
      subject <- tmp.groups
    } else {
      subject <- all.vars(findOperator(parsedForm$groupFormula,
                                       "/",
                                       side="left"))
      if (length(subject) == 0) {
        ## There is no / in the group formula, use the last element
        subject <- tmp.groups[length(tmp.groups)]
      } else if (length(subject) == 1) {
        ## There is a subject given; use it as is.
      } else {
        stop("Unknown how to handle subject definition from the fromula")
      }
    }
  } else {
    ## Ensure that the subject is part of the data definition and a
    ## scalar character string.
    if (!is.character(subject))
      stop("subject must be a character string")
    if (!(length(subject) == 1))
      stop("subject must be a scalar")
    if (!(subject %in% names(data)))
      stop("The subject parameter must map to a name in the data")
  }
  ret <- list(data=data,
              formula=formula,
              subject=subject)
  ## check and add labels and units
  if (!missing(labels))
    ret <- set.name.matching(ret, "labels", labels, data)
  if (!missing(units))
    ret <- set.name.matching(ret, "units", units, data)
  class(ret) <- c("PKNCAconc", class(ret))
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

#' @export
getDepVar.PKNCAconc <- function(x, ...) {
  x$data[, all.vars(parseFormula(x)$lhs)]
}

#' @export
getIndepVar.PKNCAconc <- function(x, ...) {
  x$data[, all.vars(parseFormula(x)$rhs)]
}

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
  grpnames <- all.vars(parseFormula(form)$groups)
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

#' Extract all the original data from a PKNCAconc or PKNCAdose object
#' @param object R object to extract the data from.
#' @export
getData.PKNCAconc <- function(object)
  object$data

#' Print and/or summarize a PKNCAconc or PKNCAdose object.
#'
#' @param x The object to print
#' @param object The object to summarize
#' @param n The number of rows of data to show (see
#' \code{\link{head}})
#' @param summarize Summarize the nested number of groups
#' @param \dots Arguments passed to \code{print.formula} and
#' \code{print.data.frame}
#' @export
print.PKNCAconc <- function(x, n=6, summarize=FALSE, ...) {
  cat(sprintf("Formula for concentration:\n "))
  print(stats::formula(x), ...)
  if (is.na(x$subject)) {
    cat("As a single-subject dataset.\n")
  } else {
    cat(sprintf("With %d subjects defined in the '%s' column.\n",
                length(unique(x$data[,x$subject])),
                x$subject))
  }
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
      cat("\nData for concentration:\n")
    } else if (n < 0) {
      cat(sprintf("\nFirst %d rows of concentration data:\n",
                  nrow(x$data)+n))
    } else {
      cat(sprintf("\nFirst %d rows of concentration data:\n",
                  n))
    }
    print.data.frame(utils::head(x$data, n=n), ..., row.names=FALSE)
  }
}

#' @rdname print.PKNCAconc
#' @export
summary.PKNCAconc <- function(object, n=0, summarize=TRUE, ...)
  print.PKNCAconc(object, n=n, summarize=summarize)

#' Plot a PKNCAconc object
#'
#' @param x The object to plot
#' @param groups The grouping variable for the plot (typically the
#' subject column)
#' @param \dots Additional arguments passed to \code{xyplot}
#' @param panel.formula The formula used for the call to xyplot
#' (defaults to the group formula of \code{x})
#' @param panel.formula.update Updates to the \code{panel.formula} to
#' simplify modifications without having to fully specify the formula.
#' @return A trellis object of the plot(s)
#' @export
plot.PKNCAconc <- function(x, ...,
                           groups=x$subject,
                           panel.formula=parseFormula(x)$groupFormula,
                           panel.formula.update) {
  ## Arguments that are set by default only overwrite if the user
  ## doesn't specify others with the same name.
  default.args <- list(lty=1, type="b", auto.key=list(columns=2),
                       scales=list(alternating=FALSE))
  call.args <- list(...)
  set.defaults <- setdiff(names(default.args), names(call.args))
  call.args[set.defaults] <- default.args[set.defaults]
  ## Update the panel formula if applicable
  if (!missing(panel.formula.update)) {
    panel.formula <- stats::update(panel.formula, panel.formula.update)
  }
  conc.formula <- parseFormula(x)
  ## If the groups are given, make sure that they are not in the
  ## panel.formula.
  if (!is.null(groups)) {
    conc.formula$groupFormula <-
      stats::update(panel.formula,
                    stats::as.formula(sprintf(".~-%s", groups)))
  }
  call.args[["x"]] <- stats::formula(conc.formula)
  call.args[["data"]] <- x$data
  ## If labels and/or units are given for the x and y variables, use
  ## them.
  xlab <- make.label("rhs", x$data, conc.formula, x$labels, x$units)
  ylab <- make.label("lhs", x$data, conc.formula, x$labels, x$units)
  if (!("xlab" %in% names(call.args)))
    call.args[["xlab"]] <- xlab
  if (!("ylab" %in% names(call.args)))
    call.args[["ylab"]] <- ylab
  do.call(lattice::xyplot, call.args)
}

