#' Create a PKNCAconc object
#'
#' @param data A data frame with concentration (or amount for
#'   urine/feces), time, and the groups defined in \code{formula}.
#' @param formula The formula defining the
#'   \code{concentration~time|groups} or \code{amount~time|groups} for
#'   urine/feces (In the remainder of the documentation, "concentration" will be
#'   used to describe concentration or amount.)  One special aspect of the
#'   \code{groups} part of the formula is that the last group is typically
#'   assumed to be the \code{subject}; see the documentation for the
#'   \code{subject} argument for exceptions to this assumption.
#' @param subject The column indicating the subject number.  If not provided,
#'   this defaults to the beginning of the inner groups: For example with
#'   \code{concentration~time|Study+Subject/Analyte}, the inner groups start
#'   with the first grouping variable before a \code{/}, \code{Subject}.  If
#'   there is only one grouping variable, it is assumed to be the subject (e.g.
#'   \code{concentration~time|Subject}), and if there are multiple grouping
#'   variables without a \code{/}, subject is assumed to be the last one.  For
#'   single-subject data, it is assigned as \code{NULL}.
#' @param time.nominal (optional) The name of the nominal time column
#'   (if the main time variable is actual time.  The \code{time.nominal}
#'   is not used during calculations; it is available to assist with
#'   data summary and checking.
#' @param exclude (optional) The name of a column with concentrations to
#'   exclude from calculations and summarization.  If given, the column
#'   should have values of \code{NA} or \code{""} for concentrations to
#'   include and non-empty text for concentrations to exclude.
#' @param volume (optional) The volume (or mass) of collection as is
#'   typically used for urine or feces measurements.
#' @param duration (optional) The duration of collection as is typically
#'   used for concentration measurements in urine or feces.
#' @param exclude_half.life,include_half.life A character scalar for the column
#'   name in the dataset of the points to exclude from the half-life calculation
#'   (still using normal curve-stripping selection rules for the other points)
#'   or to include for the half-life (using specifically those points and
#'   bypassing automatic curve-stripping point selection).  See the "Half-Life
#'   Calculation" vignette for more details on the use of these arguments.
#' @param ... Ignored.
#' @return A PKNCAconc object that can be used for automated NCA.
#' @family PKNCA objects
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
PKNCAconc.data.frame <- function(data, formula, subject,
                                 time.nominal, exclude, duration, volume,
                                 exclude_half.life, include_half.life, ...) {
  ## The data must have... data
  if (nrow(data) == 0) {
    stop("data must have at least one row.")
  }
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
        stop("Unknown how to handle subject definition from the formula") # nocov
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
  class(ret) <- c("PKNCAconc", class(ret))
  if (missing(exclude)) {
    ret <- setExcludeColumn(ret)
  } else {
    ret <- setExcludeColumn(ret, exclude=exclude)
  }
  if (missing(volume)) {
    ret <- setAttributeColumn(ret, attr_name="volume", default_value=NA_real_)
  } else {
    ret <- setAttributeColumn(ret, attr_name="volume", col_or_value=volume)
    if (!is.numeric(getAttributeColumn(ret, attr_name="volume")[[1]])) {
      stop("Volume must be numeric")
    }
  }
  if (missing(duration)) {
    ret <- setDuration.PKNCAconc(ret)
  } else {
    ret <- setDuration.PKNCAconc(ret, duration=duration)
  }
  if (!missing(time.nominal)) {
    ret <-
      setAttributeColumn(object=ret,
                         attr_name="time.nominal",
                         col_name=time.nominal)
  }
  if (!missing(exclude_half.life)) {
    ret <-
      setAttributeColumn(object=ret,
                         attr_name="exclude_half.life",
                         col_name=exclude_half.life)
  }
  if (!missing(include_half.life)) {
    ret <-
      setAttributeColumn(object=ret,
                         attr_name="include_half.life",
                         col_name=include_half.life)
  }
  ret
}

#' Extract the formula from a PKNCAconc object.
#'
#' @param x The object to extract the formula from.
#' @param \dots Unused
#' @return A formula object
#' @export
#' @importFrom stats formula
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
#' @importFrom stats model.frame
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
#' object).
#' 
#' @param object The object to extract the data from
#' @param form The formula to extract the data from (defaults to the 
#'   formula from \code{object})
#' @param level optional.  If included, this specifies the level(s) of 
#'   the groups to include.  If a numeric scalar, include the first 
#'   \code{level} number of groups.  If a numeric vector, include each 
#'   of the groups specified by the number.  If a character vector, 
#'   include the named group levels.
#' @param data The data to extract the groups from (defaults to the data
#'   from \code{object})
#' @param sep Unused (kept for compatibility with the nlme package)
#' @param ... Arguments passed to other getGroups functions
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
      if (length(level) == 1 &&
          level > 0) {
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
#' @importFrom nlme getData
getData.PKNCAconc <- function(object)
  object$data

#' @rdname getDataName
getDataName.PKNCAconc <- function(object)
  "data"

setDuration.PKNCAconc <- function(object, duration, ...) {
  if (missing(duration)) {
    object <-
      setAttributeColumn(object=object, attr_name="duration", default_value=0,
                         message_if_default="Assuming point rather than interval concentration measurement")
  } else {
    object <-
      setAttributeColumn(object=object, attr_name="duration", col_or_value=duration)
  }
  duration.val <- getAttributeColumn(object=object, attr_name="duration")[[1]]
  if (is.numeric(duration.val) &&
      !any(is.na(duration.val)) &&
      !any(is.infinite(duration.val)) &&
      all(duration.val >= 0)) {
    # It passes the test
  } else {
    stop("duration must be numeric without missing (NA) or infinite values, and all values must be >= 0")
  }
  object
}

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
#' @importFrom stats formula
#' @importFrom utils head
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
  if ("time.nominal" %in% names(x)) {
    cat("Nominal time column is: ", x$time.nominal, "\n", sep="")
  } else {
    cat("Nominal time column is not specified.\n")
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

#' Divide into groups
#' 
#' \code{split.PKNCAconc} divides data into individual groups defined by
#' \code{\link{getGroups.PKNCAconc}}.
#' 
#' @param x the object to split
#' @param f the groups to use for splitting the object
#' @param drop logical indicating if levels that do not occur should be 
#'   dropped.
#' @param ... Ignored.
#' @details If \code{x} is \code{NA} then a list with NA as the only
#'   element and a "groupid" attribute of an empty data.frame is
#'   returned.
#' @return A list of objects with an attribute of groupid consisting of 
#'   a data.frame with columns for each group.
#' @export
split.PKNCAconc <- function(x, f=getGroups(x), drop=TRUE, ...) {
  if (!drop)
    stop("drop must be TRUE")
  if (identical(x, NA)) {
    ret <- list(NA)
    groupid <- data.frame(NA)[,c()]
  } else {
    ## Do the initial separation and extract the groupid information
    f_new <-
      as.character(
        do.call(
          paste,
          append(as.list(f), list(sep="\n"))
        )
      )
    ret <- split(x=x$data, f=f_new, drop=drop, sep="\n")
    groupid <- unique(f)
    ## reorder the output to align with the input grouping order
    ret.idx <-
      factor(
        names(ret),
        levels=do.call(paste, append(as.list(groupid), list(sep="\n"))),
        ordered=TRUE
      )
    ret <- ret[order(ret.idx)]
    ## Reset the data in each split to a "data" element within a list.
    ret <-
      lapply(
        ret,
        function(y, newclass) {
          ret <- list(data=y)
          class(ret) <- newclass
          ret
        },
        newclass=class(x)
      )
    ## Add the other features back into the data
    for (n in setdiff(names(x), "data")) {
      ret <-
        lapply(
          ret,
          function(x, name, value) {
            x[[name]] <- value
            x
          },
          name=n, value=x[[n]]
        )
    }
  }
  attr(ret, "groupid") <- groupid
  ret
}
