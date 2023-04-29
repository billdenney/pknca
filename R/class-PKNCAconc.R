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
#' @param sparse Are the concentration-time data sparse PK (commonly used in
#'   small nonclinical species or with terminal or difficult sampling) or dense
#'   PK (commonly used in clinical studies or larger nonclinical species)?
#' @param ... Ignored.
#' @return A PKNCAconc object that can be used for automated NCA.
#' @family PKNCA objects
#' @export
PKNCAconc <- function(data, ...) {
  UseMethod("PKNCAconc")
}

#' @rdname PKNCAconc
#' @export
PKNCAconc.default <- function(data, ...) {
  PKNCAconc.data.frame(as.data.frame(data), ...)
}
#' @rdname PKNCAconc
#' @export
PKNCAconc.tbl_df <- function(data, ...) {
  PKNCAconc.data.frame(as.data.frame(data), ...)
}

#' @rdname PKNCAconc
#' @export
PKNCAconc.data.frame <- function(data, formula, subject,
                                 time.nominal, exclude, duration, volume,
                                 exclude_half.life, include_half.life, sparse=FALSE, ...) {
  # The data must have... data
  if (nrow(data) == 0) {
    stop("data must have at least one row.")
  }
  # Verify that all the variables in the formula are columns in the data.
  missing_vars <- setdiff(all.vars(formula), names(data))
  if (length(missing_vars) > 0) {
    stop("All of the variables in the formula must be in the data.  Missing: ", paste(missing_vars))
  }
  parsed_form_raw <- parse_formula_to_cols(form = formula)
  parsed_form_groups <-
    if (length(parsed_form_raw$groups) > 0) {
      list(
        group_vars=parsed_form_raw$groups,
        group_analyte=character()
      )
    } else {
      list(
        group_vars=parsed_form_raw$groups_left_of_slash,
        group_analyte=parsed_form_raw$groups_right_of_slash
      )
    }
  parsed_form <-
    list(
      concentration = parsed_form_raw$lhs,
      time = parsed_form_raw$rhs,
      groups = parsed_form_groups
    )
  if (length(parsed_form$concentration) != 1) {
    stop("The left hand side of the formula must have exactly one variable")
  }
  if (length(parsed_form$time) != 1) {
    stop("The right hand side of the formula (excluding groups) must have exactly one variable")
  }
  # Do some general checking of the concentration and time data to give an early
  # error if the data are not correct.  Do not check monotonic.time because the
  # data may contain information for more than one subject.
  check.conc.time(
    conc=data[[parsed_form$concentration]],
    time=data[[parsed_form$time]],
    monotonic.time=FALSE
  )
  # Values must be unique (one value per measurement)
  key_cols <- c(parsed_form$time, unlist(parsed_form$groups))
  mask_dup <- duplicated(data[,key_cols])
  if (any(mask_dup)) {
    stop("Rows that are not unique per group and time (column names: ",
         paste(key_cols, collapse=", "),
         ") found within concentration data.  Row numbers: ",
         paste(seq_along(mask_dup)[mask_dup], collapse=", "))
  }
  # Assign the subject
  if (missing(subject)) {
    subject <- parsed_form$groups$group_vars[length(parsed_form$groups$group_vars)]
  } else {
    # Ensure that the subject is part of the data definition and a scalar
    # character string.
    if (!is.character(subject))
      stop("subject must be a character string")
    if (!(length(subject) == 1))
      stop("subject must be a scalar")
    if (!(subject %in% names(data)))
      stop("The subject parameter must map to a name in the data")
  }
  parsed_form$subject <- subject
  if (sparse) {
    ret <-
      list(
        data_sparse = data,
        formula = formula,
        columns = parsed_form
      )
  } else {
    ret <-
      list(
        data = data,
        formula = formula,
        columns = parsed_form
      )
  }
  class(ret) <- c("PKNCAconc", class(ret))
  if (missing(exclude)) {
    ret <- setExcludeColumn(ret, dataname=getDataName.PKNCAconc(ret))
  } else {
    ret <- setExcludeColumn(ret, exclude=exclude, dataname=getDataName.PKNCAconc(ret))
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
formula.PKNCAconc <- function(x, ...) {
  x$formula
}

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
#' @method model.frame PKNCAconc
model.frame.PKNCAconc <- function(formula, ...) {
  formula$data[, all.vars(formula$formula), drop=FALSE]
}

#' @export
getDepVar.PKNCAconc <- function(x, ...) {
  x$data[, x$columns$concentration]
}

#' @export
getIndepVar.PKNCAconc <- function(x, ...) {
  x$data[, x$columns$time]
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
getGroups.PKNCAconc <- function(object, form=stats::formula(object), level,
                                data=as.data.frame(object), sep) {
  grpnames <- unlist(object$columns$groups)
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

#' Get grouping variables for a PKNCA object
#'
#' @param x The PKNCA object
#' @return A character vector (possibly empty) of the grouping variables
#' @exportS3Method dplyr::group_vars
group_vars.PKNCAconc <- function(x) {
  unname(unlist(x$columns$groups))
}

#' @rdname getDataName
#' @export
getDataName.PKNCAconc <- function(object) {
  if (is_sparse_pk(object)) {
    "data_sparse"
  } else {
    "data"
  }
}

setDuration.PKNCAconc <- function(object, duration, ...) {
  if (missing(duration)) {
    object <-
      setAttributeColumn(
        object=object,
        attr_name="duration",
        default_value=0,
        message_if_default="Assuming point rather than interval concentration measurement"
      )
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
print.PKNCAconc <- function(x, n=6, summarize=FALSE, ...) {
  cat(sprintf("Formula for concentration:\n "))
  print(stats::formula(x), ...)
  if (is_sparse_pk(x)) {
    data_current <- x$data_sparse
    is_sparse <- TRUE
    cat("Data are sparse PK.\n")
  } else {
    data_current <- x$data
    is_sparse <- FALSE
    cat("Data are dense PK.\n")
  }
  single_subject <- is.na(x$columns$subject) || (length(x$columns$subject) == 0)
  if (single_subject) {
    cat("As a single-subject dataset.\n")
  } else {
    cat(sprintf("With %d subjects defined in the '%s' column.\n",
                length(unique(data_current[,x$columns$subject])),
                x$columns$subject))
  }
  if ("time.nominal" %in% names(x$columns)) {
    cat("Nominal time column is: ", x$columns$time.nominal, "\n", sep="")
  } else {
    cat("Nominal time column is not specified.\n")
  }
  if (summarize) {
    cat("\n")
    grp <- getGroups(x)
    if (ncol(grp) > 0) {
      tmp_summary <- data.frame(Group.Name=names(grp),
                                Count=0)
      for (i in seq_len(ncol(grp))) {
        tmp_summary$Count[i] <- nrow(unique(grp[,1:i,drop=FALSE]))
      }
      cat("Group summary:\n")
      names(tmp_summary) <- gsub("\\.", " ", names(tmp_summary))
      print.data.frame(tmp_summary, row.names=FALSE)
    } else {
      cat("No groups.\n")
    }
  }
  if (n != 0) {
    if (n >= nrow(data_current)) {
      cat("\nData for concentration:\n")
    } else if (n < 0) {
      cat(sprintf("\nFirst %d rows of concentration data:\n",
                  nrow(data_current)+n))
    } else {
      cat(sprintf("\nFirst %d rows of concentration data:\n",
                  n))
    }
    print.data.frame(utils::head(data_current, n=n), ..., row.names=FALSE)
  }
}

#' @rdname is_sparse_pk
#' @export
is_sparse_pk.PKNCAconc <- function(object) {
  "data_sparse" %in% names(object)
}

#' @rdname print.PKNCAconc
#' @export
summary.PKNCAconc <- function(object, n=0, summarize=TRUE, ...) {
  print.PKNCAconc(object, n=n, summarize=summarize)
}

#' @export
as.data.frame.PKNCAconc <- function(x, ...) {
  if (is_sparse_pk(x)) {
    x$data_sparse
  } else {
    x$data
  }
}
