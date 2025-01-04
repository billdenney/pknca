#' @importFrom nlme getGroups
#' @export
nlme::getGroups

#' Get the dependent variable (left hand side of the formula) from a
#' PKNCA object.
#'
#' @param x The object to extract the formula from
#' @param \dots Unused
#' @returns The vector of the dependent variable from the object.
#' @family PKNCA object extractors
#' @export
getDepVar <- function(x, ...) {
  UseMethod("getDepVar", x)
}

#' Get the independent variable (right hand side of the formula) from
#' a PKNCA object.
#'
#' @param x The object to extract the formula from
#' @param \dots Unused
#' @returns The vector of the independent variable from the object.
#' @family PKNCA object extractors
#' @export
getIndepVar <- function(x, ...) {
  UseMethod("getIndepVar", x)
}

#' Get the value from a column in a data frame if the value is a column
#' there, otherwise, the value should be a scalar or the length of the
#' data.
#'
#' @param data A data.frame or similar object
#' @param value A character string giving the name of a column in the `data`, a
#'   scalar, or a vector the same length as the `data`
#' @param prefix The prefix to use if a column must be added (it will be used as
#'   the full column name if it is not already in the dataset or it will be
#'   prepended to the maximum column name if not.)
#' @returns A list with elements named "data", "name" giving the `data` with a
#'   column named "name" with the value in that column.
getColumnValueOrNot <- function(data, value, prefix="X") {
  col.name <- setdiff(c(prefix, paste(prefix, max(names(data)), sep=".")), names(data))[1]
  if (is.character(value) && length(value) == 1 && (value %in% names(data))) {
    # It was a column from the data.frame
    ret <- list(data=data, name=value)
  } else if (length(value) %in% c(1, nrow(data))) {
    data[[col.name]] <- value
    ret <- list(data=data, name=col.name)
  } else {
    stop("value was not a column name nor was it a scalar or a vector matching the length of the data.")
  }
  ret
}

#' Get the name of the element containing the data for the current
#' object.
#'
#' @param object The object to get the data name from.
#' @family PKNCA object extractors
#' @return A character scalar with the name of the data object (or `NULL` if the
#'   method does not apply).
#' @keywords Internal
getDataName <- function(object) {
  UseMethod("getDataName")
}

#' @describeIn getDataName If no data name exists, returns NULL.
getDataName.default <- function(object) {
  NULL
}

#' Add an attribute to an object where the attribute is added as a name
#' to the names of the object.
#'
#' @param object The object to set the attribute column on.
#' @param attr_name The attribute name to set
#' @param col_or_value If this exists as a column in the data, it is used as the
#'   `col_name`.  If not, this becomes the `default_value`.
#' @param col_name The name of the column within the dataset to use (if missing,
#'   uses `attr_name`)
#' @param default_value The value to fill in the column if the column does not
#'   exist (the column is filled with `NA` if it does not exist and no value is
#'   provided).
#' @param stop_if_default,warn_if_default,message_if_default A character string
#'   to provide as an error, a warning, or a message to the user if the
#'   `default_value` is used.  They are tested in order (if stop, the code
#'   stops; if warning, the message is ignored; and message last).
#' @returns The object with the attribute column added to the data.
#' @seealso [getAttributeColumn()]
setAttributeColumn <- function(object, attr_name, col_or_value, col_name, default_value,
                               stop_if_default, warn_if_default, message_if_default) {
  dataname <- getDataName(object)
  # Check inputs
  if (!is.character(attr_name) | (length(attr_name) != 1)) {
    stop("attr_name must be a character scalar.")
  }
  if (!missing(col_or_value) &
      any(!c(missing(col_name), missing(default_value)))) {
    stop("Cannot provide col_or_value and col_name or default_value")
  }
  # Apply col_or_value to col_name or to default_value
  if (!missing(col_or_value)) {
    if (all(col_or_value %in% names(object[[dataname]]))) {
      col_name <- col_or_value
    } else {
      default_value <- col_or_value
    }
  }
  # Set the column name
  if (missing(col_name)) {
    col_name <- attr_name
    if (attr_name %in% names(object[[dataname]])) {
      rlang::inform(
        message = paste0("Found column named ", attr_name, ", using it for the attribute of the same name."),
        class = paste0("pknca_foundcolumn_", attr_name)
      )
    }
  } else if (!is.character(col_name) | (length(col_name) != 1)) {
    stop("col_name must be a character scalar.")
  }
  # Set the default value
  if (missing(default_value)) {
    if (col_name %in% names(object[[dataname]])) {
      default_value <- object[[dataname]][[col_name]]
    } else {
      default_value <- NA
      # React to using the default value, if requested
      if (!missing(stop_if_default)) {
        stop(stop_if_default)
      } else if (!missing(warn_if_default)) {
        warning(warn_if_default)
      } else if (!missing(message_if_default)) {
        message(message_if_default)
      }
    }
  }
  # Check that the default_value can work
  if (!(length(default_value) %in% c(1, nrow(object[[dataname]])))) {
    stop("default_value must be a scalar or the same length as the rows in the data.")
  }
  object[[dataname]][[col_name]] <- default_value
  # Inform the object that the column exists
  if (!("columns" %in% names(object))) {
    object$columns <- list()
  }
  object$columns[[attr_name]] <- col_name
  object
}

#' Retrieve the value of an attribute column.
#'
#' @param object The object to extract the attribute value from.
#' @param attr_name The name of the attribute to extract
#' @param warn_missing Give a warning if the "attr"ibute or "column" is missing.
#'   Character vector with zero, one, or both of "attr" and "column".
#' @returns The value of the attribute (or `NULL` if the attribute is not set or
#'   the column does not exist)
getAttributeColumn <- function(object, attr_name, warn_missing=c("attr", "column")) {
  if (length(setdiff(warn_missing, c("attr", "column")))) {
    stop("warn_missing must have a valid value or be empty")
  }
  warn_missing <- warn_missing[warn_missing %in% c("attr", "column")]
  columns <- object$columns[[attr_name]]
  dataname <- getDataName(object)
  if (is.null(columns)) {
    if ("attr" %in% warn_missing)
      warning(attr_name, " is not set.")
    NULL
  } else if (length(missing_cols <- setdiff(columns, names(object[[dataname]])))) {
    if ("column" %in% warn_missing)
      warning("Columns ", paste(missing_cols, collapse=", "), " are not present.")
    NULL
  } else {
    object[[dataname]][, columns, drop=FALSE]
  }
}

#' Check for duplicate values in a dataset
#'
#' @param object A PKNCAconc or PKNCAdose object to check for duplicates
#' @param data_type The name of the type of data for error reporting
#' @returns `object` unmodified, or an error
#'
#' @keywords Internal
#' @noRd
duplicate_check <- function(object, data_type) {
  mask_excluded <- !is.na(object$data[[object$columns$exclude]])
  mask_dup <- rep(FALSE, nrow(object$data))
  key_cols <- unique(c(object$columns$time, unlist(object$columns$groups)))
  if (length(key_cols) == 0) {
    # If there are no key columns, then there can only be one data row that is
    # not excluded.
    mask_dup[!mask_excluded] <- duplicated(mask_dup[!mask_excluded])
  } else {
    # In case an excluded row is the first row of the duplicated set, do not
    # report duplication.
    mask_dup[!mask_excluded] <- duplicated(object$data[!mask_excluded, key_cols])
  }
  if (any(mask_dup)) {
    stop(
      "Rows that are not unique per group and time (column names: ",
      paste(key_cols, collapse=", "),
      ") found within ", data_type, " data.  Row numbers: ",
      paste(which(mask_dup), collapse=", ")
    )
  }
  object
}

#' Set units for a PKNCAconc or PKNCAdose object
#'
#' @param object a PKNCAconc or PKNCAdose object
#' @param units_orig unit specification that may be columns or values
#' @param units_pref unit specification that must be values
#' @returns The object with the units columns definitions optionally added an a
#'   "units" value list added.
#' @noRd
pknca_set_units <- function(object, units_orig = list(), units_pref = list()) {
  all_units <-
    list(
      orig = lapply(X = units_orig, FUN = assert_unit, data = object$data),
      pref = lapply(X = units_pref, FUN = assert_unit_value)
    )

  object$units <- list()
  for (col_units in names(all_units$orig)) {
    current_unit_type <- attr(all_units$orig[[col_units]], "unit_type")
    if (is.null(current_unit_type)) {
      # do nothing
    } else if (current_unit_type %in% "column") {
      object <-
        setAttributeColumn(
          object = object,
          attr_name = col_units,
          col_name = all_units$orig[[col_units]]
        )
    } else if (current_unit_type %in% "value") {
      object$units[[col_units]] <- all_units$orig[[col_units]]
    } else {
      stop(paste("Please report a bug. Unit setting for", col_units)) # nocov
    }
  }
  for (pref_units in names(all_units$pref)) {
    current_unit_type <- attr(all_units$pref[[pref_units]], "unit_type")
    if (is.null(current_unit_type)) {
      # do nothing
    } else if (current_unit_type %in% "value") {
      # you can only set preferred units if you set original units
      original_unit_col <- gsub(x = pref_units, pattern = "_pref", replacement = "")
      if (!(original_unit_col %in% c(names(object$columns), names(object$units)))) {
        stop("Preferred units may not be set unless original units are set: ", pref_units)
      }
      object$units[[pref_units]] <- all_units$pref[[pref_units]]
    } else {
      stop(paste("Please report a bug. Preferred unit setting for", pref_units)) # nocov
    }
  }

  object
}
