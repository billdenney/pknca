#' Add specified imputation methods to the intervals in a PKNCAdata or data.frame object.
#'
#' @param data A PKNCAdata object containing the intervals data frame, or a data frame of intervals.
#' @param target_impute A character string specifying the imputation method to be added.
#' @param after Numeric value specifying the index in which the imputation will be added (optional).
#' First is 0, last Inf. If missing, the imputation method is added at the end (Inf).
#' @param target_params A character vector specifying the parameters to be targeted (optional).
#' If missing, all TRUE in the intervals are taken.
#' @param target_groups A data frame specifying the intervals to be targeted (optional).
#' If missing, all relevant groups are considered.
#' @param ... arguments passed to `interval_add_impute`.
#' @details If already present the target_impute method will be added substituting the existing one.
#' All new intervals created will be added right after their original ones.
#' @returns A modified PKNCAdata object with specified imputation methods on the target intervals.
#' @examples
#' d_conc <- data.frame(
#'   conc = c(1, 0.6, 0.2, 0.1, 0.9, 0.4, 1.2, 0.8, 0.3, 0.2, 1.1, 0.5),
#'   time = rep(0:5, 2),
#'   ID = rep(1:2, each = 6),
#'   analyte = rep(c("Analyte1", "Analyte2"), each = 6)
#' )
#'
#' d_dose <- data.frame(
#'   dose = c(100, 200),
#'   time = c(0, 0),
#'   ID = c(1, 2)
#' )
#'
#' o_conc <- PKNCAconc(d_conc, conc ~ time | ID / analyte)
#' o_dose <- PKNCAdose(d_dose, dose ~ time | ID)
#'
#' intervals <- data.frame(
#'   start = c(0, 0, 0),
#'   end = c(3, 5, Inf),
#'   half.life = c(TRUE, TRUE, TRUE),
#'   cmax = c(TRUE, TRUE, TRUE),
#'   impute = c("start_conc0,start_predose", "start_predose", "start_conc0"),
#'   analyte = c("Analyte1", "Analyte2", "Analyte1")
#' )
#'
#' o_data <- PKNCAdata(o_conc, o_dose, intervals = intervals)
#'
#' # Apply interval_add_impute function
#' o_data <- interval_add_impute(o_data,
#'                               target_impute = "start_conc0",
#'                               target_params = "half.life",
#'                               target_groups = data.frame(analyte = "Analyte1"))
#' @export
interval_add_impute <- function(data, target_impute, after, target_params, target_groups, ...) {
  UseMethod("interval_add_impute", data)
}

#' Remove specified imputation from the intervals in a PKNCAdata or data.frame (intervals) object.
#'
#' @inheritParams interval_add_impute
#' @param target_impute A character string specifying the imputation method to remove.
#' @param ... arguments passed to `interval_remove_impute`.
#' @returns A modified object with the specified imputations removed from the targeted intervals.
#' @examples
#' d_conc <- data.frame(
#'   conc = c(1, 0.6, 0.2, 0.1, 0.9, 0.4, 1.2, 0.8, 0.3, 0.2, 1.1, 0.5),
#'   time = rep(0:5, 2),
#'   ID = rep(1:2, each = 6),
#'   analyte = rep(c("Analyte1", "Analyte2"), each = 6)
#' )
#'
#' d_dose <- data.frame(
#'   dose = c(100, 200),
#'   time = c(0, 0),
#'   ID = c(1, 2)
#' )
#'
#' o_conc <- PKNCAconc(d_conc, conc ~ time | ID / analyte)
#' o_dose <- PKNCAdose(d_dose, dose ~ time | ID)
#'
#' intervals <- data.frame(
#'   start = c(0, 0, 0),
#'   end = c(3, 5, Inf),
#'   half.life = c(TRUE, FALSE, TRUE),
#'   cmax = c(TRUE, TRUE, TRUE),
#'   impute = c("start_conc0,start_predose", "start_predose", "start_conc0"),
#'   analyte = c("Analyte1", "Analyte2", "Analyte1")
#' )
#'
#' o_data <- PKNCAdata(o_conc, o_dose, intervals = intervals)
#'
#' # Apply interval_remove_impute function
#' o_data <- interval_remove_impute(data = o_data,
#'                                  target_impute = "start_conc0",
#'                                  target_params = "half.life",
#'                                  target_groups = data.frame(analyte = "Analyte1"))
#' @export
interval_remove_impute <- function(data, target_impute, ...) {
  UseMethod("interval_remove_impute", data)
}

#' @export
interval_add_impute.PKNCAdata <- function(data, target_impute, after = Inf,
                                          target_params = NULL, target_groups = NULL, ...) {
  # If the impute column is not present, add it to the intervals
  if (!"impute" %in% names(data$intervals) && !is.null(data$impute)) {
    data$intervals$impute <- data$impute
    data$impute <- NA_character_
  }
  data$intervals <- interval_add_impute.data.frame(data$intervals, target_impute,
                                                   after, target_params, target_groups)
  data
}

#' @export
interval_remove_impute.PKNCAdata <- function(data, target_impute, target_params = NULL,
                                             target_groups = NULL, ...) {
  # If the impute column is not present in the intervals...
  if (!"impute" %in% colnames(data$intervals)) {
    if (is.null(data$impute) || is.na(data$impute)) {
      # a. If there is neither a global impute, return the input as it is (nothing to remove)
      warning("No default impute column or global method identified. No impute methods to remove")
      return(data)
    }

    # b & c. If there is a global impute..
    if (is.null(target_params) && is.null(target_groups)) {
      # b. and user changes apply to all intervals, just remove global impute
      data$impute <- remove_impute_method(data$impute, target_impute)
      return(data)
    }

    # c. but user changes are specific (target_params or target_groups), creates an impute column
    data$intervals$impute <- data$impute
    data$impute <- NA_character_
  }

  data$intervals <- interval_remove_impute.data.frame(data$intervals, target_impute,
                                                      target_params, target_groups)
  data
}

#' @export
interval_add_impute.data.frame <- function(data, target_impute, after = Inf,
                                           target_params = NULL, target_groups = NULL, ...) {
  # Validate inputs
  if (missing(data) || missing(target_impute)) {
    stop("Both 'data' and 'target_impute' must be provided.")
  }
  if (!is.character(target_impute)) {
    stop("'target_impute' must be a character string.")
  }
  if (is.na(target_impute) || target_impute == "") {
    warning("No impute method specified. No changes made.")
    return(data)
  }

  # Ensure the impute column exists and is a character column
  if (!"impute" %in% colnames(data)) {
    data$impute <- NA_character_
  } else if (!is.character(data$impute)) {
    stop("The 'impute' column in the intervals data.frame must be a character column.")
  }

  # Add an index column to preserve the original order
  index_colname <- make.unique(c("index", names(data)))[1]
  data[[index_colname]] <- seq_len(nrow(data))

  # Get all parameter column names in the data frame
  all_param_options <- setdiff(names(get.interval.cols()), c("start", "end"))
  param_cols <- intersect(names(data), all_param_options)

  # If missing, define target parameters as all parameter columns with at least one TRUE.
  if (is.null(target_params)) {
    target_params <- param_cols
  }

  assert_subset(target_params, all_param_options)

  # Identify the target interval rows based on:
  target_rows <- identify_target_rows(data, target_impute, target_params, target_groups, after)
  new_intervals <- data[target_rows, ]

  # If no target intervals are found, nothing to change
  if (nrow(new_intervals) == 0) {
    warning("No intervals found with the specified target parameters,",
            " groups, and/or after-change needed. No changes made.")
    return(data[, !names(data) %in% index_colname])

    # If target intervals are found...
  } else {
    # The new imputation should not be used non-target parameters
    new_intervals[, setdiff(param_cols, target_params)] <- NA

    # Index the new intervals to be after the original ones
    new_intervals[["impute"]] <- add_impute_method(new_intervals[["impute"]], target_impute, after)
    new_intervals[[index_colname]] <- new_intervals[[index_colname]] + 0.5
  }

  # Remove the target parameters calculation from the original target intervals
  data[target_rows, target_params] <- NA

  # Combine the new and original intervals
  data <- rbind(data, new_intervals)

  # Filter rows where all row values for param_cols are NA or FALSE
  param_data <- data[, param_cols, drop = FALSE]
  rows_no_params <- rowSums(replace(param_data, is.na(param_data), FALSE)) == 0
  data <- data[!rows_no_params, , drop = FALSE]

  # Order the intervals by the index column and then remove it
  data <- data[order(data[[index_colname]]), ]
  rownames(data) <- seq_len(nrow(data))
  data[, !names(data) %in% index_colname]
}

#' @export
interval_remove_impute.data.frame <- function(data,
                                              target_impute,
                                              target_params = NULL,
                                              target_groups = NULL,
                                              ...) {
  # Validate inputs
  if (missing(data) || missing(target_impute)) {
    stop("Both 'data' and 'target_impute' must be provided.")
  }
  if (!is.character(target_impute)) {
    stop("'target_impute' must be a character string.")
  }
  if (is.na(target_impute) || target_impute == "") {
    warning("No impute method specified. No changes made.")
    return(data)
  }

  # Ensure the impute column exists and is a character column
  if (!"impute" %in% colnames(data)) {
    warning("No default impute column identified. No impute methods to remove")
    return(data)
  } else if (!is.character(data$impute)) {
    stop("The 'impute' column in the intervals data.frame must be a character column.")
  }

  # Add an index column to preserve the original order
  index_colname <- make.unique(c("index", names(data)))[1]
  data[[index_colname]] <- seq_len(nrow(data))

  # Get all parameter column names in the data frame
  all_param_options <- setdiff(names(get.interval.cols()), c("start", "end"))
  param_cols <- intersect(names(data), all_param_options)

  # Handle target_params
  if (is.null(target_params)) {
    target_params <- param_cols
  }

  assert_subset(target_params, all_param_options)

  # Identify the interval rows that need to be changed
  target_rows <- identify_target_rows(data, target_impute, target_params, target_groups)
  new_intervals <- data[target_rows, ]

  # If no target intervals are found, nothing to change
  if (nrow(new_intervals) == 0) {
    warning(paste0("No intervals found with the specified target parameters,",
                   " groups and/or impute method. No changes made."))
    return(data[, !names(data) %in% index_colname])

    # If target intervals are found...
  } else {
    # The new imputation should not involve non-target parameters
    new_intervals[, setdiff(param_cols, target_params)] <- NA

    # Index the new intervals to be after the original ones
    new_intervals[["impute"]] <- remove_impute_method(new_intervals[["impute"]], target_impute)
    new_intervals[[index_colname]] <- new_intervals[[index_colname]] + 0.5
  }

  # Remove the target parameters calculation from the original target intervals
  data[target_rows, target_params] <- NA

  # Combine the new and original intervals
  data <- rbind(data, new_intervals)

  # Filter rows where all row values for param_cols are NA/FALSE
  param_data <- data[, param_cols, drop = FALSE]
  rows_no_params <- rowSums(replace(param_data, is.na(param_data), FALSE)) == 0
  data <- data[!rows_no_params, , drop = FALSE]

  # Order the intervals by the index column and then remove it
  data <- data[order(data[[index_colname]]), ]
  rownames(data) <- seq_len(nrow(data))
  data[, !names(data) %in% index_colname]
}

#' Add impute method to the impute column
#'
#' This is an internal helper function used to add an impute method to the impute column.
#'
#' @param impute_vals Character vector of impute methods.
#' @param target_impute The imputation method to be added.
#' @param after Numeric value specifying the index position in which to add the impute.
#' @returns A character string or vector with the added impute method.
#' @keywords internal
add_impute_method <- function(impute_vals, target_impute, after) {
  # Make sure the character vector has length
  if (length(impute_vals) == 0) return(impute_vals)

  # Remove the impute from the other methods in each value
  impute_vals <- ifelse(is.na(impute_vals), "", impute_vals)
  strsplit(impute_vals, split = "[ ,]+") |>
    lapply(FUN = setdiff, target_impute) |>
    lapply(FUN = append, values = target_impute, after = after) |>
    vapply(FUN = paste, collapse = ",", FUN.VALUE = "")
}

#' Remove impute method from the impute column
#'
#' This is an internal helper function used to remove an impute method from the impute column.
#'
#' @param impute_vals Character vector of impute methods.
#' @param target_impute The imputation method to be removed.
#' @returns A character string or vector without the specified impute method.
#' @details Resulting empty string values are replaced with NA_character_.
#' @keywords internal
remove_impute_method <- function(impute_vals, target_impute) {
  # Make sure the character vector has length
  if (length(impute_vals) == 0) return(impute_vals)

  # Remove the impute from the other methods in each value
  impute_vals <- ifelse(is.na(impute_vals), "", impute_vals)
  impute_vals <- strsplit(impute_vals, split = "[ ,]+") |>
    lapply(FUN = setdiff, target_impute) |>
    vapply(FUN = paste, collapse = ",", FUN.VALUE = "")

  # Replace empty strings with NA_character_
  ifelse(impute_vals == "", NA_character_, impute_vals)
}

#' Identify target rows based on groups, parameters, and impute method
#'
#' This is an internal helper function used to identify the target rows in the data frame
#' based on the specified groups, parameters, and impute method.
#'
#' @param data A data frame containing the intervals.
#' @param target_impute The imputation method to be added or removed.
#' @param target_params A character vector specifying the parameters to be targeted.
#' @param target_groups A data frame specifying the intervals to be targeted.
#' @param after Numeric value specifying the index position in which to add the impute (optional).
#' @returns A logical vector indicating the target rows.
#' @keywords internal
identify_target_rows <- function(data, target_impute, target_params, target_groups, after = NULL) {
  # Identify the target interval rows based on:
  ## 1. The target groups (perfect match)
  is_target_group <- {
    if (!is.null(target_groups)) {
      sapply(data[, names(target_groups), drop = FALSE], paste0) %in% sapply(target_groups, paste0)
    } else {
      rep(TRUE, nrow(data))
    }
  }

  ## 2. The target parameters (at least one calculated: not-FALSE/not-NA)
  target_params_data <- data[, target_params, drop = FALSE]
  is_target_param <- rowSums(replace(target_params_data, is.na(target_params_data), FALSE)) > 0

  ## 3. The target impute method is not present and correctly positioned (if after is provided)
  if (!is.null(after)) {
    after_vals <- sapply(strsplit(data$impute, "[ ,]+"), \(x) {
      after_x <- which(x == target_impute)
      if (length(after_x) == 0) return(TRUE)
      if (after_x == length(x)) Inf else after_x
    })
    is_after <- after_vals != after | is.na(after_vals)
  } else {
    is_after <- grepl(target_impute, data$impute, fixed = TRUE)
  }

  is_target_group & is_target_param & is_after
}

#' Checks if a vector is a subset of another. If there are any values in `a` that are not present
#' in `b`, throws an error.
#' @param a Vector to check.
#' @param b Vector with possible values.
#' @noRd
assert_subset <- function(a, b) {
  if (!all(a %in% b)) {
    stop(
      "The following parameters are invalid interval columns: ",
      paste0(setdiff(a, b), collapse = ", ")
    )
  }
}
