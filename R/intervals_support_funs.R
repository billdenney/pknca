#' Add specified imputation methods to the intervals in a PKNCAdata or data.frame object.
#'
#' @param data A PKNCAdata object containing the intervals and data components, or a data frame of intervals.
#' @param target_impute A character string specifying the imputation method to be added.
#' @param after Numeric value specifying the position after which the imputation method should be added (optional). 
#' First is 0, last Inf. If missing, the imputation method is added at the end (Inf).
#' @param target_params A character vector specifying the parameters to be targeted (optional). 
#' If missing, all TRUE in the intervals are taken.
#' @param target_groups A data frame specifying the intervals to be targeted (optional). 
#' If missing, all relevant groups are considered.
#' @param ... arguments passed to `interval_add_impute`.
#' @details. If already present the target_impute method will be added substituting the existing one. All new intervals 
#' created will be added right after their original ones.
#' @return A modified PKNCAdata object with the specified imputation methods added to the targeted intervals.
#' @examples
#' d_conc <- data.frame(
#'   conc = c(1, 0.6, 0.2, 0.1, 0.9, 0.4, 1.2, 0.8, 0.3, 0.2, 1.1, 0.5),
#'   time = rep(0:5, 2),
#'   ID = rep(1:2, each = 6),
#'   analyte = rep(c("Analyte1", "Analyte2"), each = 6)
#' )

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
interval_add_impute <- function(data, ...) {
  UseMethod("interval_add_impute", data)
}

#' @export
interval_add_impute.PKNCAdata <- function(data, target_impute, after = Inf, target_params = NULL, target_groups = NULL, ...) {
  # If the impute column is not present, add it to the intervals
  if (!"impute" %in% names(data$intervals) && !is.null(data$impute)) {
    data$intervals$impute <- data$impute
    data$impute <- NA_character_
  }
  data$intervals <- interval_add_impute.data.frame(data$intervals, target_impute, after, target_params, target_groups)
  data
}

#' Add impute method to the impute column
#'
#' This is an internal helper function used to add an impute method to the impute column.
#'
#' @param impute_vals Character vector of impute methods.
#' @param target_impute The imputation method to be added.
#' @param after Numeric value specifying the position after which the imputation method should be added.
#' @return A character string or vector with the added impute method.
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
#' @export
interval_add_impute.data.frame <- function(intervals, target_impute, after = Inf, target_params = NULL, target_groups = NULL, ...) {
  # Validate inputs
  if (missing(intervals) || missing(target_impute)) {
    stop("Both 'data' and 'target_impute' must be provided.")
  }
  if (!is.character(target_impute)) {
    stop("'target_impute' must be a character string.")
  }

  # Ensure the impute column exists and is a character column
  if (!"impute" %in% colnames(intervals)) {
    intervals$impute <- NA_character_
  } else if (!is.character(intervals$impute)) {
    stop("The 'impute' column in the intervals data.frame must be a character column.")
  }

  # Add an index column to preserve the original order
  index_colname <- make.unique(c("index", names(intervals)))[1]
  intervals[[index_colname]] <- seq_len(nrow(intervals))

  # Get all parameter column names in the data frame
  all_param_options <- setdiff(names(get.interval.cols()), c("start", "end"))
  param_cols <- intersect(names(intervals), all_param_options)

  # If missing, define the target parameters as all parameter columns. Filter based on at least one TRUE value.
  if (is.null(target_params)) {
    target_params <- param_cols
  } else {
    checkmate::assert_subset(target_params, choices = all_param_options, empty.ok = TRUE)
  }

  # Ifentify the target interval rows based on:
  ## 1. The target groups (perfect match)
  target_rows <- rep(TRUE, nrow(intervals))
  if (!is.null(target_groups)) {
    target_groups_data <- intervals[, names(target_groups), drop = FALSE]
    target_rows <- target_rows & (do.call(paste0, target_groups_data) %in% do.call(paste0, target_groups))
  }
  ## 2. The target parameters (at least one calculated: not-FALSE/not-NA)
  target_params_data <- intervals[, target_params, drop = FALSE]
  target_rows <- target_rows & (rowSums(!is.na(replace(target_params_data, target_params_data == FALSE, NA))) > 0)
  ## 3. The target impute method is not already present and correctly positioned
  after_vals <- lapply(strsplit(intervals$impute, "[ ,]+"), function(x) {
    after.x <- which(x == target_impute)
    if (length(after.x) == 0) return(NA)
    if (after.x == length(x)) Inf else after.x
  }) |> unlist()
  target_rows <- target_rows & (after_vals != after | is.na(after_vals))

  new_intervals <- intervals[target_rows, ]

  # If no target intervals are found, nothing to change
  if (nrow(new_intervals) == 0) {
    warning("No intervals found with the specified target parameters, groups, and/or after-change needed. No changes made.")
    return(intervals[, !names(intervals) %in% index_colname])
    
    # If target intervals are found...
  } else {
    # The new imputation should not be used non-target parameters
    new_intervals[, setdiff(param_cols, target_params)] <- NA
    
    # Index the new intervals to be after the original ones
    new_intervals[["impute"]] <- add_impute_method(new_intervals[["impute"]], target_impute, after)
    new_intervals[[index_colname]] <- new_intervals[[index_colname]] + 0.5
  }

  # Remove the target parameters calculation from the original target intervals
  intervals[target_rows, target_params] <- NA

  # Combine the new and original intervals
  intervals <- rbind(intervals, new_intervals)

  # Filter rows where all row values for param_cols are NA or FALSE
  param_data <- intervals[, param_cols, drop = FALSE]
  rows_no_params <- rowSums(!is.na(replace(param_data, param_data == FALSE, NA))) == 0
  intervals <- intervals[!rows_no_params, , drop = FALSE]

  # Order the intervals by the index column and then remove it
  intervals <- intervals[order(intervals[, index_colname]), ]
  rownames(intervals) <- seq_len(nrow(intervals))
  intervals[, !names(intervals) %in% index_colname]
}

#' Remove specified imputation methods from the intervals in a PKNCAdata or data.frame (intervals) object.
#'
#' @inheritParams interval_add_impute
#' @param target_impute A character string specifying the imputation method to remove.
#' @param ... arguments passed to `interval_remove_impute`.
#' @return A modified object with the specified imputation methods removed from the targeted intervals.
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
interval_remove_impute <- function(data, ...) {
  UseMethod("interval_remove_impute", data)
}

#' @export
interval_remove_impute.PKNCAdata <- function(data, target_impute, target_params = NULL, target_groups = NULL) {
  # If the impute column is not present in the intervals...
  ## a. If neither is in the global impute method, return the data as it is
  if (!"impute" %in% colnames(data$intervals) && (is.null(data$impute) || is.na(data$impute))) {
    warning("No default impute column or global method identified. No impute methods to remove")
    return(data)
  } else if (!"impute" %in% names(data$intervals) && (!is.null(data$impute) || !is.na(data$impute))) {
    if (is.null(target_params) && is.null(target_groups)) {
      ## b. If it is in the global impute and no target parameters or groups are specified, remove the global impute method
      data$impute <- remove_impute_method(data$impute, target_impute)
      return(data)
    } else {
      ## c. If it is in the global impute but target parameters or groups are specified, create a impute column in the intervals
      data$intervals$impute <- data$impute
      data$impute <- NA_character_
    }
  }
  data$intervals <- interval_remove_impute.data.frame(data$intervals, target_impute, target_params, target_groups)
  data
}

#' Remove impute method from the impute column
#'
#' This is an internal helper function used to remove an impute method from the impute column.
#'
#' @param impute_vals Character vector of impute methods.
#' @param target_impute The imputation method to be removed.
#' @return A character string or vector without the specified impute method.
#' @details Resulting empty string values are replaced with NA_character_.
#' @keywords internal
remove_impute_method <- function(impute_vals, target_impute){
  # Make sure the character vector has length
  if (length(impute_vals) == 0) return(impute_vals)
  
  # Remove the impute from the other methods in each value
  impute_vals <- ifelse(is.na(impute_vals), "", impute_vals)
  strsplit(impute_vals, split = "[ ,]+") |>
    lapply(FUN = setdiff, target_impute) |>
    vapply(FUN = paste, collapse = ",", FUN.VALUE = "")
}

#' @export
interval_remove_impute.data.frame <- function(intervals, target_impute, target_params = NULL, target_groups = NULL) {
  # Validate inputs
  if (missing(intervals) || missing(target_impute)) {
    stop("Both 'data' and 'target_impute' must be provided.")
  }
  if (!is.character(target_impute)) {
    stop("'target_impute' must be a character string.")
  }

  # Ensure the impute column exists and is a character column
  if (!"impute" %in% colnames(intervals)) {
    warning("No default impute column identified. No impute methods to remove")
    return(intervals)
  } else if (!is.character(intervals$impute)) {
    stop("The 'impute' column in the intervals data.frame must be a character column.")
  }

  # Add an index column to preserve the original order
  index_colname <- make.unique(c("index", names(intervals)))[1]
  intervals[[index_colname]] <- seq_len(nrow(intervals))

  # Get all parameter column names in the data frame
  all_param_options <- setdiff(names(get.interval.cols()), c("start", "end"))
  param_cols <- intersect(names(intervals), all_param_options)

  # Handle target_params
  if (is.null(target_params)) {
    target_params <- param_cols
  } else {
    checkmate::assert_subset(target_params, choices = all_param_options)
  }

  # Ifentify the target interval rows based on:
  ## 1. The target groups (perfect match)
  target_rows <- rep(TRUE, nrow(intervals))
  if (!is.null(target_groups)) {
    target_groups_data <- intervals[, names(target_groups), drop = FALSE]
    target_rows <- target_rows & (do.call(paste0, target_groups_data) %in% do.call(paste0, target_groups))
  }
  ## 2. The target parameters (at least one calculated: not-FALSE/not-NA)
  target_params_data <- intervals[, target_params, drop = FALSE]
  target_rows <- target_rows & (rowSums(!is.na(replace(target_params_data, target_params_data == FALSE, NA))) > 0)
  ## 3. The target impute method to be removed (contained in the string)
  target_rows <- target_rows & grepl(target_impute, intervals$impute, fixed = TRUE)
  new_intervals <- intervals[target_rows, ]
  
  # If no target intervals are found, nothing to change
  if (nrow(new_intervals) == 0) {
    warning("No intervals found with the specified target parameters, groups and/or impute method. No changes made.")
    return(intervals[, !names(intervals) %in% index_colname])

  # If target intervals are found...
  } else {
    # The new imputation should not involve non-target parameters
    new_intervals[, setdiff(param_cols, target_params)] <- NA
    
    # Index the new intervals to be after the original ones
    new_intervals[["impute"]] <- remove_impute_method(new_intervals[["impute"]], target_impute)
    new_intervals[[index_colname]] <- new_intervals[[index_colname]] + 0.5
  }

  # Remove the target parameters calculation from the original target intervals
  intervals[target_rows, target_params] <- NA

  # Combine the new and original intervals
  intervals <- rbind(intervals, new_intervals)

  # Filter rows where all row values for param_cols are NA or FALSE
  param_data <- intervals[, param_cols, drop = FALSE]
  rows_no_params <- rowSums(!is.na(replace(param_data, param_data == FALSE, NA))) == 0
  intervals <- intervals[!rows_no_params, , drop = FALSE]

  # Order the intervals by the index column and then remove it
  intervals <- intervals[order(intervals[, index_colname]), ]
  rownames(intervals) <- seq_len(nrow(intervals))
  intervals[, !names(intervals) %in% index_colname]
}
