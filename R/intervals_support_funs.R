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
#' @param impute_column A character string specifying the name of the impute column (optional). 
#' If missing, the default name "impute" is used.
#' @details. If already present the target_impute method will be added substituting the existing one. All new intervals 
#' created will be added right after their original ones.
#' @return A modified PKNCAdata object with the specified imputation methods added to the targeted intervals.
#' @examples
#' d_conc <- data.frame(
#'   conc = c(1, 0.6, 0.2, 0.1, 0.9, 0.4, 1.2, 0.8, 0.3, 0.2, 1.1, 0.5),
#'   time = rep(0:5, 2),
#'   ID = rep(1:2, each = 6),
#'   analyte = rep(c("Analyte1", "Analyte2"), each = 6),
#'   include_hl = c(FALSE, NA, TRUE, TRUE, TRUE, TRUE, FALSE, NA, TRUE, TRUE, TRUE, TRUE)
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
#' # Apply interval_add_impute function
#' o_data <- interval_add_impute(o_data, target_impute = "start_conc0", target_params = c("half.life"), target_groups = data.frame(analyte = "Analyte1"))
#' @export
interval_add_impute <- function(data, target_impute, ...) {
  UseMethod("interval_add_impute", data)
}

#' @export
interval_add_impute.PKNCAdata <- function(data, target_impute, after = Inf, target_params = NULL, target_groups = NULL) {
  # If the impute column is not present, add it to the intervals
  if (!"impute" %in% names(data$intervals) && !is.null(data$impute)) {
    data$intervals$impute <- data$impute
    data$impute <- NA_character_
  } else if (class(data$intervals$impute) != "character") {
    stop("The 'impute' column in the intervals must be a character column.")
  }
  data$intervals <- interval_add_impute.data.frame(data$intervals, target_impute, after, target_params, target_groups)
  data
}

# Helper function to process impute methods
add_impute_method <- Vectorize(function(impute_col_value, target_impute, after) {
  impute_methods <- unlist(strsplit(ifelse(is.na(impute_col_value), "", impute_col_value), "[ ,]+")) |>
    setdiff(target_impute) |>
    append(target_impute, after) |>
    paste(collapse = ",")
}, vectorize.args = "impute_col_value", USE.NAMES = FALSE)

#' @export
interval_add_impute.data.frame <- function(intervals, target_impute, after = Inf, target_params = NULL, target_groups = NULL) {
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
    stop("The 'impute' column in the data frame must be a character column.")
  }
  
  # Add an index column to preserve the original order
  index_colname <- make.unique(c("index", names(intervals)))[1]
  intervals[[index_colname]] <- 1:nrow(intervals)
  
  # Get all parameter column names in the data frame
  all_param_options <- setdiff(names(get.interval.cols()), c("start", "end"))
  param_cols <- intersect(names(intervals), all_param_options)
  
  # If missing, define the target parameters as all parameter columns. Filter based on at least one TRUE value.
  if (is.null(target_params)) {
    target_params <- param_cols
  } else {
    checkmate::assert_subset(target_params, choices = all_param_options, empty.ok = TRUE)
  }
  target_param_rows <- rowSums(!is.na(replace(intervals[, target_params, drop = FALSE], FALSE, NA))) > 0
  
  # If missing, define the target groups as all intervals. Filter based on a perfect row match.
  if (!is.null(target_groups)) {
    target_group_rows <- do.call(paste0, intervals[, names(target_groups), drop = FALSE]) %in% do.call(paste0, target_groups)
  } else {
    target_group_rows <- rep(TRUE, nrow(intervals))
  }
  
  # Combine the two conditions to get the final targeted rows (filter for intervals)
  target_rows <- target_group_rows & target_param_rows
  
  # Create new intervals for the target parameters including the impute method and indexed right after the original intervals
  new_intervals <- intervals[target_rows, ]
  new_intervals[, setdiff(param_cols, target_params)] <- NA
  new_intervals[["impute"]] <- add_impute_method(new_intervals[["impute"]], target_impute, after)
  new_intervals[[index_colname]] <- new_intervals[[index_colname]] + 0.5
  
  # Keep the original intervals without the target parameters and without the impute method
  original_intervals <- intervals[target_rows, ]
  original_intervals[, target_params] <- NA
  
  # Combine non-modified intervals, new intervals, and original intervals
  intervals <- rbind(intervals[!target_rows, ], original_intervals, new_intervals)
  
  # Filter rows where all row values for param_cols are NA or FALSE
  param_data <- intervals[, param_cols, drop = FALSE]
  rows_no_params <- rowSums(!is.na(replace(param_data, param_data == FALSE, NA))) == 0
  intervals <- intervals[!rows_no_params, , drop = FALSE]
  
  # Order the data by the index column and then remove it
  intervals <- intervals[order(intervals[, index_colname]), ]
  rownames(intervals) <- 1:nrow(intervals)
  intervals[, !names(intervals) %in% index_colname]
}

#' Remove specified imputation methods from the intervals in a PKNCAdata or data.frame object.
#'
#' @inheritParams interval_add_impute
#' @param target_impute A character string specifying the imputation method to be removed.
#' @return A modified object with the specified imputation methods removed from the targeted intervals.
#' @examples
#' d_conc <- data.frame(
#'   conc = c(1, 0.6, 0.2, 0.1, 0.9, 0.4, 1.2, 0.8, 0.3, 0.2, 1.1, 0.5),
#'   time = rep(0:5, 2),
#'   ID = rep(1:2, each = 6),
#'   analyte = rep(c("Analyte1", "Analyte2"), each = 6),
#'   include_hl = c(FALSE, NA, TRUE, TRUE, TRUE, TRUE, FALSE, NA, TRUE, TRUE, TRUE, TRUE)
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
#' o_data <- interval_remove_impute(data = o_data, target_impute = "start_conc0", target_params = c("half.life"), target_groups = data.frame(analyte = "Analyte1"))
#' @export
interval_remove_impute <- function(data, target_impute, target_params = NULL, target_groups = NULL, impute_column = NULL, new_rows_after_original = TRUE) {
  if (missing(data) || missing(target_impute)) {
    stop("Both 'data' and 'target_impute' must be provided.")
  }
  if (!inherits(data, "PKNCAdata") && !is.data.frame(data)) {
    stop("The 'data' object must be a PKNCAdata object or a data frame.")
  }
  if (!is.character(target_impute)) {
    stop("'target_impute' must be a character string.")
  }
  UseMethod("interval_remove_impute")
}

#' @export
interval_remove_impute.PKNCAdata <- function(data, target_impute, target_params = NULL, target_groups = NULL, impute_column = NULL, new_rows_after_original = TRUE) {
  if (is.null(impute_column) && !is.na(data$impute)) {
    impute_column <- data$impute
  }
  data$intervals <- interval_remove_impute(data$intervals, target_impute, target_params, target_groups, impute_column, new_rows_after_original)
  data
}

#' @export
interval_remove_impute.data.frame <- function(data, target_impute, target_params = NULL, target_groups = NULL, impute_column = NULL, new_rows_after_original = TRUE) {
  
  # Add an index column to preserve the original order
  data <- dplyr::mutate(data, index = dplyr::row_number())
  
  # Get all parameter column names in the data frame
  all_param_options <- setdiff(names(get.interval.cols()), c("start", "end"))
  param_cols <- intersect(names(data), all_param_options)
  
  # Handle target_params
  if (is.null(target_params)) {
    # Take all parameter columns present with at least one TRUE value
    target_params <- param_cols[colSums(data[param_cols]) > 0]
  } else {
    # Check that all target_params are present in the data frame
    missing_params <- setdiff(target_params, param_cols)
    if (length(missing_params) > 0) {
      stop("The following target_params are not interval columns and/or known PKNCA parameters: ", paste(missing_params, collapse = ", "))
      target_params <- intersect(target_params, param_cols)
    }
  }
  
  # Determine the name of the impute column
  impute_col <- if (!is.null(impute_column)) {
    if (!impute_column %in% colnames(data)) {
      stop("The 'data' object does not contain the specified impute column.")
    }
    impute_column
  } else if ("impute" %in% colnames(data)) {
    "impute"
  } else {
    warning("No default impute column identified. No impute methods to remove. If there is an impute column, please specify it in argument 'impute_column'")
    return(data %>% dplyr::select(-index))
  }
  
  # Identify the targeted intervals based on the groups
  if (!is.null(target_groups)) {
    target_intervals <- dplyr::inner_join(data, target_groups, by = names(target_groups))
  } else {
    target_intervals <- data
  }
  
  # Identify the targeted intervals based on the impute method and parameters
  target_intervals <- target_intervals %>%
    dplyr::filter(rowSums(dplyr::across(dplyr::any_of(target_params), ~ . == TRUE)) > 0) %>%
    dplyr::filter(grepl(
      pattern = paste0(".*(", paste0(target_impute, collapse = ")|("), ").*"),
      .data[[impute_col]]
    ))
  
  # Create the new version intervals only for the target parameters
  new_intervals_without_impute <- target_intervals %>%
    dplyr::mutate(dplyr::across(dplyr::any_of(param_cols), ~FALSE)) %>%
    dplyr::mutate(dplyr::across(dplyr::any_of(target_params), ~TRUE)) %>%
    dplyr::rowwise() %>%
    # Eliminate the target impute method from the impute column
    dplyr::mutate(dplyr::across(dplyr::any_of(impute_col), ~ paste0(setdiff(unlist(strsplit(.data[[impute_col]], ",")), target_impute),
                                                                    collapse = ","
    ))) %>%
    dplyr::mutate(dplyr::across(dplyr::any_of(impute_col), ~ ifelse(.data[[impute_col]] == "", NA_character_, .data[[impute_col]]))) %>%
    dplyr::ungroup() %>%
    # Make sure the class of the impute_col remains the same
    dplyr::mutate(dplyr::across(any_of(impute_col), ~ as.character(.data[[impute_col]]))) %>%
    as.data.frame()
  
  # Eliminate from the old intervals the target parameters
  old_intervals_with_impute <- target_intervals %>%
    dplyr::mutate(dplyr::across(dplyr::any_of(target_params), ~FALSE)) %>%
    dplyr::mutate(index = if (new_rows_after_original) index + 0.5 else index + max(index))
  
  # Make parameters FALSE in original intervals and join the new ones
  data <- data %>%
    dplyr::anti_join(target_intervals, by = names(data)) %>%
    dplyr::bind_rows(old_intervals_with_impute, new_intervals_without_impute) %>%
    dplyr::filter(rowSums(dplyr::across(dplyr::any_of(param_cols), as.numeric)) > 0) %>%
    dplyr::arrange(index) %>%
    dplyr::select(-index)
  
  data
}



