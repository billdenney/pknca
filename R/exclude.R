#' Exclude data points or results from calculations or summarization.
#' 
#' @param object The object to exclude data from.
#' @param reason The reason to add as a reason for exclusion.
#' @param mask A logical vector or numeric index of values to exclude 
#'   (see details).
#' @param FUN A function to operate on the data (one group at a time) to
#'   select reasons for exclusions (see details).
#' @return The object with updated information in the exclude column. 
#'   The exclude column will contain the \code{reason} if \code{mask} or
#'   \code{FUN} indicate.  If a previous reason for exclusion was given,
#'   then subsequent reasons for exclusion will be added to the first 
#'   with a semicolon space ("; ") separator.
#'   
#' @details Only one of \code{mask} or \code{FUN} may be given.  If 
#'   \code{FUN} is given, it will be called with two arguments:  a 
#'   data.frame (or similar object) that consists of a single group of 
#'   the data and the full object (e.g. the PKNCAconc object), 
#'   \code{FUN(current_group, object)}, and it must return a logical 
#'   vector equivalent to \code{mask} or a character vector with the 
#'   reason text given when data should be excluded or 
#'   \code{NA_character_} when the data should be included (for the
#'   current exclusion test).
#' @examples
#' myconc <- PKNCAconc(data.frame(subject=1,
#'                                time=0:6,
#'                                conc=c(1, 2, 3, 2, 1, 0.5, 0.25)),
#'                     conc~time|subject)
#' exclude(myconc,
#'         reason="Carryover",
#'         mask=c(TRUE, rep(FALSE, 6)))
#' @export
#' @family Result exclusions
#' @importFrom dplyr "%>%"
#' @importFrom dplyr n
#' @importFrom rlang syms
exclude <- function(object, reason, mask, FUN)
  UseMethod("exclude")

utils::globalVariables(c("exclude_current_group_XXX", "row_number_XXX", "exclude_current_group_XXX_row_num"))
#' @describeIn exclude The general case for data exclusion
#' @export
exclude.default <- function(object, reason, mask, FUN) {
  dataname <- getDataName(object)
  # Check inputs
  if (missing(mask) & !missing(FUN)) {
    # operate on one group at a time
    groupnames <-
      unique(c(
        names(getGroups(object)),
        intersect(names(object[[dataname]]),
                  c("start", "end"))
      ))
    mask_df <-
      object[[dataname]] %>%
      dplyr::mutate(row_number_XXX=1:n()) %>%
      dplyr::group_by(!!! rlang::syms(groupnames)) %>%
      dplyr::mutate(
        exclude_current_group_XXX_row_num=row_number_XXX,
        exclude_current_group_XXX=
          do.call(
            FUN,
            list(
              as.data.frame(., stringsAsFactors=FALSE)[.$row_number_XXX %in% row_number_XXX,,drop=FALSE],
              object
            )
          )
      ) %>%
      dplyr::mutate(
        exclude_lengths_match=length(exclude_current_group_XXX) ==
          length(exclude_current_group_XXX_row_num)
      )
    # Extract the output and ensure that the output order equals the input order
    mask <- mask_df$exclude_current_group_XXX[order(mask_df$exclude_current_group_XXX_row_num)]
    if (is.character(mask)) {
      reason <- mask
      mask <- !is.na(reason)
    }
  } else if (!xor(missing(mask), missing(FUN))) {
    stop("Either mask for FUN must be given (but not both).")
  }
  if (!(length(reason) %in% c(1, nrow(object[[dataname]])))) {
    stop("reason must be a scalar or have the same length as the data.")
  } else if (!is.character(reason)) {
    stop("reason must be a character string.")
  }
  if (!("exclude" %in% names(object))) {
    stop("object must have an exclude column specified.")
  } else if (!(object$exclude %in% names(object[[dataname]]))) {
    stop("exclude column must exist in object[['", dataname, "']].")
  }
  # Make a scalar reason a vector
  if (length(reason) == 1)
    reason <- rep(reason, length(mask))
  # Find the original value of the 'exclude' column.
  orig <- object[[dataname]][[object$exclude]]
  if (length(mask) != length(orig)) {
    stop("mask must match the length of the data.")
  }
  # No current value for exclude
  mask.none <- orig %in% c(NA, "")
  # Replace the empty value with the reason
  mask.one <- mask & mask.none
  # Add the new reason to an existing reason
  mask.multiple <- mask & (!mask.one)
  ret <- orig
  if (any(mask.one)) {
    ret[mask.one] <- reason[mask.one]
  }
  if (any(mask.multiple)) {
    ret[mask.multiple] <- paste(ret[mask.multiple], reason[mask.multiple], sep="; ")
  }
  object[[dataname]][,object$exclude] <- ret
  object
}

#' Set the exclude parameter on an object
#' 
#' This function adds the exclude column to an object.  To change the 
#' exclude value, use the \code{\link{exclude}} function.
#' @param object The object to set the exclude column on.
#' @param exclude The column name to set as the exclude value.
#' @param dataname The name of the data.frame within the object to add
#'   the exclude column to.
#' @return The object with an exclude column and attribute
setExcludeColumn <- function(object, exclude, dataname="data") {
  add.exclude <- FALSE
  if (missing(exclude)) {
    # Exclude is not provided.
    if ("exclude" %in% names(object)) {
      # If exclude is already given, then do nothing.
    } else {
      add.exclude <- TRUE
    }
  } else if ("exclude" %in% names(object)) {
    # If exclude is already in the object, then make sure it matches
    # (and do nothing).
    if (!(object$exclude == exclude)) {
      stop("exclude is already set for the object.")
    }
  } else {
    # If exclude is not already in the object and it is given, then add
    # the column.
    add.exclude <- TRUE
  }
  if (add.exclude) {
    if (missing(exclude)) {
      # Generate the column name
      exclude <-
        setdiff(c("exclude", paste0("exclude.", max(names(object[[dataname]])))),
                names(object[[dataname]]))[1]
      object[[dataname]][,exclude] <- rep(NA_character_, nrow(object[[dataname]]))
    } else if (nrow(object[[dataname]]) == 0) {
      object[[dataname]][,exclude] <- rep(NA_character_, nrow(object[[dataname]]))
    } else if (!(exclude %in% names(object[[dataname]]))) {
      stop("exclude, if given, must be a column name in the input data.")
    } else {
      if (is.factor(object[[dataname]][,exclude])) {
        object[[dataname]][,exclude] <- as.character(object[[dataname]][,exclude])
      } else if (is.logical(object[[dataname]][,exclude]) &
                 all(is.na(object[[dataname]][,exclude]))) {
        object[[dataname]][,exclude] <- rep(NA_character_, nrow(object[[dataname]]))
      } else if (!is.character(object[[dataname]][,exclude])) {
        stop("exclude column must be character vector or something convertable to character without loss of information.")
      }
    }
    object[["exclude"]] <- exclude
  }
  object
}

#' Normalize the exclude column by setting blanks to NA
#'
#' @param object The object to extract the exclude column from
#' @return The exclude vector where \code{NA} indicates not to exclude
#'   and anything else indicates to exclude.
normalize_exclude <- function(object) {
  dataname <- getDataName(object)
  if (is.null(dataname)) {
    ret <- object
  } else {
    ret <- object[[dataname]][[object[["exclude"]]]]
  }
  mask_blank <- ret %in% ""
  if (any(mask_blank)) {
    ret[mask_blank] <- NA
  }
  ret
}
