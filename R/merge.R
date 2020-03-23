#' Merge two or more lists with a data.frame 'groupid' attribute 
#' defining the matching.
#' 
#' @param ... lists with 'groupid' attributes
#' @return A list of lists with elements for the matching items between 
#'   the 'groupid's of each of the input lists.  The output list will 
#'   have a new 'groupid' attribute added with additional columns to 
#'   indicate the of each input to its output location.  If the inputs 
#'   are named, then each list item will be named the same as the input
#'   name.
#' @details The merge is equivalent to a \code{full_join} where items 
#'   missing from one or the other item will be missing, but the 
#'   element(s) will exist.
#' @export
merge.splitlist <- function(...) {
  args <- list(...)
  # Check inputs
  if (length(args) < 2) {
    stop("At least two lists must be given")
  }
  if (any(!sapply(args, is.list))) {
    stop("All arguments must be lists")
  }
  groups <- lapply(args, function(x) attr(x, "groupid"))
  if (any(sapply(groups, is.null))) {
    stop("All arguments must have a 'groupid' attribute")
  }
  if (!all(sapply(groups, nrow) == sapply(args, length))) {
    stop("The number of rows in the 'groupid' attribute must match the length of the list.")
  }
  # Prepare the names for the columns that will identify the index from each list.
  idxname <- max(unlist(lapply(groups, names)))
  idxname <- paste0(idxname, seq_along(args))
  for (i in seq_along(groups)) {
    groups[[i]][[idxname[i]]] <- seq_len(nrow(groups[[i]]))
    if (i == 1) {
      finalgroupid <- groups[[i]]
    } else {
      finalgroupid <-
        dplyr::full_join(
          finalgroupid,
          groups[[i]],
          by=intersect(names(finalgroupid), names(groups[[i]]))
        )
    }
  }
  # Combine elements
  ret <- list()
  for (i in seq_len(nrow(finalgroupid))) {
    ret[[i]] <- list()
    for (j in seq_along(idxname)) {
      curname <- idxname[j]
      if (is.na(finalgroupid[[curname]][i])) {
        ret[[i]][j] <- list(NULL)
      } else {
        ret[[i]][[j]] <- args[[j]][[finalgroupid[[curname]][i]]]
      }
    }
    names(ret[[i]]) <- names(args)
  }
  attr(ret, 'groupid') <- finalgroupid[,setdiff(names(finalgroupid), idxname), drop=FALSE]
  attr(ret, 'sourcemap') <- finalgroupid[,idxname, drop=FALSE]
  names(attr(ret, 'sourcemap')) <- names(args)
  ret
}