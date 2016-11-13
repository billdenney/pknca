#' Merge lists of data to make a list of lists.
#'
#' @param \dots Lists of class \code{splitByData}.  If the arguments
#' are named, the names will be used for the output lists.
#' @param missing.value The value to place when there is not a
#' matching value.  If not given, set to NULL.
#' @return A list the length of merging all the \code{groupid}
#' attributes where each element is a list with the first element
#' coming from the first input list, the second from the second imput
#' list, and proceeding until all inputs are completed.  The order of
#' the output changes most slowly with the first and faster with each
#' subsequent argument.
#' @export
merge.splitByData <- function(..., missing.value=NULL) {
  args <- list(...)
  if (length(args) < 2)
    stop("Must give at least two lists to merge")
  ## Extract all the groupid values
  all.groupids <- lapply(args, FUN=function(x) attr(x, "groupid"))
  ## Get the maximum name to ensure that temporary column names for
  ## ordering the inputs don't interfere.
  if (any(sapply(all.groupids, is.null)))
    stop("All inputs must have a groupid attribute")
  col.order <-
    paste(max(do.call(c, lapply(X=all.groupids, FUN=names))),
          seq_along(args), sep=".")
  for (i in seq_len(length(args)))
    all.groupids[[i]][,col.order[i]] <- 1:nrow(all.groupids[[i]])
  intermediate.merge <- all.groupids[[1]]
  for (i in 2:length(args))
    intermediate.merge <- merge(intermediate.merge,
                                all.groupids[[i]],
                                all=TRUE)
  ## Sort the output so that it is in the order of the inputs (first
  ## arg changing most slowly).
  intermediate.merge <-
    intermediate.merge[do.call(order,
                               as.list(intermediate.merge[,col.order])),]
  ## Put together the names for all the items to go in the sub-lists.
  listnames <- 1:length(args)
  if (!is.null(names(args))) {
    mask.rename <- !(names(args) %in% "")
    listnames[mask.rename] <- names(args)[mask.rename]
  }
  ## Generate the output
  ret <- list()
  for (i in seq_len(nrow(intermediate.merge))) {
    ## For each row in the intermediate.merge (which is also each
    ## element of the outer list for the output), put a new list in as
    ## its main element.
    ret[[i]] <- list()
    for (j in seq_len(length(col.order))) {
      ## For each inner list, give it the name from the names of the
      ## input arguments and keeping them in order.  Then, pull the
      ## sub-list from the input arguments as chosen by the merging.
      list.idx <- intermediate.merge[i,col.order[j]]
      if (is.na(list.idx)) {
        ## Note the different number of brackets on the left side to
        ## allow assignment of list(NULL)
        ret[[i]][listnames[j]] <- list(missing.value)
      } else {
        ret[[i]][[listnames[j]]] <- 
          args[[j]][[intermediate.merge[i,col.order[j]]]]
      }
    }
  }
  ## Give back the data.
  attr(ret, "groupid") <- intermediate.merge[,
                                             setdiff(names(intermediate.merge), col.order),
                                             drop=FALSE]
  class(ret) <- c("mergedSplitByData", class(ret))
  ret
}

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
          finalgroupid, groups[[i]],
          by=intersect(names(finalgroupid), names(groups[[i]])))
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