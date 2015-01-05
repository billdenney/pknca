#' Merge lists of data to make a list of lists.
#'
#' @param \dots Lists of class \code{splitByData}
#' @param missing.value The value to place when there is not a
#' matching value.  If not given, set to NULL.
#' @return A list the length of merging all the \code{groupid}
#' attributes where each element is a list with the first element
#' coming from the first input list, the second from the second imput
#' list, and proceeding until all inputs are completed.  Note that the
#' order of the output will be defined by the order returned from
#' merging the \code{groupid} attributes in order.
#' @export
merge.splitByData <- function(..., missing.value) {
  args <- list(...)
  if (length(args) < 2)
    stop("Must give at least two lists to merge")
  ## Extract all the groupid values
  all.groupids <- lapply(args, FUN=function(x) attr(x, "groupid"))
  ## Get the maximum name to ensure that temporary column names for
  ## ordering the inputs don't interfere.
  col.order <-
    paste(max(sapply(all.groupids, names, USE.NAMES=FALSE)),
          1:length(args), sep=".")
  for (i in seq_len(length(args)))
    all.groupids[[i]][,col.order[i]] <- 1:nrow(all.groupids[[i]])
  intermediate.merge <- all.groupids[[1]]
  for (i in 2:length(args))
    intermediate.merge <- merge(intermediate.merge,
                                all.groupids[[i]],
                                all=TRUE)
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
    for (j in length(col.order)) {
      ## For each inner list, give it the name from the names of the
      ## input arguments and keeping them in order.  Then, pull the
      ## sub-list from the input arguments as chosen by the merging.
      ret[[i]][[listnames[j]]] <-
        args[[j]][[intermediate.merge[i,col.order[j]]]]
    }
  }
  ## Give back the data.
  ret
}
