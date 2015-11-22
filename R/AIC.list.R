#' Assess the AIC for all models in a list of models
#'
#' @param object the list of models
#' @param \dots parameters passed to the underlying AIC function
#' (typically the parameter k)
#' @param assess.best determine which model is the best (by lowest
#' AIC)
#' @seealso \code{\link{get.best.model}}
#' @return a data frame with row names matching the names of the list
#' \code{x} and columns for degrees of freedom (\code{df}) and
#' \code{AIC}.  If \code{assess.best} is true, then there will be
#' another column \code{isBest}.
#' @export
AIC.list <- function(object, ..., assess.best=TRUE) {
  allAICs <-
    lapply(object, FUN=function(subobject, ...) {
      ## Return the AIC of the new model relative to the reference model
      if (identical(NA, subobject)) {
        ret <- data.frame(AIC=NA, df=NA, indentation=0)
      } else {
        ret <- stats::AIC(subobject, ...)
        if (is.numeric(ret)) {
          ret <- data.frame(AIC=ret,
                            df=attr(stats::logLik(subobject), "df"),
                            indentation=0)
        } else if (is.data.frame(ret)) {
          if ("indentation" %in% names(ret)) {
            ret$indentation <- ret$indentation + 1
          } else {
            stop("Unknown way to get a data.frame without indentation set.  This is likely a bug.")
          }
        }
      }
      ret
    })
  retnames <- names(allAICs)
  if (is.null(retnames))
    retnames <- rep("", length(allAICs))
  ret <- data.frame()
  for (i in seq_len(length(allAICs))) {
    tmpAICs <- allAICs[[i]]
    ## If the best value has already been established, drop it for
    ## assessment later.
    tmpAICs$isBest <- NULL
    ## Assign the corret rownames to tmpAICs
    if (!(retnames[i] %in% ""))
      if (nrow(tmpAICs) > 1 |
          !identical(rownames(tmpAICs), as.character(1:nrow(tmpAICs)))) {
        rownames(tmpAICs) <- paste(retnames[i], rownames(tmpAICs))
      } else {
        rownames(tmpAICs) <- retnames[i]
      }
    ## Add tmpAICs to the data frame to return
    ret <- rbind(ret, tmpAICs)
  }
  if (assess.best) {
    ret$isBest <- ""
    ## The next row prevents warnings about no data when na.rm=TRUE
    if (any(!is.na(ret$AIC)))
      ret$isBest[ret$AIC %in% min(ret$AIC, na.rm=TRUE)] <- "Best Model"
  }
  ret
}

#' Extract the best model from a list of models using AIC.list.
#'
#' @param object the list of models
#' @param \dots Parameters passed to AIC.list
#' @seealso \code{\link{AIC.list}}
#' @return The model which is assessed as best.  If more than one are
#' equal, the first is chosen.
#' @export
get.best.model <- function(object, ...)
  object[stats::AIC(object, ...)$isBest %in% "Best Model"][[1]]

#' Get the first model from a list of models
#'
#' @param object the list of (lists of, ...) models
#' @return The first item in the \code{object} that is not a list or
#' \code{NA}.  If \code{NA} is passed in or the list (of lists) is all
#' \code{NA}, then \code{NA} is returned.
get.first.model <- function(object) {
  ret <- NA
  if (inherits(object, "list")) {
    idx <- 0
    while (identical(NA, ret) & idx < length(object)) {
      idx <- idx + 1
      if (identical(NA, object[[idx]])) {
        ## Do nothing
      } else if (inherits(object[[idx]], "list")) {
        ret <- get.first.model(object[[idx]])
      } else {
        ## It is neither NA or a list, it's our first usable object;
        ## return it.
        ret <- object[[idx]]
      }
    }
  } else {
    ret <- object
  }
  ret
}
