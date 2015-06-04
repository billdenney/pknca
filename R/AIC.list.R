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
AIC.list <- function(object, ..., assess.best=TRUE,
                     reference.model=get.first.model(object)) {
  ret <- data.frame()
  for (i in 1:length(object)) {
    if (identical(object[[i]], NA)) {
      tmp <- data.frame(df=NA, AIC=NA, indentation=0)
      rownames(tmp) <- names(object)[i]
      ret <- rbind(ret, tmp)
    } else {
      tmp <- AIC(object[[i]], ..., reference.model=reference.model)
      if ("indentation" %in% names(tmp)) {
        ## If the object was another list, then add to the indentation
        ## level
        tmp$indentation <- tmp$indentation + 1
      } else {
        ## Otherwise, this is the first time that the indentation is
        ## being assessed.
        tmp$indentation <- 0
      }
      ret <- rbind(ret, tmp[-1,])
    }
  }
  if (assess.best) {
    ret$isBest <- ""
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
  object[AIC(object, ...)$isBest %in% "Best Model"][[1]]

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
