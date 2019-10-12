#' Add a hash and associated information to enable checking object 
#' provenance.
#' 
#' @param object The object to add provenance
#' @param replace Replace provenance if the object already has a 
#'   provenance attribute.  (If the object already has provenance and
#'   \code{replace} is \code{FALSE}, then an error will be raised.)
#' @return The object with provenance as an added item
#' @seealso \code{\link{checkProvenance}}
#' @export
#' @importFrom digest digest
#' @importFrom utils sessionInfo
addProvenance <- function(object, replace=FALSE) {
  if (replace) {
    attr(object, "provenance") <- NULL
  }
  if (is.null(attr(object, "provenance", exact=TRUE))) {
    # Get most of the provenance added
    tmp.prov <- list(
      sessionInfo=utils::sessionInfo(),
      datetime=Sys.time(),
      sysInfo=Sys.info(),
      hash=NA)
    class(tmp.prov) <- c("provenance", class(tmp.prov))
    attr(object, "provenance") <- tmp.prov
    attr(object, "provenance")$hash <-
      digest::digest(as.character(object), serialize=FALSE)
  } else {
    stop("object already has provenance and the option to replace it was not selected.")
  }
  object
}

#' Check the hash of an object to confirm its provenance.
#' 
#' @param object The object to check provenance for
#' @return \code{TRUE} if the provenance is confirmed to be consistent, 
#'   \code{FALSE} if the provenance is not consistent, or \code{NA} if
#'   provenance is not present.
#' @seealso \code{\link{addProvenance}}
#' @export
checkProvenance <- function(object) {
  tmp.prov <- attr(object, "provenance", exact=TRUE)
  if (is.null(tmp.prov)) {
    NA
  } else {
    hash <- tmp.prov$hash
    attr(object, "provenance")$hash <- NA
    (hash == digest::digest(as.character(object), serialize=FALSE))
  }
}

#' Print the summary of a provenance object
#' 
#' @param x The object to be printed
#' @param ... Ignored
#' @return invisible text of the printed information
#' @export
print.provenance <- function(x, ...) {
  ret <- sprintf("Provenance hash %s generated on %s with %s.",
                 x$hash,
                 x$datetime,
                 x$sessionInfo$R.version[["version.string"]])
  cat(ret, "\n", sep="")
  invisible(ret)
}
