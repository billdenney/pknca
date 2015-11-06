#' Compute noncompartmental pharmacokinetics
#'
#' PKNCA computes pharmacokinetic (PK) noncompartmental analysis (NCA)
#' parameters.
#'
#' Typically, you will start with the \code{\link{pk.nca}}
#' function and then use the \code{\link{summary}} function to
#' generate tables for reporting.
#'
#' @docType package
#' @name PKNCA
#' @importFrom doBy renameCol splitBy summaryBy
#' @importFrom lattice xyplot
#' @importFrom nlme fixef getData getGroups gnls intervals lme nlme ranef
#' @importFrom parallel mclapply
#' @importFrom plyr rbind.fill
#' @importFrom digest digest
NULL
