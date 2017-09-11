#' Compute noncompartmental pharmacokinetics
#' 
#' Compute pharmacokinetic (PK) noncompartmental analysis (NCA) 
#' parameters.
#' 
#' A common workflow would load data from a file or database into a 
#' data.frame then run the following
#' @examples
#' \dontrun{
#' # Load concentration-time data into a data.frame called d.conc
#' # with columns named "conc", "time", and "subject".
#' my.conc <- PKNCAconc(d.conc, conc~time|subject)
#' # Load dose-time data into a data.frame called d.dose
#' # with columns named "dose", "time", and "subject".
#' my.dose <- PKNCAdose(d.dose, dose~time|subject)
#' # Combine the concentration-time and dose-time data into an object
#' # ready for calculations.
#' my.data <- PKNCAdata(my.conc, my.dose)
#' # Perform the calculations
#' my.results <- pk.nca(my.data)
#' # Look at summary results
#' summary(my.results)
#' # Look at a listing of results
#' as.data.frame(my.results)
#' }
#' @docType package
#' @name PKNCA
#' @importFrom digest digest
#' @importFrom graphics plot
#' @importFrom lattice xyplot
#' @importFrom nlme fixef getData gnls intervals lme nlme ranef
#' @importFrom parallel mclapply
#' @importFrom plyr rbind.fill
#' @importFrom stats AIC as.formula coef confint formula glm logLik lm 
#'   median model.frame na.exclude na.omit predict sd update 
#'   update.formula
#' @importFrom utils head sessionInfo
NULL

# To work with the use of dplyr's pipe within the exclude function
utils::globalVariables(".")
