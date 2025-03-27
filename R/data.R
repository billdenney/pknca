#' SDTM mapping of parameter to its code and description
#'
#' @format
#' \describe{
#'   \item{PPTESTCD}{The code describing the PK parameter}
#'   \item{PPTEST}{The human-readable description of the PK parameter}
#'   \item{CDISC_definition}{The long description of the parameter}
#' }
#' @source <https://evs.nci.nih.gov/ftp1/CDISC/SDTM/SDTM%20Terminology.xls>
"sdtm_paramcd"

#' SDTM mapping of parameter to its code and description
#'
#' @format
#' \describe{
#'   \item{submission}{The preferred unit for PK parameters, may be duplicated if there are multiple synonyms}
#'   \item{synonym}{Synonyms for the preferred unit, may be `NA` if there is no synonym}
#' }
#' @source <https://evs.nci.nih.gov/ftp1/CDISC/SDTM/SDTM%20Terminology.xls>
"sdtm_pkunit"
