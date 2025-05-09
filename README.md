<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/PKNCA)](https://CRAN.R-project.org/package=PKNCA)
[![R-CMD-check](https://github.com/humanpred/pknca/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/humanpred/pknca/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/humanpred/pknca/graph/badge.svg)](https://app.codecov.io/gh/humanpred/pknca)
[![CodeFactor](https://www.codefactor.io/repository/github/humanpred/pknca/badge?style=plastic)](https://www.codefactor.io/repository/github/humanpred/pknca)
[![CII Best Practices](https://bestpractices.coreinfrastructure.org/projects/2054/badge)](https://bestpractices.coreinfrastructure.org/projects/2054)
<!-- badges: end -->

The PKNCA R Package <img src="man/figures/pknca-hex.png" align="right" />
=====

The PKNCA R package is designed to perform all noncompartmental
analysis (NCA) calculations for pharmacokinetic (PK) data.  The
package is broadly separated into two parts (calculation and summary)
with some additional housekeeping functions.

The primary and secondary goals of the PKNCA package are to 1) only
give correct answers to the specific questions being asked and 2)
automate as much as possible to simplify the task of the analyst. When
automation would leave ambiguity or make a choice that the analyst may
have an alternate preference for, it is either not used, is possible
to override or raises an error message.

Note that backward compatibility will not be guaranteed until version
1.0.  Argument and function changes will continue until then.  These
will be especially noticeable around the inclusion of IV NCA parameters and additional specifications of the drug, which include dose amount and route of administration.

# Citation

Citation information for the PKNCA package is available with a call to
`citation(package="PKNCA")`.  The preferred citation until publication
of version 1.0 is below:

Denney W, Duvvuri S and Buckeridge C (2015). "Simple, Automatic
Noncompartmental Analysis: The PKNCA R Package." _Journal of
Pharmacokinetics and Pharmacodynamics_, *42*(1), pp. 11-107,S65. ISSN
1573-8744, doi: 10.1007/s10928-015-9432-2, <URL:
https://github.com/humanpred/pknca>.

# Installation

## From CRAN

The current stable version of PKNCA is available on CRAN.  You can
install it and its dependencies using the following command:

    install.packages("PKNCA")

## From GitHub

To install the development version from GitHub, type the following commands:

    install.packages("remotes")
    remotes::install_github("humanpred/pknca")

# Calculating parameters

    # Load the package
    library(PKNCA)
    # Set the business rule options with the PKNCA.options() function
    # Load your concentration-time data
    conc_raw <- read.csv("myconc.csv", stringsAsFactors=FALSE)
    # Load your dose data
    dose_raw <- read.csv("mydose.csv", stringsAsFactors=FALSE)
    # Put your concentration data into a PKNCAconc object
    o_conc <- PKNCAconc(data=conc_raw,
                        formula=conc~time|treatment+subject/analyte)
    # Put your dose data into a PKNCAdose object
    o_dose <- PKNCAdose(data=dose_raw,
                        formula=dose~time|treatment+subject)
    # Combine the two (and automatically determine the intervals of
    # interest
    o_data <- PKNCAdata(o_conc, o_dose)
    # Compute the NCA parameters
    o_results <- pk.nca(o_data)
    # Summarize the results
    summary(o_results)

More help is available in the function help files. Be sure to look at the PKNCA.options function for choices on making PKNCA conform to your company’s business rules on calculation and summarization.

# Feature requests

Please use the github issues page
(https://github.com/humanpred/pknca/issues) to make feature requests
and bug reports.
