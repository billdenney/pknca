Note that backward compatibility will not be guaranteed until version
1.0.  Argument and function changes will continue until then.  These
will be especially noticable around the inclusion of IV NCA parameters
and additional specifications of the dosing including dose amount and
route.

# PKNCA 0.8

This release is not backward compatible.  The switch to observed and
predicted-related NCA parameters (like aucinf.obs and aucinf.pred)
changed the format of the intervals specification.

* Added Clast.pred related NCA calculations
* Added N to summary of PKNCAresults
* Added parameter selection between Clast,observed and Clast,predicted across all parameters
* Enabled PKNCAdose to be specified with one-sided formula
* Improved error reporting so that the group and time (interval specification) is reported in addition to the error.
* PKNCAdose now allows route of administration and IV infusion parameters of rate/duration to be specified

# PKNCA 0.7.1

* Updated vignettes
* Standardize rounding and significance with missing values in signifString and roundString
* Enable wide data output with as.data.frame(PKNCAresults, out.format="wide")
* Correct calculation of Vz
* Various CRAN-related cleanups

# PKNCA 0.7

* Features added
  * Additional PK parameters to support IV dosing added
  * Fix #11, Intervals can be specified manually, and will apply across appropriate parts of the grouping variables
  * Enable dose and dose.time as parameters to NCA calculations
  * More NCA parameters are calculated, especially related to IV dosing
  * Fix #8, Reporting times for time-based parameters are now within the current interval rather than since first dose (e.g. Tmax on day 14 should be between 0 and 24 not 2*7*24+c(0, 24))
  * Added several vignettes
* Bugs fixed
  * Dosing without concentration is probably placebo; warn and continue
  * Fix #6, make merge.splitByData work with more than one dosing level
  * Export some generic classes that were not previously exported to simplify their use
  * Superposition extensions when lambda.z cannot be calculated
  * Significance rounding into character strings works when the rounding moves up one order of magnitude.
  * Fix #9, summarization of parameters that are not calculated show not calculated instead of missing.

# PKNCA 0.6

First release targeting CRAN
