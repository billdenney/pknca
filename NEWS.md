Note that backward compatibility of function arguments will not be
(near-)guaranteed until version 1.0.  Argument and function changes
will continue until then.  These will be especially noticable around
the inclusion of IV NCA parameters and additional specifications of
the dosing including dose amount and route.

# PKNCA 0.9.0

* BACKWARD INCOMPATIBILITY: The function supplied to the exclude
  argument 'FUN' now requires two arguments and operates on the level
  of a single group rather than the full object.  The function can
  also return the reason as a character string instead of a logical
  mask of when to exclude data.
* BACKWARD INCOMPATIBILITY: Added back-end functionality to only
  require one function to handle many NCA parameters that are related
  (e.g. combine pk.calc.aucpext, pk.calc.aucpext.obs,
  pk.calc.aucpext.pred, etc.).  If your current code calls a specific
  function (like pk.calc.aucpext.pred), you must change to using the
  generic function (like pk.calc.aucpext)
* BACKWARD INCOMPATIBILITY: Functions that previously may have
  returned Infinity due to dividing by zero (e.g. when AUC=0
  calculating clearance) now return NA.
* Add manual half-life point selection (Fixes issue #18)
* Add parameters for Ceoi and intravenous MRT
* Updated vignettes to improve clarity
* Added dose-normalized PK parameters (Fixes issue #41)
* Added checks to confirm that concentration and time are numeric
  (Fixes feature request #40)
* Improved test coverage

# PKNCA 0.8.1

* A PKNCAdose is no longer required for calculations.
* Data may now be excluded from calculations.

# PKNCA 0.8

This release is not backward compatible.  The switch to observed and
predicted-related NCA parameters (like aucinf.obs and aucinf.pred)
changed the format of the intervals specification.

* Remove dependency on doBy library
* Dose-aware interpolation and extrapolation was added with the interp.extrap.conc.dose function.
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
