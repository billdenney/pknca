# PKNCA 0.7

* Features added
  * Additional PK parameters to support IV dosing added
  * Fixes #11, Intervals can be specified manually, and will apply across appropriate parts of the grouping variables
* Bugs fixed
  * Dosing without concentration is probably placebo; warn and continue

# PKNCA 0.6

First release targeting CRAN

Note that backward compatibility will not be guaranteed until version
1.0.  Argument and function changes will continue until then.  These
will be especially noticable around the inclusion of IV NCA parameters
and additional specifications of the dosing including dose amount and
route.
