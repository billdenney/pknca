Design of the PKNCA R Package
=====

The PKNCA R package is designed to perform all noncompartmental
analysis (NCA) calculations for pharmacokinetic (PK) data.  The
package is broadly separated into two parts (calculation and summary)
with some additional housekeeping functions.

The primary and secondary goals of the PKNCA package are to 1) only
give correct answers to the specific questions being asked and 2)
automate as much as possible to simplify the task of the analyst. When
automation would leave ambiguity or make a choice that the analyst may
have an alternate preference for, it is either not used or is possible
to override.

# Function Groups

## Housekeeping and Option Setting

So that each session can be simplified, many of the options passed to
a function on calculation or summary choices have their defaults set
in the `PKNCA.options` function.  Any option that is set in
`PKNCA.options` can be overridden on an individual function call, but
consistency is maximized by setting them once.

Two data cleaning functions handle missing (concentration=NA) and
below limit of quantification (BLQ, concentration=0) concentrations.
These enable 0 and NA to be handled equivalently across all
calculations and summaries.

# Calculation Functions

Each calculation function is a generic function that can accept
numeric vectors, a data frame and a grouped formula, or groupedData
(from the nlme package).  Each of the functions that accepts numeric
vectors assumes that it is for a single subject in a single interval
for summary (e.g. the pk.calc.tmax function for numeric vectors
assumes that you want the Tmax of all measurements passed in, and it
does not subset the data further).

The groupedData and grouped formula with data frame methods separate
the subjects based on the groups specified and then passes the data
into the numeric vector methods.  The grouped formulae are written as
`concentration~time|groups`.  The groups are specified where a plus
sign (`+`) separates parameters between subjects and a front slash
(`/`) separates parameters within a subject.  The parameters within a
subject should always start with the subject identifier column.  For
example, the groups may be written as `study+subject/analyte/day`.

Dosing data is specified equivalently to `concentration~time` data
with an addition if there are multiple analytes and/or drugs given.
With multiple analytes and dosing of more than one drug (e.g. in a
drug interaction study), a data frame mapping the drug name to the
analyte name must be given if the names are not identical.  An example
of this when running a midazolam drug interaction study where
midazolam and 5-OH midazolam are measured with DRUG1 can be
data.frame(drug=c("midazolam", "midazolam", "DRUG1"),
analyte=c("midazolam", "5-OH midazolam", "DRUG1")).

The numeric vector methods return a number or data frame, as
appropriate.  The other methods return a data frame that has the
additional class of PKNCAresults to assist with summaries.

# Summary Functions

The summary functions take input of parameters to summarize on the
left hand side (assuming all if none are specified), and the right
hand side indicates which groups to remove (with minus signs) or which
groups to keep (with plus signs).  Two equivalent summaries from a
calculation input of "concentration~time|study+subject/analyte/day"
are "~-subject" and "~study+analytes/day", or to just summarize the
AUC0_24, write "AUC0_24~-subject".  The right hand side for removal
and keeping groups cannot be mixed.  If no removal or addition is
specified, then the output is designed to be appropriate for
generation of listings of all parameters as calculated; keeping
everything can be written as ".~.".

# Calculation Details

## Cmax

## Cmin

## Clast observed

## Clast predicted

## Tmax

## Tlast

## Tfirst

## Half-Life

## Area Under the Curve (AUC)

### Automatically Choosing AUC Intervals

### AUC with a Specific Interval

### AUCall

### AUClast

### AUCinf

### AUC percent extrapolated

## Observed Systemic Clearance (CL)

## Mean Residence Time (MRT)

## Volume of Distribution

### Terminal Volume of Distribution (Vz)

### Steady-State Volume of Distribution (Vss)

## Time to Steady-State (TSS)

### Monoexponential Time to Steady-State

### Stepwise-Linear Time to Steady-State

## Concentration Interpolation

# Summary Details

# Business Rules
