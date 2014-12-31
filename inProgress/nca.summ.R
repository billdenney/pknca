############################################################
##  This Script is for Summarizing the NCA Data
############################################################

## d.ind.nca <- read.csv(filename.output[["indiv.nca"]])

## Find various column names that will be used for removing values due
## to business rules.

## AUC
auc.cols <- grep("^AUC[0-9]", names(d.ind.nca), value=TRUE)
## AUMC
aumc.cols <- grep("^AUMC[0-9]", names(d.ind.nca), value=TRUE)
## Concentration
conc.cols <- grep("^Conc.[0-9]", names(d.ind.nca), value=TRUE)
## AUC*Inf
aucinf.cols <- grep("^AUM?C.*Inf.*$", names(d.ind.nca), value=TRUE)

## Find various column names that will be used for standard reporting
## and list them by averaging type.
report.auc.cols <-
  grep("^AUC[0-9]+_([0-9]+|last|Inf)$", names(d.ind.nca), value=TRUE)
report.conc.cols <-
  grep("^Conc.[0-9]+$", names(d.ind.nca), value=TRUE)

## Functions use for summarization in the report table
report.method.am <- function(param, data) {
  base <- data[,paste(param, "am", sep=".")]
  ret <- sprintf("%s (%s)",
                 signif(base, 3),
                 signif(data[,paste(param, "am.cv", sep=".")], 3))
  ret[is.na(base)] <- "NC"
  ret
}
report.method.gm <- function(param, data) {
  base <- data[,paste(param, "gm", sep=".")]
  ret <- sprintf("%s (%s)",
                 signif(base, 3),
                 signif(data[,paste(param, "gm.cv", sep=".")], 3))
  ret[is.na(base)] <- "NC"
  ret
}
report.method.median <- function(param, data) {
  base <- data[,paste(param, "median", sep=".")]
  ret <- sprintf("%s (%s-%s)",
                 signif(base, 3),
                 signif(data[,paste(param, "min", sep=".")], 3),
                 signif(data[,paste(param, "max", sep=".")], 3))
  ret[is.na(base)] <- "NC"
  ret
}
report.method.count <- function(param, data) {
  base <- data[,paste(param, "n.tot", sep=".")]
  ret <- sprintf("%d", base)
  ret[is.na(base)] <- "NC"
  ret
}

## Setup all estimation methods for the table and the column order
report <-
  data.frame(param=c("CMAX", "TMAX",
               report.auc.cols,
               report.conc.cols,
               "THALF_T",
               "CMAX"),
             method=c("gm", "median",
               rep("gm", length(report.auc.cols)),
               rep("gm", length(report.conc.cols)),
               "am",
               "count"), stringsAsFactors=FALSE)
## Make sure that the report ends up in the right order
report$order <- 1:nrow(report)
## Only report on columns that are given in the d.ind.nca data frame.
report <- subset(report, param %in% names(d.ind.nca))
## Make column headers look nicer
report$pretty.name <- report$param
report$pretty.name[report$param %in% "CMAX"] <- "Cmax"
report$pretty.name[report$param %in% "TMAX"] <- "Tmax"
report$pretty.name[report$param %in% "THALF_T"] <- "t1/2"
## Pretty names for methods
method.frame <-
  data.frame(method=c("am", "gm", "median", "count"),
             pretty.method=c("[Mean (%CV)]",
               "[Geo Mean (%CV)]",
               "[median (range)]",
               "[N]"))
report <- merge(method.frame, report, all=TRUE)
report$pretty.name <- with(report, paste(pretty.name, pretty.method))
report <- report[order(report$order),]

############################################################
## Apply business rules for what data to include.
## Is the AUCInf extrapolating too far?
d.ind.nca$Exclusion.Reason <- ""
half.life.dependent.cols <-
  c("LAMBDA.Z", "THALF_T", "TFIRST.LZ", "N.POINTS.LZ", "YINT",
    "r.squared", "adj.r.squared", "AUCPEXT", "AUMCPEXT", "MRT", "KEL",
    "THALF_EFF", "CL", "VZ", "VSS", aucinf.cols)
if ("AUCPEXT" %in% names(d.ind.nca)) {
  mask.excessive.extrapolation <- (!is.na(d.ind.nca$AUCPEXT) &
                                   d.ind.nca$AUCPEXT >= pds.pk$max.AUCInf.pext)
  d.ind.nca[mask.excessive.extrapolation, aucinf.cols] <- NA
  d.ind.nca$Exclusion.Reason[mask.excessive.extrapolation] <-
    sprintf("AUC %% Extrapolated > %g", pds.pk$max.AUCInf.pext)
}
## Is the span ratio sufficient?
if ("SPANR" %in% names(d.ind.nca)) {
  mask.insufficient.span <- (!is.na(d.ind.nca$SPANR) &
                             d.ind.nca$SPANR < pds.pk$min.span.ratio)
  d.ind.nca[mask.insufficient.span, half.life.dependent.cols] <- NA
  d.ind.nca$Exclusion.Reason[mask.insufficient.span] <-
    sprintf("Span ratio < %g", pds.pk$min.span.ratio)
}
## Is the r-squared value sufficient?
if ("r.squared" %in% names(d.ind.nca)) {
  mask.insufficient.rsquared <-
    (!is.na(d.ind.nca$r.squared) &
     d.ind.nca$r.squared < pds.pk$min.hl.r.squared)
  d.ind.nca[mask.insufficient.rsquared, half.life.dependent.cols] <- NA
  d.ind.nca$Exclusion.Reason[mask.insufficient.rsquared] <-
    sprintf("r-squared value < %g", pds.pk$min.hl.r.squared)
}

## Remove subjects without CMAX because they have not had any samples
## analyzed.
d.ind.nca <- subset(d.ind.nca, !is.na(CMAX))

############################################################
## Do the summarization

## What parameters are used for sorting the output into different
## rows?
SORTV <- c("TRTG", "DAY", "NTAFD.END")
SORTV <- intersect(SORTV, names(d.ind.nca))

## What parameters do we need to summarize?  Summarize everything
## imaginable so that the data are available.  The report will only
## have items listed in the report above, but all data should be
## available in case it is desired to examine.
text.lhs <-
  c("TMAX", "CMAX", "CMAX_D", "CMIN", "CMIN_D", "CLAST",
    ## All "Concentration at time" measures
    conc.cols,
    ## All AUCs
    auc.cols,
    ## All AUMCs
    aumc.cols,
    "LAMBDA.Z", "r.squared", "adj.r.squared", "SPANR",
    "THALF_T", "MRT", "KEL", "THALF_EFF", "CL", "VZ",
    "VSS")
text.lhs <- paste(text.lhs[text.lhs %in% names(d.ind.nca)], collapse="+")
text.rhs <- paste(SORTV, collapse="+")

## A function to summarize all measures of interest with protections
## to prevent large numbers of errors when all data are missing.
detailed.summary <- function(x, ...) {
  c(min=business.min(x, ...),
    median=business.median(x, ...),
    max=business.max(x, ...),
    am=business.mean(x, ...),
    sd=business.sd(x, ...),
    am.cv=business.cv(x, ...),
    gm=business.geomean(x, ...),
    gm.cv=business.geocv(x, ...),
    n.tot=length(x),
    n.missing=sum(is.na(x)))
}

summaryfun <-
  parse(text=sprintf("summaryBy(%s~%s, data=d.ind.nca, FUN=detailed.summary, na.rm=TRUE)",
          text.lhs, text.rhs))
pk.summary.values <- eval(summaryfun)
pk.summary.values$TRTG <- as.factor.smarttext(pk.summary.values$TRTG)

## Generate the output PK summary table
RSUM <- pk.summary.values[,SORTV,drop=FALSE]
for (i in 1:nrow(report)) {
  ## For each of the desired columns,
  RSUM[,report$pretty.name[i]] <-
    do.call(paste("report.method", report$method[i], sep="."),
            list(report$param[i], pk.summary.values))
}
RSUM <- RSUM[order(RSUM$TRTG),]
attr(RSUM, "caption") <-
  paste("NC - not calculated due to business rules;",
        "All values reported as BLQ have been replaced with zero for all calculations.",
        sep="\n")

## Write the PK table out
output.data[["unblinded"]][["pk"]] <- RSUM
## Write the output csvs.
output.data[["csv"]][["nca.parsum"]] <- pk.summary.values
## Writing the cleaned individual NCA results (after applying the
## business rules)
output.data[["csv"]][["individual.nca"]] <- d.ind.nca
