## Compute the max allowable dose based on PK stopping criteria

## Get the other types that have diversity in values
simplify.other.types <- function(other.types, data) {
  if (length(other.types) > 0) {
    ## Only keep a column if there is diversity in its values.
    mask.keep <- sapply(data[,other.types, drop=FALSE],
                        FUN=function(x) {length(unique(x)) != 1})
    other.types <- other.types[mask.keep]
  }
  other.types
}

get.other.types.txt <- function(other.types) {
  if (length(other.types) > 0) {
    other.types.txt <-
      paste(c("", sprintf("factor(%s)", other.types)),
            collapse="+")
  } else {
    other.types.txt <- ""
  }
  other.types.txt
}

get.other.types.txt.plt <- function(other.types) {
  if (length(other.types) > 0) {
    other.types.txt.plt <-
      paste(c("", sprintf("factor(%s)", other.types)),
            collapse="|")
  } else {
    other.types.txt.plt <- ""
  }
  other.types.txt.plt
}


##@@@@@@@@@@@@@@@@@@@@@@@@@@
## Create space for new plots for simple comparison
output.data[["unblinded"]][["PK Stopping Limit"]] <- list()
##@@@@@@@@@@@@@@@@@@@@@@@@@@

if (identical(d.pk.stopping, NA)) {
  print("No PK stopping rules given, no maximum dose prediction made.")
} else {
  ## Merge the randomization code data into the individual NCA data.
  d.ind.nca <- my.merge(d.ind.nca, randcodes)
  output.data[["blinded"]][["PK.Stopping"]] <- d.pk.stopping

  ## Find potential alternate treatment methods (other types of
  ## dosing) that may not be included.
  other.types <-
    simplify.other.types(c("DAY", setdiff(names(d.trtcd),
                                          c("TRT.CDE", "TRTG","DOSE"))),
                         d.ind.nca)
  ## Setup the formula portion for the other types and the prediction
  ## table for the other types of dosing.
  if (length(other.types) > 0) {
    for (n in other.types) {
      output.data[["blinded"]][["PK.Stopping"]] <-
        merge(output.data[["blinded"]][["PK.Stopping"]],
              unique(d.ind.nca[,n,drop=FALSE]))
    }
  }

  ## Determine the number of observations per subject
  n.sub <- length(unique(d.ind.nca$RANDID))
  n.obs <- nrow(unique(d.ind.nca[,c("RANDID", "TRTG", other.types)]))

  for (s in names(d.pk.stopping)) {
    if (!s %in% names(d.ind.nca)) {
      warning(sprintf("Cannot estimate maximum allowed dose for %s, the parameter was not computed in the individual nca parameters.  See output individual parameter columns for allowed values or add the value to your AUC intervals.", s))
    } else {
      d.tmp.pk.limit <- d.ind.nca
      ## Remove rows with NA in the prediction column
      d.tmp.pk.limit <- d.tmp.pk.limit[!is.na(d.tmp.pk.limit[,s]),]
      ## If there's only one dose level, the prediction will be
      ## singular.  Add a row of placebo is 0 response for any of the
      ## other types of dosing (if applicable) to assume a linear effect.
      model.type <- "exp"
      if (length(unique(d.tmp.pk.limit$DOSE)) == 1) {
        model.type <- "linear"
        new.rows <-
          (nrow(d.tmp.pk.limit)+1):(nrow(d.tmp.pk.limit)+nrow(output.data[["blinded"]][["PK.Stopping"]]))
        d.tmp.pk.limit[new.rows,] <- NA
        if (length(other.types) > 0) {
          ## This may not work...
          d.tmp.pk.limit[new.rows,other.types] <-
            output.data[["blinded"]][["PK.Stopping"]][,other.types]
        }
        d.tmp.pk.limit[new.rows,c("DOSE", s)] <- 0
        ## Set the RANDID for these new rows to some ficticious value.
        if (is.factor(d.tmp.pk.limit$RANDID)) {
          ## It always should be, but it's better to test.
          new.randid <- min(as.numeric(as.character(levels(d.tmp.pk.limit$RANDID))),
                            na.rm=TRUE)
          levels(d.tmp.pk.limit$RANDID) <-
            c(levels(d.tmp.pk.limit$RANDID), new.randid)
        } else {
          ## Assuming that this is not used
          new.randid <- -99999
        }
        d.tmp.pk.limit$RANDID[new.rows] <- new.randid
      }

      ## Note that this formula looks like it's predicting the opposite
      ## way than we want, but we're wanting to predict the dose that
      ## gives the result.
      if (model.type %in% "exp") {
        mask.zero.rows <- d.tmp.pk.limit[,s] %in% 0
        if (any(mask.zero.rows)) {
          warning(sprintf("When estimating PK stopping limits, dropped %d zero values.", sum(mask.zero.rows)))
          d.tmp.pk.limit <- d.tmp.pk.limit[!mask.zero.rows,]
        }
        tmp.other.types <- simplify.other.types(other.types, d.tmp.pk.limit)
        other.types.txt <- get.other.types.txt(tmp.other.types)
        my.formula <-
          as.formula(sprintf("log(DOSE)~log(%s)%s", s, other.types.txt))

        ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        my.plot.formula<- as.formula(sprintf("pred+DOSE~%s%s", s,
                                             get.other.types.txt.plt(tmp.other.types)))

      } else if (model.type %in% "linear") {
        other.types.txt <- get.other.types.txt(other.types)
        my.formula <-
          as.formula(sprintf("DOSE~%s%s-1", s, other.types.txt))
      } else {
        stop("There is a bug around assignment of the model type for PK prediction.")
      }
      print(my.formula)
      pred.dose <- rep(NA, nrow(output.data[["blinded"]][["PK.Stopping"]]))

      if ((n.obs > n.sub) & (model.type %in% "exp")) {
        ## crossover study
        try({
          print(sprintf("Using linear mixed effects model to predict %s dose limit.", s))
          d.mod <- lme(my.formula,
                       random=~1|RANDID,
                       data=d.tmp.pk.limit,
                       na.action=na.omit)
          pred.dose <-
            predict(d.mod, output.data[["blinded"]][["PK.Stopping"]], level=0)
        })
      } else {
        ## parallel study
        try({
          print(sprintf("Using linear model to predict %s dose limit.", s))
          d.mod <- gls(my.formula,
                       data=d.tmp.pk.limit,
                       na.action=na.omit)
          pred.dose <-
            predict(d.mod, newdata=output.data[["blinded"]][["PK.Stopping"]])
        })
      }

      ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      ##        PLOTS
      ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      ## Create plots only if number of doses >1
      if (length(unique(d.ind.nca$DOSE)) > 1) {
        ## Create extensions to data so that predicted lines stretch
        ## to the end of panel
        my.data.extn <-
          renameCol(data.frame(x=c(0.2*min(d.tmp.pk.limit[,s]),
                                 max(5*max(d.tmp.pk.limit[,s]),
                                     d.pk.stopping[,s]))),
                    "x", s)
        if(length(unique(d.tmp.pk.limit[,other.types]))!=0){
          my.data.extn <- merge(my.data.extn,
                                unique(d.tmp.pk.limit[,other.types]))
        }
        my.data.extn$DOSE <- NA

        ## Predict data using the established lme model and also
        ## calulate the limits
        my.pred <- cbind(rbind(d.tmp.pk.limit[!is.na(d.tmp.pk.limit[,s]),
                                              c("DOSE",s,other.types)],
                               my.data.extn),
                         pred=exp(predict(d.mod,
                           rbind(d.tmp.pk.limit[!is.na(d.tmp.pk.limit[,s]),
                                                c("DOSE",s,other.types)],
                                 my.data.extn), level=0))
                         )
        ## Create a column with the original stopping criteria
        my.pred$pk.stop<-rep(d.pk.stopping[,s])

        my.limits <- as.data.frame(cbind(
          limits=exp(predict(d.mod, output.data[["blinded"]][["PK.Stopping"]],
            level=0)),
          output.data[["blinded"]][["PK.Stopping"]][, other.types]))
        names(my.limits) <- c("limits", other.types)

        if(length(other.types>0)){
          my.pred<-merge(my.pred,my.limits) }else{
            my.pred$limits<-rep(unique(my.limits$limits),nrow(my.pred))}

        ## Append plot
        output.data[["unblinded"]][["PK Stopping Limit"]]<-
          append(output.data[["unblinded"]][["PK Stopping Limit"]],
                 list(xyplot(my.plot.formula,
                             data=my.pred[order(my.pred[,s]),],
                             type=c("l","p"), distribute.type=TRUE,
                             LIMITS=my.pred$limits,
                             PKSTOP=my.pred$pk.stop,
                             scale=list(alternating=FALSE,
                               x=list(log=TRUE),
                               y=list(log=TRUE)),
                             xscale.components=xscale.components.log10ticks,
                             yscale.components=yscale.components.log10ticks,
                             xlim=c(0.7*min(my.pred[,s], na.rm=TRUE),
                               2*max(max(my.pred[!is.na(my.pred$DOSE),s],
                                         na.rm=TRUE),
                                     unique(my.pred$pk.stop))) ,
                              ylim=c(0.7*min(my.pred$DOSE, na.rm=TRUE),
                                2*max(max(my.pred$DOSE, na.rm=TRUE),
                                      max(my.pred$pred))),
                             xlab=sprintf("%s",s),
                             ylab="Dose",
                             main=paste(my.formula[2],
                               my.formula[1], my.formula[3]),
                             panel=function(x, y, ..., LIMITS, PKSTOP) {
                               panel.xyplot(x, y, ...,
                                            pch=19, col="black", cex=1.1);
                               panel.abline(h=unique(log(LIMITS,10))[panel.number()], lty=2) ;
                               panel.abline(v=unique(log(PKSTOP,10)),
                                            col="red", lwd=2, lty=2)})
                      ))
      }

      ##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      if (model.type %in% "exp")
        pred.dose <- exp(pred.dose)
      output.data[["blinded"]][["PK.Stopping"]][,paste("DOSE", s, sep=".")] <-
        pred.dose
    }
  }
}
