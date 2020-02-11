## ---- warning = FALSE, echo = FALSE, message = FALSE---------------------
knitr::opts_chunk$set(collapse = T, comment = "#>")
options(datatable.print.nrows = 8)
library(tidyvpc)
library(ggplot2)
set.seed(1014)

## ----message=FALSE-------------------------------------------------------
exampleobs <- tidyvpc::exampleobs
examplesim <- tidyvpc::examplesim
head(exampleobs)


## ------------------------------------------------------------------------
exampleobs <- exampleobs[MDV == 0]
examplesim <- examplesim[MDV == 0]

## ------------------------------------------------------------------------
exampleobs$PRED <- examplesim[REP == 1, PRED]

## ------------------------------------------------------------------------
exampleobs$LLOQ <- exampleobs[, ifelse(ISM == 0, 100, 25)]

## ----message = FALSE-----------------------------------------------------
library(tidyvpc)
vpc <- observed(exampleobs, x = TIME, y = DV)

## ------------------------------------------------------------------------
vpc <- observed(exampleobs, x = TIME, y = DV) %>%
  simulated(examplesim, y = DV)

## ----message=FALSE, fig.width = 9, fig.height = 6, out.width=640---------
vpc <- observed(exampleobs, x=TIME, y=DV) %>%
    simulated(examplesim, y=DV) %>%
    binning(bin = TIME)

## ------------------------------------------------------------------------
vpc <- observed(exampleobs, x = TIME, y = DV) %>%
  simulated(examplesim, y = DV) %>%
  binning(bin = "ntile", nbins = 9)

## ------------------------------------------------------------------------
vpc <- observed(exampleobs, x = TIME, y = DV) %>%
  simulated(examplesim, y = DV) %>%
  binning(bin = "breaks", breaks = c(1,5,7,9,10))

## ------------------------------------------------------------------------
vpc <- observed(exampleobs, x=TIME, y=DV) %>%
    simulated(examplesim, y=DV) %>%
    binless(optimize = TRUE)

## ------------------------------------------------------------------------
vpc <- observed(exampleobs, x=TIME, y=DV) %>%
    simulated(examplesim, y=DV) %>%
    binless(qpred = c(.1, .5, .9), optimize = FALSE, lambda = c(1,3,2))

## ----fig.width = 9, fig.height = 6, out.width=640------------------------
vpc <- observed(exampleobs, x=TIME, y=DV) %>%
    simulated(examplesim, y=DV) %>%
    binning(bin = TIME) %>%
    vpcstats()

plot(vpc)

## ----fig.width = 9, fig.height = 6, out.width=640------------------------
vpc <- observed(exampleobs, x=TIME, y=DV) %>%
    simulated(examplesim, y=DV) %>%
    stratify(~ ISM) %>%
    binning(bin = "pam", nbins = 8) %>%
    vpcstats()

plot(vpc)

## ----warning = FALSE, fig.width = 9, fig.height = 6, out.width=640-------
#Include dummy variable in data for example.
exampleobs$SEX <- rep(c("F", "M"), len=nrow(exampleobs))

vpc <- observed(exampleobs, x=TIME, y=DV) %>%
    simulated(examplesim, y=DV) %>%
    stratify(~ ISM + SEX) %>%
    binless() %>%
    vpcstats()

plot(vpc)

## ----fig.width = 9, fig.height = 6, out.width=640------------------------
vpc <- observed(exampleobs, x=TIME, y=DV) %>%
    simulated(examplesim, y=DV) %>%
    censoring(blq=(DV < LLOQ), lloq=25) %>%
    binning(TIME) %>%
    vpcstats()

plot(vpc)

## ----fig.width = 9, fig.height = 6, out.width=640------------------------
vpc <- observed(exampleobs, x=TIME, y=DV) %>%
    simulated(examplesim, y=DV) %>%
    censoring(blq=(DV < LLOQ), lloq=LLOQ) %>%
    stratify(~ ISM) %>%
    binless() %>%
    vpcstats()

plot(vpc)

## ----fig.width = 9, fig.height = 6, out.width=640------------------------
vpc <- observed(exampleobs, x=TIME, y=DV) %>%
    simulated(examplesim, y=DV) %>%
    binning(bin = "jenks", nbins = 7) %>%
    predcorrect(pred=PRED) %>%
    vpcstats()

plot(vpc)

## ----warning=FALSE, fig.width = 9, fig.height = 6, out.width=640---------
vpc <- observed(exampleobs, x=TIME, y=DV) %>%
    simulated(examplesim, y=DV) %>%
    stratify(~ ISM) %>%
    predcorrect(pred=PRED) %>%
    binless(qpred = c(0.1, 0.5, 0.9), optimize = TRUE, loess.ypc = TRUE) %>%
    vpcstats()

plot(vpc)

## ----fig.width = 9, fig.height = 6, out.width=640------------------------
vpc <- observed(exampleobs, x=TIME, y=DV) %>%
    simulated(examplesim, y=DV) %>%
    predcorrect(pred=PRED) %>%
    binless(optimize = FALSE, lambda = c(1,3, 1.5), loess.ypc = TRUE, span = .4) %>%
    vpcstats()

plot(vpc)

## ----fig.width = 9, fig.height = 6, out.width=640------------------------
vpc <- observed(exampleobs, x=TIME, y=DV) %>%
    simulated(examplesim, y=DV) %>%
    stratify(~ ISM + SEX) %>%
    binning(stratum = list(ISM = 0), bin = "jenks", nbins = 7, by.strata = T) %>%
    binning(stratum = list(ISM = 1), bin = "centers", centers = c(1,3,5,7), by.strat = T) %>%
    binning(stratum = list(SEX = "M"), bin = TIME, by.strata = T) %>%
    binning(stratum = list(SEX = "F"), bin = "pam", nbins = 10, by.strata = T) %>%
    predcorrect(pred=PRED) %>%
    vpcstats()

plot(vpc)

## ----warning = FALSE, fig.width = 9, fig.height = 6, out.width=640-------
new_lambda = data.frame(ISM0 = c(1,3,5), ISM1 = c(2,4,1))

vpc <- observed(exampleobs, x=TIME, y=DV) %>%
    simulated(examplesim, y=DV) %>%
    stratify(~ ISM) %>%
    predcorrect(pred=PRED) %>%
    binless(optimize = FALSE, lambda = new_lambda, loess.ypc = TRUE, span = c(.43, .85)) %>%
    vpcstats()

plot(vpc)

