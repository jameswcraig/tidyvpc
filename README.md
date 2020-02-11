tidyvpc
========
<!-- badges: start -->
  [![Travis build status](https://travis-ci.org/jameswcraig/tidyvpc.svg?branch=master)](https://travis-ci.org/jameswcraig/tidyvpc)
  <!-- badges: end -->
  
### Installation and Running information
```
# Install the development version from GitHub:
# install.packages("devtools")
devtools::install_github("jameswcraig/tidyvpc")

```

### Usage

``` r
library(vpc)
library(magrittr)
library(ggplot2)
library(tidyvpc)

exampleobs <- as.data.table(vpc::simple_data$obs)[MDV == 0]
examplesim <- as.data.table(vpc::simple_data$sim)[MDV == 0]
exampleobs$PRED <- examplesim[REP == 1, PRED]
exampleobs$LLOQ <- exampleobs[, ifelse(ISM == 0, 100, 25)]

# Binning Method on x-variable (TIME)
vpc <- observed(exampleobs, x=TIME, y=DV) %>%
    simulated(examplesim, y=DV) %>%
    censoring(blq=(DV < LLOQ), lloq=LLOQ) %>%
    stratify(~ ISM) %>%
    binning(TIME) %>%
    predcorrect(pred=PRED) %>%
    vpcstats()

```
![Example](./inst/img/snapshot1.png)

Plot Code:

``` r
ggplot(vpc$stats, aes(x=xbin)) +
    facet_grid(~ ISM) +
    geom_ribbon(aes(ymin=lo, ymax=hi, fill=qname, col=qname, group=qname), alpha=0.1, col=NA) +
    geom_line(aes(y=md, col=qname, group=qname)) +
    geom_line(aes(y=y, linetype=qname), size=1) +
    geom_hline(data=unique(exampleobs[, .(ISM, LLOQ)]),
        aes(yintercept=LLOQ), linetype="dotted", size=1) +
    geom_text(data=unique(exampleobs[, .(ISM, LLOQ)]),
        aes(x=10, y=LLOQ, label=paste("LLOQ", LLOQ, sep="="),), vjust=-1) +
    scale_colour_manual(
        name="Simulated Percentiles\nMedian (lines) 95% CI (areas)",
        breaks=c("q0.05", "q0.5", "q0.95"),
        values=c("red", "blue", "red"),
        labels=c("5%", "50%", "95%")) +
    scale_fill_manual(
        name="Simulated Percentiles\nMedian (lines) 95% CI (areas)",
        breaks=c("q0.05", "q0.5", "q0.95"),
        values=c("red", "blue", "red"),
        labels=c("5%", "50%", "95%")) +
    scale_linetype_manual(
        name="Observed Percentiles\n(black lines)",
        breaks=c("q0.05", "q0.5", "q0.95"),
        values=c("dotted", "solid", "dashed"),
        labels=c("5%", "50%", "95%")) +
    guides(
        fill=guide_legend(order=2),
        colour=guide_legend(order=2),
        linetype=guide_legend(order=1)) +
    theme(
        legend.position="top",
        legend.key.width=grid::unit(2, "cm")) +
    labs(x="Time (h)", y="Concentration (ng/mL)")
```
Or use the built in `plot()` function from the tidyvpc package.

```{r}
# Binless method using 10%, 50%, 90% qunatiles and LOESS Prediction Corrected
vpc <- observed(exampleobs, x=TIME, y=DV) %>%
    simulated(examplesim, y=DV) %>%
    stratify(~ ISM) %>%
    predcorrect(pred=PRED) %>%
    binless(qpred = c(0.1, 0.5, 0.9), optimize = TRUE, loess.ypc = TRUE) %>%
    vpcstats()

plot(vpc)
```

![Example](./inst/img/snapshot2.png)