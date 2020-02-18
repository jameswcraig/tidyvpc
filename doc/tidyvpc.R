## ---- warning = FALSE, echo = FALSE, message = FALSE---------------------
knitr::opts_chunk$set(collapse = T, comment = "#>")
options(datatable.print.nrows = 8)
library(tidyvpc)
library(ggplot2)
set.seed(1014)

## ----message=FALSE-------------------------------------------------------
obs_data <- as.data.table(tidyvpc::obs_data)
sim_data <- as.data.table(tidyvpc::sim_data)
head(obs_data)


## ------------------------------------------------------------------------
obs_data <- obs_data[MDV == 0]
sim_data <- sim_data[MDV == 0]

## ------------------------------------------------------------------------
obs_data$PRED <- sim_data[REP == 1, PRED]

## ----message = FALSE-----------------------------------------------------
library(tidyvpc)
vpc <- observed(obs_data, x = TIME, y = DV)

## ------------------------------------------------------------------------
vpc <- observed(obs_data, x = TIME, y = DV) %>%
  simulated(sim_data, y = DV)

## ----message=FALSE, fig.width = 9, fig.height = 6, out.width=640---------
vpc <- observed(obs_data, x=TIME, y=DV) %>%
    simulated(sim_data, y=DV) %>%
    binning(bin = NTIME)

## ------------------------------------------------------------------------
vpc <- observed(obs_data, x = TIME, y = DV) %>%
  simulated(sim_data, y = DV) %>%
  binning(bin = "ntile", nbins = 9)

## ------------------------------------------------------------------------
vpc <- observed(obs_data, x = TIME, y = DV) %>%
  simulated(sim_data, y = DV) %>%
  binning(bin = "breaks", breaks = c(1,5,7,9,10))

## ----warning=FALSE-------------------------------------------------------
vpc <- observed(obs_data, x=TIME, y=DV) %>%
    simulated(sim_data, y=DV) %>%
    binless(optimize = TRUE)

## ----warning = FALSE-----------------------------------------------------
vpc <- observed(obs_data, x=TIME, y=DV) %>%
    simulated(sim_data, y=DV) %>%
    binless(qpred = c(.1, .5, .9), optimize = FALSE, lambda = c(1,3,2))

## ----fig.width = 9, fig.height = 6, out.width=640------------------------
vpc <- observed(obs_data, x=TIME, y=DV) %>%
    simulated(sim_data, y=DV) %>%
    binning(bin = NTIME) %>%
    vpcstats()

plot(vpc)

## ----fig.width = 9, fig.height = 6, out.width=640------------------------
vpc <- observed(obs_data, x=TIME, y=DV) %>%
    simulated(sim_data, y=DV) %>%
    stratify(~ GENDER) %>%
    binning(bin = "pam", nbins = 7) %>%
    vpcstats()

plot(vpc)

## ----warning = FALSE, fig.width = 9, fig.height = 6, out.width=640-------

vpc <- observed(obs_data, x=TIME, y=DV) %>%
    simulated(sim_data, y=DV) %>%
    stratify(~ GENDER + GROUP) %>%
    binless() %>%
    vpcstats()

plot(vpc)

## ----fig.width = 9, fig.height = 6, out.width=640, warning = FALSE-------
vpc <- observed(obs_data, x=TIME, y=DV) %>%
    simulated(sim_data, y=DV) %>%
    censoring(blq=(DV < 25), lloq=25) %>%
    binning(bin = "ntile", nbins = 9) %>%
    vpcstats()

plot(vpc)

## ------------------------------------------------------------------------
obs_data$LLOQ <- 50

## ----fig.width = 9, fig.height = 6, out.width=640, warning = FALSE-------
vpc <- observed(obs_data, x=TIME, y=DV) %>%
    simulated(sim_data, y=DV) %>%
    censoring(blq=(DV < LLOQ), lloq=LLOQ) %>%
    binless(optimize = FALSE, lambda = c(1.5, 2.5, 1)) %>%
    vpcstats()

plot(vpc)

## ----fig.width = 9, fig.height = 6, out.width=640, warning = FALSE-------
obs_data$LLOQ <- obs_data[, ifelse(GENDER == "F", 100, 25)]

vpc <- observed(obs_data, x=TIME, y=DV) %>%
    simulated(sim_data, y=DV) %>%
    censoring(blq=(DV < LLOQ), lloq=LLOQ) %>%
    stratify(~ GENDER) %>%
    binless() %>%
    vpcstats()

plot(vpc)

## ----fig.width = 9, fig.height = 6, out.width=640------------------------
vpc <- observed(obs_data, x=TIME, y=DV) %>%
    simulated(sim_data, y=DV) %>%
    binning(bin = "jenks", nbins = 7) %>%
    predcorrect(pred=PRED) %>%
    vpcstats()

plot(vpc)

## ----warning=FALSE, fig.width = 9, fig.height = 6, out.width=640---------
vpc <- observed(obs_data, x=TIME, y=DV) %>%
    simulated(sim_data, y=DV) %>%
    stratify(~ GENDER) %>%
    predcorrect(pred=PRED) %>%
    binless(qpred = c(0.1, 0.5, 0.9), optimize = TRUE, loess.ypc = TRUE) %>%
    vpcstats()

plot(vpc)

## ----fig.width = 9, fig.height = 6, out.width=640, warning = FALSE-------
vpc <- observed(obs_data, x=TIME, y=DV) %>%
    simulated(sim_data, y=DV) %>%
    predcorrect(pred=PRED) %>%
    binless(optimize = FALSE, lambda = c(1,3, 1.5), loess.ypc = TRUE, span = .4) %>%
    vpcstats()

plot(vpc)

## ----fig.width = 9, fig.height = 6, out.width=640------------------------
vpc <- observed(obs_data, x=TIME, y=DV) %>%
    simulated(sim_data, y=DV) %>%
    stratify(~ GENDER + GROUP) %>%
    binning(stratum = list(GENDER = "M"), bin = "jenks", nbins = 7, by.strata = T) %>%
    binning(stratum = list(GENDER = "F"), bin = "centers", centers = c(1,3,5,7), by.strata = T) %>%
    binning(stratum = list(GROUP = "Concomitant"), bin = "kmeans", by.strata = T) %>%
    binning(stratum = list(GROUP = "No Concomitant"), bin = "pam", nbins = 10, by.strata = T) %>%
    predcorrect(pred=PRED) %>%
    vpcstats()

plot(vpc)

## ----warning = FALSE, fig.width = 9, fig.height = 6, out.width=640-------
new_lambda = data.frame(GENDERM = c(1,3,.8), GENDERF = c(2,4,1))

vpc <- observed(obs_data, x=TIME, y=DV) %>%
    simulated(sim_data, y=DV) %>%
    stratify(~ GENDER) %>%
    predcorrect(pred=PRED) %>%
    binless(optimize = FALSE, lambda = new_lambda, loess.ypc = TRUE, span = c(.43, .85)) %>%
    vpcstats()

plot(vpc)

## ----fig.width = 9, fig.height = 6, out.width=640, warning = FALSE-------
vpc <- observed(obs_data, x=TIME, y=DV) %>%
    simulated(sim_data, y=DV) %>%
    binning(bin = "jenks", nbins = 7) %>%
    vpcstats()

plot(vpc, show.stats = FALSE)


## ------------------------------------------------------------------------
vpc <- observed(obs_data, x=TIME, y=DV) %>%
    simulated(sim_data, y=DV) %>%
    binning(bin = "hclust", nbins = 8) %>%
    vpcstats()

bininfo(vpc)


## ----fig.width = 9, fig.height = 6, out.width=640------------------------
vpc <- observed(obs_data, x=TIME, y=DV) %>%
    simulated(sim_data, y=DV) %>%
    censoring(blq=(DV < 25), lloq=25) %>%
    binning(bin = "ntile", nbins = 9) %>%
    vpcstats()

ggplot(vpc$pctblq, aes(x=xbin, y=md))+
  geom_line()+
  geom_point() + 
  ylab("%blq")

## ----fig.width = 9, fig.height = 6, out.width=640, warning=FALSE---------
obs_data$LLOQ <- 50

vpc <- observed(obs_data, x = TIME, y = DV) %>%
  simulated(sim_data, y = DV) %>%
  censoring(blq = DV < LLOQ, lloq = LLOQ) %>%
  stratify(~GENDER) %>%
  binning(bin = NTIME) %>%
  vpcstats()

ggplot(vpc$stats, aes(x = xbin)) + 
  facet_grid(~GENDER, scales = "free", as.table = FALSE) + 
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = qname, col = qname, group = qname),alpha = 0.1, col = NA) + 
  geom_line(aes(y = md, col = qname, group = qname)) +
  geom_line(aes(y = y, linetype = qname), size = 1) + 
  geom_hline(data = unique(vpc$data[, .(LLOQ), by = eval("GENDER")]), 
            aes(yintercept = !!as.symbol("LLOQ")), linetype = "dotted", size = 1) +
  geom_text(data = unique(vpc$data[, .(LLOQ), by = eval("GENDER")]), 
            aes(x = 10, y = LLOQ, label = paste("LLOQ", LLOQ, sep = "="), ), vjust = 1, hjust = 1) +
  scale_colour_manual(name = "Simulated Percentiles\nMedian (lines) 95% CI (areas)",
                    breaks = c("q0.05", "q0.5", "q0.95"), 
                    values = c("red", "blue", "red"), 
                    labels = c("5%", "50%", "95%")) + 
  scale_fill_manual(name = "Simulated Percentiles\nMedian (lines) 95% CI (areas)", 
                    breaks = c("q0.05", "q0.5", "q0.95"), 
                    values = c("red", "blue", "red"), 
                    labels = c("5%", "50%", "95%")) + 
  scale_linetype_manual(name = "Observed Percentiles\nMedian (lines) 95% CI (areas)", 
                    breaks = c("q0.05", "q0.5", "q0.95"), 
                    values = c("dotted", "solid", "dashed"), 
                    labels = c("5%", "50%", "95%")) + 
  guides(fill = guide_legend(order = 2), colour = guide_legend(order = 2), linetype = guide_legend(order = 1)) + 
  theme(legend.position = "top", legend.key.width = grid::unit(2, "cm")) + 
  labs(x = "NTIME", y = "Concentration") + 
  geom_point(data = vpc$obs, aes(x = x, y = y), size = 1, alpha = 0.1, show.legend = F) + 
  geom_vline(data = bininfo(vpc)[, .(x = sort(unique(c(xleft, xright)))), by = names(vpc$strat)],aes(xintercept = x), size = rel(0.5), col = "gray80") + 
  theme(panel.grid = element_blank()) + 
  geom_rug(data = bininfo(vpc)[, .(x = sort(unique(c(xleft, xright)))), by = names(vpc$strat)],aes(x = x), sides = "t", size = 1)


## ----fig.width = 9, fig.height = 6, out.width=640, warning = FALSE-------
ggplot(vpc$pctblq) + 
  facet_grid(~GENDER) +
  geom_ribbon(aes(x = xbin, ymin= lo, ymax = hi), fill = "red", alpha = .2) + 
  geom_line(aes(x = xbin, y = y)) + 
  geom_line(aes(x = xbin, y = md), color = "red") + 
  labs(x= "NTIME", y= "% BLQ")

