% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vpcstats.R
\name{binless}
\alias{binless}
\alias{binless.tidyvpcobj}
\title{binless}
\usage{
binless(o, ...)

\method{binless}{tidyvpcobj}(
  o,
  qpred = c(0.05, 0.5, 0.95),
  optimize = TRUE,
  optimization.interval = c(0, 7),
  conf.level = 0.95,
  loess.ypc = FALSE,
  lambda = NULL,
  span = NULL,
  ...
)
}
\arguments{
\item{o}{tidyvpc object}

\item{...}{other arguments}

\item{qpred}{numeric vector of length 3 specifying quantiles (lower, median, upper) i.e. \code{c(0.1, 0.5, 0.9)}}

\item{optimize}{logical indicating whether lambda and span should be optimized using AIC}

\item{optimization.interval}{numeric vector of length 2 specifying interval for lambda optimization}

\item{conf.level}{numeric confidence level for binless fit}

\item{loess.ypc}{logical indicating loess precition corrected. Must first use \code{predcorrect()} if \code{loess.ypc = TRUE}}

\item{lambda}{numeric vector of length 3 specifying lambda values for each quantile}

\item{span}{numeric number between 0,1 specying smoothing paramter for loess prediction corrected}
}
\value{
Updates \code{tidyvpcobj} with additive quantile regression fits for observed and simulated data for quantiles specified in \code{qpred} argument.
  If \code{optimize = TRUE} argument is specified, the resulting \code{tidyvpcobj} will contain optimized lambda values according to AIC.  For prediction
  corrected VPC (pcVPC), specifying \code{loess.ypc = TRUE} will return optimized span value for LOESS smoothing.
}
\description{
Perform binless Visual Predictive Check (VPC)
}
\details{
Use this function in subsitute of traditional binning methods to derive VPC using additive quantile regression and loess for pcVPC.
}
\examples{
\donttest{

 vpc <- observed(obs_data, y = DV, x = TIME) \%>\%
      simulated(sim_data, y = DV) \%>\%
      binless() \%>\%
      vpcstats()
      
 # Binless example with LOESS prediction correction
 
 obs_data$PRED <- sim_data[REP == 1, PRED]
 
 vpc <- observed(obs_data, y = DV, x = TIME) \%>\%
      simulated(sim_data, y = DV) \%>\%
      predcorrect(pred = PRED) \%>\%
      binless(optimize = TRUE, loess.ypc = TRUE) \%>\%
      vpcstats()
      
# Binless example with user specified lambda values stratified on 
# "GENDER" with 2 levels ("M", "F"), 10\%, 50\%, 90\% quantiles.
 
 lambda_strat <- data.table(
 GENDER_M = c(3,5,2),
 GENDER_F = c(1,3,4),
 )
 
 vpc <- observed(obs_data, y = DV, x = TIME) \%>\%
      simulated(sim_data, y = DV) \%>\%
      stratify(~ GENDER) \%>\%
      binless(qpred = c(0.1, 0.5, 0.9), optimize = FALSE, lambda = lambda_strat) \%>\%
      vpcstats()
}
}
\seealso{
\code{\link{observed}} \code{\link{simulated}} \code{\link{censoring}} \code{\link{predcorrect}} \code{\link{stratify}} \code{\link{binning}} \code{\link{vpcstats}}
}
