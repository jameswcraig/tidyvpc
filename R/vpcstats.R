#' Perform a Visual Predictive Check (VPC) computation
#'
#' These functions work together to calculate the statistics that are plotted
#' in a VPC. They would typically be chained together using the "pipe" operator
#' (see Examples).
#'
#' @param o An object.
#' @param ... Additional arguments.
#'
#' @import data.table
#' @import magrittr
#' @import quantreg
#' @importFrom stats median model.frame quantile setNames update AIC fitted loess na.exclude optimize resid time
#' @importFrom utils install.packages installed.packages
#' @name generics
NULL

#' Specify observed dataset and variables for VPC
#' 
#' The observed function is the first function in the vpc piping chain and is used for specifying observed data and variables for VPC
#' 
#' @title observed
#' @param o data.frame or data.table of observation data
#' @param x numeric x-variable, typically named TIME
#' @param yobs numeric y-variable, typically named DV
#' @param pred population prediction variable, typically named PRED
#' @param blq logical variable indicating below limit of quantification 
#' @param lloq number or numeric variable in data indicating the lower limit of quantification
#' @param alq logical variable indicating above limit of quantification 
#' @param uloq number or numeric variable in data indicating the upper limit of quantification
#' @param ... other arguments
#' @examples
#' \dontrun{
#' vpc <- observed(exampleobs, x=TIME, y=DV) %>%
#'     simulated(examplesim, y=DV) %>%
#'     binning(bin = "ntile", nbins = 8) %>%
#'     vpcstats()
#'     }
#' @seealso \code{\link{simulated}} \code{\link{censoring}} \code{\link{stratify}} \code{\link{predcorrect}} \code{\link{binning}} \code{\link{binless}} \code{\link{vpcstats}}
#' @export
observed <- function(o, ...) UseMethod("observed")

#' @rdname observed
#' @export
observed.data.frame <- function(o, x, yobs, pred=NULL, blq, lloq=-Inf, alq, uloq=Inf, ...) {
  data <- o
  x    <- rlang::eval_tidy(rlang::enquo(x),    data)
  yobs <- rlang::eval_tidy(rlang::enquo(yobs), data)
  pred <- rlang::eval_tidy(rlang::enquo(pred), data)
  lloq <- rlang::eval_tidy(rlang::enquo(lloq), data)
  uloq <- rlang::eval_tidy(rlang::enquo(uloq), data)
  lloq <- as.numeric(lloq)
  uloq <- as.numeric(uloq)
  
  if (missing(blq)) {
    blq <- (yobs < lloq)
  } else {
    blq <- rlang::eval_tidy(rlang::enquo(blq),  data)
  }
  blq  <- as.logical(blq)
  
  if (missing(alq)) {
    alq <- (yobs > uloq)
  } else {
    alq <- rlang::eval_tidy(rlang::enquo(alq),  data)
  }
  alq <- as.logical(alq)

  obs <- data.table(x, y=yobs, blq, lloq, alq, uloq)
  
  o <- structure(list(data=data), class="tidyvpcobj")
  update(o, obs=obs, pred=pred)
}

#' Specify simulated dataset and variables for VPC
#' 
#' The simulated function is the second function in the vpc piping chain and is used for specifying simulated data and variables for VPC
#' 
#' @title observed
#' @param o tidyvpcobj
#' @param data data.frame or data.table of simulated data
#' @param ysim numeric y-variable, typically named DV
#' @param ... other arguments
#' @examples
#' \dontrun{
#' vpc <- observed(exampleobs, x=TIME, y=DV) %>%
#'     simulated(examplesim, y=DV) %>%
#'     binning(bin = "pam", nbins = 8) %>%
#'     vpcstats()
#'     }
#' @seealso \code{\link{observed}} \code{\link{censoring}} \code{\link{stratify}} \code{\link{predcorrect}} \code{\link{binning}} \code{\link{binless}} \code{\link{vpcstats}}

#' @export
simulated <- function(o, ...) UseMethod("simulated")

#' @rdname simulated
#' @export
simulated.tidyvpcobj <- function(o, data, ysim, ...) {
  ysim <- rlang::eval_tidy(rlang::enquo(ysim), data)
  obs  <- o$obs
  x    <- obs$x
  nrep <- length(ysim)/nrow(obs)
  repl <- rep(1:nrep, each=nrow(obs))
  
  sim <- data.table(x, y=ysim, repl)
  update(o, sim=sim)
}

#' Censoring observed data for Visual Predictive Check (VPC)
#' 
#' Specify censoring variables or censoring value for VPC using this function
#' 
#' @title censoring
#' @param o tidyvpc object
#' @param blq blq variable if present in observed data 
#' @param lloq lloq variable if present in observed data. Use numeric to specify lloq value
#' @param alq logical variable indicating above limit of quantification 
#' @param uloq number or numeric variable in data indicating the upper limit of quantification
#' @param data observed data supplied in \code{observed()} function
#' @param ... Other arguments to include
#' @examples
#' \dontrun{
#' vpc <- observed(exampleobs, x=TIME, y=DV) %>%
#'     simulated(examplesim, y=DV) %>%
#'     censoring(blq=(DV < LLOQ), lloq=50) %>%
#'     binning(TIME) %>%
#'     vpcstats()
#' 
#' #Using LLOQ variable in data:
#' 
#' exampleobs$LLOQ <- exampleobs[, ifelse(ISM == 0, 100, 25)]
#' 
#' vpc <- observed(exampleobs, x=TIME, y=DV) %>%
#'     simulated(examplesim, y=DV) %>%
#'     censoring(blq=(DV < LLOQ), lloq=LLOQ) %>%
#'     stratify(~ ISM) %>%
#'     binning(TIME) %>%
#'     vpcstats()
#'
#' } 
#' @seealso \code{\link{observed}} \code{\link{simulated}} \code{\link{stratify}} \code{\link{predcorrect}} \code{\link{binning}} \code{\link{binless}} \code{\link{vpcstats}}

#' @export 
censoring <- function(o, ...) UseMethod("censoring")

#' @rdname censoring
#' @export
censoring.tidyvpcobj <- function(o, blq, lloq, alq, uloq, data=o$data, ...) {
  if (missing(blq)) {
    blq <- o$obs$blq
  } else {
    blq <- rlang::eval_tidy(rlang::enquo(blq), data)
  }
  if (missing(lloq)) {
    lloq <- o$obs$lloq
  } else {
    lloq <- rlang::eval_tidy(rlang::enquo(lloq), data)
  }
  if (missing(alq)) {
    alq <- o$obs$alq
  } else {
    alq <- rlang::eval_tidy(rlang::enquo(alq), data)
  }
  if (missing(uloq)) {
    uloq <- o$obs$uloq
  } else {
    uloq <- rlang::eval_tidy(rlang::enquo(uloq), data)
  }

  if (is.null(blq)) {
    stop("No blq specified")
  }
  if (is.null(alq)) {
      browser()
    stop("No alq specified")
  }
  if (!is.null(blq) & is.null(lloq)) {
    stop("No lloq specified for blq")
  }
  if (!is.null(alq) & is.null(uloq)) {
    stop("No uloq specified for alq")
  }
  
  .blq <- blq
  .lloq <- lloq
  .alq <- alq
  .uloq <- uloq
  o$obs[, blq := rep(.blq, len=.N)]
  o$obs[, lloq := rep(.lloq, len=.N)]
  o$obs[, alq := rep(.alq, len=.N)]
  o$obs[, uloq := rep(.uloq, len=.N)]
  
  update(o, censoring=TRUE)
}

#' Stratification for Visual Predictive Check (VPC)
#' 
#' specify stratification variables for VPC using this function
#' 
#' @title stratify
#' @param o tidyvpc object
#' @param formula formula for stratfication
#' @param data Observed data supplied in \code{observed()} function
#' @param ... Other arguments to include
#' @examples 
#' \dontrun{
#' vpc <- observed(exampleobs, x=TIME, y=DV) %>%
#'     simulated(examplesim, y=DV) %>%
#'     stratify(~ ISM) %>%
#'     binning(TIME) %>%
#'     vpcstats()
#'
#' # Example with 2-way stratification
#' # Add in SEX dummy variable to observed data
#' 
#' exampleobs$SEX <- rep(c("F", "M"), len=nrow(exampleobs))
#'
#' vpc <- vpc %>%
#'     stratify(~ ISM + SEX) %>%
#'     binning(TIME) %>%
#'     vpcstats()
#'}
#' @seealso \code{\link{observed}} \code{\link{simulated}} \code{\link{censoring}} \code{\link{predcorrect}} \code{\link{binning}} \code{\link{binless}} \code{\link{vpcstats}}

#' @export 
stratify <- function(o, ...) UseMethod("stratify")

#' @method stratify tidyvpcobj
#' @rdname stratify
#' @export
stratify.tidyvpcobj <- function(o, formula, data=o$data, ...) {
  if (!inherits(formula, "formula")) {
    stop("Expecting a formula")
  }
  flist <- as.list(formula)
  if (length(flist) == 3) {
    lhs <- as.call(c(flist[[1]], flist[[2]]))
    rhs <- as.call(c(flist[[1]], flist[[3]]))
    if (flist[[2]] == as.symbol(".")) {
      lhsmf <- NULL
    } else {
      lhsmf <- as.data.table(model.frame(lhs, data))
    }
    if (flist[[3]] == as.symbol(".")) {
      rhsmf <- NULL
    } else {
      rhsmf <- as.data.table(model.frame(rhs, data))
    }
    if (is.null(lhsmf) && is.null(rhsmf)) {
      stop("Invalid stratification formula: no variables specified")
    }
    strat <- cbind(lhsmf, rhsmf)
  } else {
    strat <- as.data.table(model.frame(formula, data))
  }
  
  reserved.names <- c("x", "y", "ypc", "pred", "blq", "lloq", "repl", "bin", "xbin", "qname", "lo", "md", "hi",
                      "nobs", "xmedian", "xmean", "xmin", "xmax", "xmid", "xleft", "xright", "xcenter")
  if (any(names(strat) %in% reserved.names)) {
    stop(paste0("The names of used for stratification must not include: ",
                paste0(reserved.names, collapse=", ")))
  }
  
  o$obs[, names(strat) := strat]
  
  strat.split <- split(o$obs, strat)
  
  update(o, strat=strat, strat.split = strat.split, strat.formula=formula)
}

#' Binning methods for Visual Predictive Check (VPC)
#' 
#' This function executes binning methods available in classInt i.e. "jenks", "kmeans", "sd", "pretty", "pam", "kmeans", "hclust", "bclust", "fisher", and "dpih".
#' You may also bin directly on x-variable or alternatively specify "centers" or "breaks". For explanation of binning methods see \code{\link[classInt]{classIntervals}}
#' 
#' @title binning
#' @param o tidyvpc object
#' @param bin Character string indicating binning method or unquoted variable name if binning on x-variable. 
#' @param data Observed data supplied in \code{observed()} function
#' @param xbin Character string indicating midpoint type for binning
#' @param centers Numeric vector of centers for binning. Use \code{bin = "centers"} if supplying centers
#' @param breaks Numeric vector of breaks for binning. Use \code{bin = "breaks"} if supplying breaks
#' @param nbins Numeric number indicating the number of bins to use
#' @param altx  Unquoted variable name in observed data for elternative x-variable binning
#' @param stratum List indicating the name of stratification variable and level if using different binning methods by strata
#' @param by.strata Logical indicating whether binning should be performed by strata
#' @param ... Other arguments to include
#' @seealso \code{\link{observed}} \code{\link{simulated}} \code{\link{censoring}} \code{\link{predcorrect}} \code{\link{stratify}} \code{\link{binless}} \code{\link{vpcstats}}
#' @examples 
#' \dontrun{
#'  # Binning on x-variable
#' vpc <- observed(exampleobs, x=TIME, y=DV) %>%
#'     simulated(examplesim, y=DV) %>%
#'     binning(bin = TIME) %>%
#'     vpcstats()
#'     
#'  # Binning using ntile and xmean for midpoint
#' vpc <- observed(exampleobs, x=TIME, y=DV) %>%
#'     simulated(examplesim, y=DV) %>%
#'     binning(bin = "ntile", nbins = 8, xbin = "xmean") %>%
#'     vpcstats()
#'     
#'  # Binning using centers
#' vpc <- observed(exampleobs, x=TIME, y=DV) %>%
#'     simulated(examplesim, y=DV) %>%
#'     binning(bin = "centers", centers = c(1,3,5,7)) %>%
#'     vpcstats()
#'         
#'  # Different Binning for each level of Strata
#' vpc <- observed(exampleobs, x=TIME, y=DV) %>%
#'     simulated(examplesim, y=DV) %>%
#'     stratify(~ ISM) %>%
#'     binning(stratum = list(ISM = 0), bin = "pam", nbins = 9, by.strata = T) %>%
#'     binning(stratum = list(ISM = 1), bin = "breaks", breaks = c(2,5,7,10), by.strata = T) %>%
#'     vpcstats()
#' }
#' @export
binning <- function(o, ...) UseMethod("binning")

#' @method binning tidyvpcobj
#' @rdname binning
#' @export
binning.tidyvpcobj <- function(o, bin, data=o$data, xbin="xmedian", centers, breaks, nbins, altx, stratum=NULL, by.strata=T,  ...) {
  keep <- i <- NULL
  . <- list
  
  # If xbin is numeric, then that is the bin
  xbin <- rlang::eval_tidy(rlang::enquo(xbin), data)
  if (is.numeric(xbin)) {
    if (length(xbin) != nrow(o$obs)) {
      stop("A numeric xbin be of length equal to the number of observations")
    }
    bin <- xbin
  } else {
    if (missing(bin) && !missing(centers)) {
      bin <- "centers"
    } else if (missing(bin) && !missing(breaks)) {
      bin <- "breaks"
    } else {
      bin <- rlang::eval_tidy(rlang::enquo(bin), data)
    }
  }
  
  # If you don't want to bin on the observed x, you can specify an alternate x for binning
  if (missing(altx)) {
    x <- o$obs$x
  } else {
    x <- rlang::eval_tidy(rlang::enquo(altx), data)
  }
  
  if (!missing(nbins)) {
    nbins <- rlang::eval_tidy(rlang::enquo(nbins), data)
    if (is.numeric(nbins) && !is.null(o$strat) && (length(nbins) == nrow(o$strat))) {
      nbins <- data.table(nbins)[, .(nbins=unique(nbins)), by=o$strat]
    }
  }
  
  args <- lapply(rlang::enquos(...), rlang::eval_tidy, data=data)
  
  by.strata <- isTRUE(by.strata)
  
  # Check if specific stratum selected (can be more than 1), setup filter
  if (!is.null(stratum)) {
    if (!is.list(stratum)) {
      stop("stratum must be a list, data.frame or data.table")
    }
    if (is.null(o$strat)) {
      stop("No stratification has been specified")
    }
    if (!by.strata) {
      stop("by.strata must be TRUE when stratum is specified")
    }
    filter <- copy(o$strat)[, keep := F]
    filter[as.data.table(stratum), keep := T, on=names(stratum)]
    filter <- filter$keep
  } else {
    filter <- rep(T, nrow(o$obs))  # Keep all
  }
  
  # If xbin is numeric, then that is the bin
  if (is.numeric(xbin)) {
    if (length(xbin) != nrow(o$obs)) {
      stop("A numeric xbin be of length equal to the number of observations")
    }
    bin <- xbin
  } else if (is.character(bin) && length(bin) == 1) {
    
    known.classInt.styles <- c("fixed", "sd", "equal", "pretty", "quantile",
                               "kmeans", "hclust", "bclust", "fisher", "jenks", "dpih")
    
    if (bin == "centers") {
      if (missing(centers)) {
        stop("centers must be specified to use this binning method")
      }
      if (!by.strata && is.data.frame(centers)) {
        stop("by.strata must be TRUE when centers is a data.frame")
      }
      bin <- nearest(centers)
    } else if (bin == "breaks") {
      if (missing(breaks)) {
        stop("breaks must be specified to use this binning method")
      }
      if (!by.strata && is.data.frame(breaks)) {
        stop("by.strata must be TRUE when breaks is a data.frame")
      }
      bin <- cut_at(breaks)
    } else if (bin == "ntile") {
      if (missing(nbins)) {
        stop("nbins must be specified to use this binning method")
      }
      bin <- bin_by_ntile(nbins)
    } else if (bin == "eqcut") {
      if (missing(nbins)) {
        stop("nbins must be specified to use this binning method")
      }
      bin <- bin_by_eqcut(nbins)
    } else if (bin == "pam") {
      if (missing(nbins)) {
        stop("nbins must be specified to use this binning method")
      }
      bin <- bin_by_pam(nbins)
    } else if (bin %in% known.classInt.styles) {
      if (missing(nbins)) {
        nbins <- NULL
      }
      bin <- bin_by_classInt(bin, nbins)
    } else {
      stop(sprintf("Unknown binning method: %s", bin))
    }
  }
  
  if (is.function(bin)) {
    xdat <- data.table(i=1:nrow(o$obs), x=x)
    if (any(is.na(xdat[filter]$x))) {
      warning("x contains missing values, which could affect binning")
    }
    if (any(is.infinite(xdat[filter]$x))) {
      warning("x contains non-finite values, which could affect binning")
    }
    if (by.strata && !is.null(o$strat)) {
      sdat <- copy(o$strat)
      temp <- xdat[filter, .(i=i, j=do.call(bin, c(list(x), args, .BY))), by=sdat[filter]]
      j <- temp[order(i), j]
    } else {
      j <- xdat[filter, do.call(bin, c(list(x), args))]
    }
    if (length(j) != sum(filter)) {
      stop("The binning function did not return the right number of elements")
    }
  } else if (length(bin) == nrow(o$obs)) {
    j <- bin[filter]
  } else {
    stop("Incorrect binning specification")
  }
  o$obs[filter, bin := as.factor(j)]
  bin <- o$obs$bin
  
  if (!is.null(o$strat)) {
    stratbin <- data.table(o$strat, bin)
  } else {
    stratbin <- data.table(bin)
  }
  o <- update(o, .stratbin=stratbin, bin.by.strata=by.strata)
  
  # Assign an x value to each bin
  if (is.numeric(xbin)) {
    xbin <- data.table(xbin=xbin)[, .(xbin = unique(xbin)), by=stratbin]
  } else if (is.character(xbin) && length(xbin) == 1) {
    bi <- bininfo(o)
    xbin <- data.table(bi[, names(stratbin), with=F], xbin=bi[[xbin]])
  } else if (is.function(xbin)) {
    xbin <- data.table(x=x)[, .(xbin = xbin(x)), by=stratbin]
  } else {
    stop("Invalid xbin")
  }
  update(o, xbin=xbin)
}

#' Perform binless Visual Predictive Check (VPC)
#' 
#' Use this function in subsitute of traditional binning methods to derive VPC using additive quantile regression and loess for pcVPC.
#' 
#' @title binless
#' @param o tidyvpc object
#' @param qpred numeric vector of length 3 specifying quantiles (lower, median, upper) i.e. \code{c(0.1, 0.5, 0.9)}
#' @param optimize logical indicating whether lambda and span should be optimized using AIC
#' @param optimization.interval numeric vector of length 2 specifying interval for lambda optimization
#' @param conf.level numeric confidence level for binless fit
#' @param loess.ypc logical indicating loess precition corrected. Must first use \code{predcorrect()} if \code{loess.ypc = TRUE}
#' @param lambda numeric vector of length 3 specifying lambda values for each quantile
#' @param span numeric number between 0,1 specying smoothing paramter for loess prediction corrected
#' @param ... other arguments
#' @seealso \code{\link{observed}} \code{\link{simulated}} \code{\link{censoring}} \code{\link{predcorrect}} \code{\link{stratify}} \code{\link{binning}} \code{\link{vpcstats}}
#' @examples 
#' \dontrun{
#' 
#'  vpc <- observed(dataObs, y = DV, x = TIME) %>%
#'       simulated(dataSim, y = DV) %>%
#'       binless() %>%
#'       vpcstats()
#'       
#'  #Binless example with LOESS prediction correctection
#'  
#'  vpc <- observed(dataObs, y = DV, x = TIME) %>%
#'       simulated(dataSim, y = DV) %>%
#'       stratify(~ ISM) %>%
#'       predcorrect(pred = PRED) %>%
#'       binless(optimize = TRUE, loess.ypc = TRUE) %>%
#'       vpcstats()
#'       
#' # Binless example with user specified lambda values stratified on 
#' # "ISM" with 2 levels (0, 1), 5%, 50%, 95% quantiles.
#'  
#'  lambda_strat <- data.table(
#'  ISM0 = c(3,5,2),
#'  ISM1 = c(1, 3, 4),
#'  )
#'  
#'  vpc <- observed(dataObs, y = DV, x = TIME) %>%
#'       simulated(dataSim, y = DV) %>%
#'       stratify(~ ISM) %>%
#'       binless(qpred = c(0.1, 0.5, 0.9), optimize = FALSE, lambda = lambda_strat) %>%
#'       vpcstats()
#' }
#' @export 
binless <- function(o, ...) UseMethod("binless")

#' @rdname binless
#' @export
binless.tidyvpcobj <- function(o, qpred = c(0.05, 0.50, 0.95), optimize = TRUE, optimization.interval = c(0,7), conf.level = .95, loess.ypc = FALSE,  lambda = NULL, span = NULL, ...) {
  
  if(class(o) != "tidyvpcobj") {
    stop("No tidyvpcobj found, observed(...) %>% simulated(...) must be called prior to binless()")
  }
  
  if(!optimize && is.null(lambda)) {
    stop("Set optimize = TRUE if no lambda specified")
  }
  
  if(loess.ypc && is.null(o$predcor)) {
    stop("Use predcorrect() before binless() in order to use LOESS prediction corrected")
  }
  
  if(!is.null(span) && !loess.ypc) {
    stop("Set loess.ypc = TRUE and optimize = FALSE if setting span smoothing parameter for LOESS prediction corrected")
  }
  
  o %>%
    binlessaugment(qpred = qpred, interval =  optimization.interval, loess.ypc = loess.ypc) %>%
    binlessfit(conf.level = conf.level, llam.quant = lambda, span = span)
  
}

#' Prediction corrected Visual Predictive Check (pcVPC)
#' 
#' Specify prediction variable for pcVPC
#' 
#' @title predcorrect
#' @param o tidyvpc object
#' @param pred prediction variable in observed data 
#' @param data observed data supplied in \code{observed()} function
#' @param ... Other arguments to include
#' @param log logical indicating whether DV was modeled in logarithimic scale
#' @seealso \code{\link{observed}} \code{\link{simulated}} \code{\link{censoring}} \code{\link{stratify}} \code{\link{binning}} \code{\link{binless}} \code{\link{vpcstats}}

#' @export
predcorrect <- function(o, ...) UseMethod("predcorrect")

#' @rdname predcorrect
#' @export
predcorrect.tidyvpcobj <- function(o, pred, data=o$data, ..., log=FALSE) { 
  
  ypc <- y <- NULL
  
  if (missing(pred)) {
    pred <- o$pred
  } else {
    pred <- rlang::eval_tidy(rlang::enquo(pred), data)
  }
  if (is.null(pred)) {
    stop("No pred specified")
  }
  
  stratbin <- o$.stratbin #create loess predcorrect argument in function if want to use below stop because binless comes after predcorrect
  # if (is.null(stratbin)) {
  #     stop("Need to specify binning before pred correction. For binless method set argument loess.ypc = TRUE.")
  # }
  
  mpred <- data.table(stratbin, pred)
  mpred <- mpred[, mpred := median(pred), by=stratbin]
  mpred <- mpred$mpred
  
  if (log) {
    o$obs[, ypc := (mpred - pred) + y]
    o$sim[, ypc := (mpred - pred) + y]
  } else {
    o$obs[, ypc := (mpred/pred)*y]
    o$sim[, ypc := (mpred/pred)*y]
  }
  
  
  update(o, predcor=TRUE, predcor.log=log, pred=pred )
}

#' No pred correction for Visual Predictive Check (VPC)
#' 
#' Optional function to use indicating no pred correction for VPC. 
#' 
#' @title nopredcorrect
#' @param o tidyvpcobj
#' @param ... other arguments to include
#' @export
nopredcorrect <- function(o, ...) UseMethod("nopredcorrect")

#' @rdname nopredcorrect
#' @export
nopredcorrect.tidyvpcobj <- function(o, ...) {
  update(o, predcor=FALSE)
}

#' Compute VPC statistics
#' 
#' Compute predictional interval statistics for VPC
#' 
#' @title vpcstats
#' @param o tidyvpc object
#' @param qpred Numeric vector of length 3 specifying quantile prediction interval 
#' @param ... Other arguments to include
#' @param conf.level Numeric specifying confidence level
#' @param quantile.type Numeric indicating quantile type. See \code{\link[stats]{quantile}} 
#' @seealso \code{\link{observed}} \code{\link{simulated}} \code{\link{censoring}} \code{\link{stratify}} \code{\link{binning}} \code{\link{binless}} \code{\link{predcorrect}}
#' @export
vpcstats <- function(o, ...) UseMethod("vpcstats")

#' @rdname vpcstats
#' @export
vpcstats.tidyvpcobj <- function(o, qpred=c(0.05, 0.5, 0.95), ..., conf.level=0.95, quantile.type=7) {
  
  if(!is.null(o$rqss.obs.fits)) {
    .binlessvpcstats(o)
  } else {
    repl <- ypc <- y <- blq <- lloq <- alq <- uloq <- NULL
    . <- list
    
    obs      <- o$obs
    sim      <- o$sim
    predcor  <- o$predcor
    stratbin <- o$.stratbin
    xbin     <- o$xbin
    
    
    if (is.null(stratbin)) {
      stop("Need to specify binning before calling vpcstats.")
    }
    if (any(is.na(stratbin$bin))) {
      warning("There are bins missing. Has binning been specified for all strata?", call.=F)
    }
    
    .stratbinrepl <- data.table(stratbin, sim[, .(repl)])
    
    myquant1 <- function(y, probs, qname=paste0("q", probs), type=quantile.type, blq=F, alq=F) {
      y <- y + ifelse(blq, -Inf, 0) + ifelse(alq, Inf, 0)
      y <- quantile(y, probs=probs, type=type, names=F, na.rm=T)
      y[y == -Inf] <- NA
      data.frame(qname, y)
    }
    
    myquant2 <- function(y, probs, qname=paste0("q", probs), type=quantile.type) {
      y <- quantile(y, probs=probs, type=type, names=F, na.rm=T)
      setNames(as.list(y), qname)
    }
    
    if (isTRUE(predcor)) {
      qobs <- obs[, myquant1(ypc, probs=qpred, blq=blq, alq=alq), by=stratbin]
      qsim <- sim[, myquant1(ypc, probs=qpred, blq=F, alq=F),     by=.stratbinrepl]
    } else {
      qobs <- obs[, myquant1(y, probs=qpred, blq=blq, alq=alq), by=stratbin]
      qsim <- sim[, myquant1(y, probs=qpred, blq=F, alq=F),     by=.stratbinrepl]
    }
    
    .stratbinquant <- qsim[, !c("repl", "y")]
    qconf <- c(0, 0.5, 1) + c(1, 0, -1)*(1 - conf.level)/2
    qqsim <- qsim[, myquant2(y, probs=qconf, qname=c("lo", "md", "hi")), by=.stratbinquant]
    stats <- qobs[qqsim, on=names(.stratbinquant)]
    stats <- xbin[stats, on=names(stratbin)]
    setkeyv(stats, c(names(o$strat), "xbin"))
    
    if (!is.null(obs$blq) && any(obs$blq)) {
      sim[, lloq := rep(obs$lloq, len=.N)]
      sim[, blq := (y < lloq)]
      pctblqobs <- obs[, .(y=100*mean(blq)), by=stratbin]
      pctblqsim <- sim[, .(y=100*mean(blq)), by=.stratbinrepl]
      .stratbinpctblq <- pctblqsim[, !c("repl", "y")]
      qpctblqsim <- pctblqsim[, myquant2(y, probs=qconf, qname=c("lo", "md", "hi")), by=.stratbinpctblq]
      pctblq <- pctblqobs[qpctblqsim, on=names(.stratbinpctblq)]
      pctblq <- xbin[pctblq, on=names(stratbin)]
      setkeyv(pctblq, c(names(o$strat), "xbin"))
    } else {
      pctblq <- NULL
    }
    
    if (!is.null(obs$alq) && any(obs$alq)) {
      sim[, uloq := rep(obs$uloq, len=.N)]
      sim[, alq := (y > uloq)]
      pctalqobs <- obs[, .(y=100*mean(alq)), by=stratbin]
      pctalqsim <- sim[, .(y=100*mean(alq)), by=.stratbinrepl]
      .stratbinpctalq <- pctalqsim[, !c("repl", "y")]
      qpctalqsim <- pctalqsim[, myquant2(y, probs=qconf, qname=c("lo", "md", "hi")), by=.stratbinpctalq]
      pctalq <- pctalqobs[qpctalqsim, on=names(.stratbinpctalq)]
      pctalq <- xbin[pctalq, on=names(stratbin)]
      setkeyv(pctalq, c(names(o$strat), "xbin"))
    } else {
      pctalq <- NULL
    }

    update(o, stats=stats, pctblq=pctblq, pctalq=pctalq, conf.level=conf.level)
  }
}

#' @export
update.tidyvpcobj <- function(object, ...) {
  args <- list(...)
  for (i in names(args)) {
    object[[i]] <- args[[i, exact=TRUE]]
  }
  object
}

#' Print a \code{tidyvpcobj}.
#' @param x An object.
#' @param ... Further arguments can be specified but are ignored.
#' @return Returns \code{x} invisibly.
#' @export
print.tidyvpcobj <- function(x, ...) {
  if (!is.null(x$sim)) {
    nrep <- nrow(x$sim)/nrow(x$obs)
    if (isTRUE(x$predcor)) {
      cat("Prediction corrected ")
    }
    cat(sprintf("VPC with %d replicates", nrep), "\n")
  }
  cat(sprintf("Stratified by: %s", paste0(names(x$strat), collapse=", ")), "\n")
  if (!is.null(x$stats)) {
    print(x$stats)
  }
  invisible(x)
}

.binlessvpcstats <-  function(o, qpred=c(0.05, 0.5, 0.95), ..., conf.level=0.95, quantile.type=7){
  y <- x <- blq <- fit <- . <- repl <- cprop <- rqssmed <- llam.med <- c.rqssmed <-  NULL
  
  obs.fits <- o$rqss.obs.fits
  sim.fits <- o$rqss.sim.fits
  obs      <- o$obs
  sim      <- o$sim
  predcor  <- o$predcor
  xbinless <- o$obs$x
  
  if(!is.null(o$strat)) {
    stratx <- obs.fits[, list(x, o$strat)]
    x.binless <-  c("x", "qname", names(o$strat))
  } else {
    x.binless <- c("x", "qname")
  }
  
  qpred <- o$qpred
  conf.level <- o$conf.level
  qconf <- c(0, 0.5, 1) + c(1, 0, -1)*(1 - conf.level)/2
  
  if(!is.null(obs$blq) && any(obs$blq)) {
    if(!is.null(o$strat)) {
      stratlloq <- c(names(o$strat), "lloq")
      lloq <- obs[, stratlloq, with = FALSE] 
      lloq <- unique(lloq)
      obs.fits <- obs.fits[lloq, on = names(o$strat)]
    } else {
      obs.fits[, lloq := rep(obs$lloq, len=.N)]    
    }
    obs.fits[, blq := ifelse(fit < lloq, TRUE, FALSE)]
  }
  
  obs.fits <- setnames(obs.fits[, lapply(.SD, median), by = x.binless], "fit", "y")
  sim.fits <- setnames(sim.fits, c("fit", "fit.lb", "fit.ub"), c("md", "lo", "hi"))
  
  if(!is.null(obs$blq) && any(obs$blq)) {
    obs.fits[, blq := ifelse(y < lloq, TRUE, FALSE)]
    obs.fits[, y := ifelse(blq == TRUE, NA, y)]
  }
  
  if (!is.null(o$strat)) {
    stats <- obs.fits[sim.fits, on = c("x", "qname", names(o$strat))]
  } else {
    stats <- obs.fits[sim.fits, on = c("x", "qname")]
  }
  
  if (!is.null(obs$blq) && any(obs$blq)) {
    if(is.null(o$strat)) {
    sim[, lloq := rep(obs$lloq, len=.N)]
    sim[, blq := (y < lloq)]
    setorder(obs, cols = x)
    cprop.obs <- obs[, cprop := cumsum(blq) / 1:length(blq)]
    
    sic.cprop <- function(llam) {
      a <- AIC(
        rqss(
          cprop.obs$cprop ~ 
            qss(cprop.obs$x, lambda=exp(llam)), 
          tau=0.5, na.action=na.exclude
        ),
        k=-1
      )
    }
    llam.med.cprop <- optimize(sic.cprop, interval=c(0, 7))$min
    
    med.obs.cprop <- rqss(
      cprop.obs$cprop ~ qss(cprop.obs$x, lambda=exp(llam.med.cprop)), 
      tau=0.50
    )
    cprop.obs$med <- fitted(med.obs.cprop)
    
    setorder(sim, repl, x)[, cprop := cumsum(blq) / 1:length(blq), by=list(repl)]
    
    suppressWarnings(sim[, rqssmed := fitted(rqss(cprop ~ qss(x, lambda = exp(llam.med.cprop)),
                                                      tau = 0.5, na.action = na.exclude, .SD)), by = .(repl)])  
    
    u.x <- unique(cprop.obs$x) #%#
    med.obs.cprop <- tapply(cprop.obs$med, cprop.obs$x, median)
    med.sims.bql    <- tapply(sim$rqssmed, sim$x, median)
    med.sims.bql.lb <- tapply(sim$rqssmed, sim$x, quantile, probs=c(qconf[[1]]))
    med.sims.bql.ub <- tapply(sim$rqssmed, sim$x, quantile, probs=c(qconf[[3]]))
    pctblq <- data.table(cbind(u.x,med.obs.cprop, med.sims.bql.lb, med.sims.bql, med.sims.bql.ub))
    
    setnames(pctblq, c("x", "y", "lo", "md", "hi"))
    } else {
        strat <- o$strat
        strat.split <- split(obs, strat)
        x.strat <- c("x", names(strat))
        sim[, lloq := rep(obs$lloq, len=.N), by = names(strat)]
        sim[, blq := (y < lloq)]
        stratx.binless <- obs[, list(x, o$strat)]
        stratxrepl <- data.table(stratx.binless, sim[, .(repl)])
        #sim.strat <- sim[, c(names(strat)) := rep(strat, len = .N), by = .(repl)]
        strat.split.sim <- split(sim, strat)    
        sic.strat.cprop <- function(llam){
          a <- AIC(
            rqss(
              cprop ~
                qss(x, lambda=exp(llam)),
              tau=.5, na.action=na.exclude, data = strat.split[[i]]
            ),
            k=-1
          )
        }
        llam.strat.med.cprop <- vector("list", length(strat.split))
        for (i in seq_along(strat.split)) {
          setorder(strat.split[[i]], cols = x)
          strat.split[[i]] <- strat.split[[i]][, cprop := cumsum(blq) / 1:length(blq)]
          llam.strat.med.cprop[[i]]   <- strat.split[[i]][, .(llam.med = optimize(sic.strat.cprop,  interval=c(0, 7))$min)][,.(med = unlist(llam.med))]
          strat.split[[i]][, c.rqssmed := fitted(rqss(cprop ~ qss(x, lambda = exp(llam.strat.med.cprop[[i]][[1]])),tau= .5, na.action = na.exclude))]
        }
        
        obs.cprop <- rbindlist(strat.split)
        obs.cprop <- setnames(obs.cprop[, lapply(.SD, median, na.rm = TRUE), by = x.strat, .SDcols = "c.rqssmed"], "c.rqssmed", "y")
        
        for (i in seq_along(strat.split.sim)) {
          setorder(strat.split.sim[[i]], cols = repl, x)
          strat.split.sim[[i]] <- strat.split.sim[[i]][, cprop := cumsum(blq) / 1:length(blq), by = .(repl)]
          strat.split.sim[[i]][, c.rqssmed := fitted(rqss(cprop ~ qss(x, lambda = exp(llam.strat.med.cprop[[i]][[1]])),tau= .5, na.action = na.exclude)), by = .(repl)]
        }
        
        sim.cprop <- rbindlist(strat.split.sim)
        sim.med <- setnames(sim.cprop[, lapply(.SD, median, na.rm = TRUE), by = x.strat, .SDcols = "c.rqssmed"], "c.rqssmed", "md")
        sim.lb <- setnames(sim.cprop[, lapply(.SD, quantile, probs = qconf[[1]]), by = x.strat, .SDcols = "c.rqssmed"], "c.rqssmed", "lo")
        sim.ub <- setnames(sim.cprop[, lapply(.SD, quantile, probs = qconf[[3]]), by = x.strat, .SDcols = "c.rqssmed"], "c.rqssmed", "hi")
        
        pctblq <- obs.cprop[sim.med, on = x.strat]
        pctblq <- pctblq[sim.lb, on = x.strat]
        pctblq <- pctblq[sim.ub, on = x.strat]
    }
  } else {
    pctblq <- NULL
  }
  
  update(o, stats = stats, pctblq = pctblq)
}

#' Run Shiny app for tidyvpc
#' 
#' Use this function to run Shiny application to parameterize VPC from a GUI and generate corresponding tidyvpc code to derive VPC.
#' 
#' @title runShinyVPC
#' @seealso \href{https://github.com/jameswcraig/shiny-vpc/blob/master/README.md/}{Shiny-VPC GitHub}
#' @export

runShinyVPC <- function() {
 
  packagesCRAN <- c("remotes", "shiny", "backports", "DT", "ggplot2", "rlang", "shinyAce", "shinydashboard", "shinydashboardPlus", "shinyjs", "shinycssloaders", "shinyWidgets")
  if (length(setdiff(packagesCRAN, rownames(installed.packages()))) > 0) {
    install.packages(setdiff(packagesCRAN, rownames(installed.packages())))  
  }
  
  if(!is.element("shinymeta", installed.packages()[,1])) {
    if (requireNamespace("remotes", quietly = TRUE)) {
      remotes::install_github("rstudio/shinymeta")
    }
  }
  
  if (requireNamespace("shiny", quietly = TRUE)) {
  shiny::runGitHub("shiny-vpc", "jameswcraig")
  }
  
}

#' Obtain information about the bins from a VPC object.
#' @param o An object.
#' @param ... Additional arguments.
#' @return A `data.table` containing the following columns:
#' \itemize{
#'   \item \code{nobs}: the number of observed data points in the bin
#'   \item \code{xmedian}: the median x-value of the observed data points in the bin
#'   \item \code{xmean}: the mean x-value of the observed data points in the bin
#'   \item \code{xmax}: the maximum x-value of the observed data points in the bin
#'   \item \code{xmin}: the minimum x-value of the observed data points in the bin
#'   \item \code{xmid}: the value halfway between `xmin` and `xmax`.
#'   x-value of the observed data points in the bin
#'   \item \code{xleft}: the value halfway between the minimum x-value of the
#'   current bin and the maximum x-value of the previous bin to the left (for
#'   the left-most bin it is the minimum x-value).
#'   \item \code{xright}: the value halfway between the maximum x-value of the
#'   current bin and the minimum x-value of the next bin to the right (for the
#'   right-most bin it is the maximum x-value).
#'   \item \code{xcenter}: the value halfway between `xleft` and `xright`.
#' }
#' In addition, if statification was performed, the stratification columns will
#' be included as well.
#' @export
bininfo <- function(o, ...) UseMethod("bininfo")

#' @describeIn bininfo Method for \code{tidyvpcobj}.
#' @param by.strata Should the calculations be done by strata? Defaults to what
#'  was specified when the binning was done.
#' @export
bininfo.tidyvpcobj <- function(o, by.strata=o$bin.by.strata, ...) {
  
  x <- xmin <- xmax <- bin <- NULL
  
  f1 <- function(x) {
    nobs    <- sum(!is.na(x))
    xmedian <- median(x, na.rm=T)
    xmean   <- mean(x, na.rm=T)
    xmin    <- min(x, na.rm=T)
    xmax    <- max(x, na.rm=T)
    xmid    <- 0.5*(xmin + xmax)
    data.table(nobs, xmedian, xmean, xmin, xmax, xmid)
  }
  
  # Compute xleft and xright
  f2 <- function(xmin, xmax) {
    xmin    <- c(xmin, xmax[length(xmax)])
    xmax    <- c(xmin[1], xmax)
    breaks  <- 0.5*(xmin + xmax)
    xleft   <- breaks[-length(breaks)]
    xright  <- breaks[-1]
    xcenter <- 0.5*(xleft + xright)
    data.table(xleft, xright, xcenter)
  }
  if (by.strata && !is.null(o$strat)) {
    bi <- o$obs[, f1(x), by=o$.stratbin]
    setkeyv(bi, c(names(o$strat), "xmin"))
    bi[, c(.SD, f2(xmin, xmax)), by=names(o$strat)]
  } else {
    bi <- o$obs[, f1(x), by=bin]
    setkeyv(bi, "xmin")
    bi <- cbind(bi, bi[, f2(xmin, xmax)])
    bi <- bi[unique(o$.stratbin), on="bin"]
    setkeyv(bi, "xmin")
    bi[, c(names(o$.stratbin), setdiff(names(bi), names(o$.stratbin))), with=F]
  }
}

#' Different functions that perform binning.
#'
#' @param breaks A numeric vector of values that designate cut points between bins.
#' @param centers A numeric vector of values that designate the center of each bin.
#' @param nbins The number of bins to split the data into.
#' @param style a binning style (see ?classInt::classIntervals for details).
#' @return Each of these functions returns a function of a single numeric
#' vector `x` that assigns each value of `x` to a bin.
#' @examples
#'
#' x <- c(rnorm(10, 1, 1), rnorm(10, 3, 2), rnorm(20, 5, 3))
#' centers <- c(1, 3, 5)
#' nearest(centers)(x)
#'
#' breaks <- c(2, 4)
#' cut_at(breaks)(x)
#'
#' bin_by_eqcut(nbins=4)(x)
#' bin_by_ntile(nbins=4)(x)
#'
#' \dontrun{
#' bin_by_pam(nbins=4)(x)
#' bin_by_classInt("pretty", nbins=4)(x)
#' }
#'
#' @name binningfunctions
NULL

#' @rdname binningfunctions
#' @export
cut_at <- function(breaks) {
  breaks <- .check_breaks(breaks)
  function(x, ..., right=F) {
    breaks <- .resolve_breaks(breaks, ...)
    breaks <- sort(unique(breaks))
    if (min(x) < min(breaks)) {
      breaks <- c(min(x), breaks)
    }
    if (max(x) > max(breaks)) {
      breaks <- c(breaks, max(x))
    }
    as.character(cut(x, breaks, include.lowest=T, right=right))
  }
}

#' @rdname binningfunctions
#' @export
nearest <- function(centers) {
  centers <- .check_centers(centers)
  function(x, ...) {
    centers <- .resolve_centers(centers, ...)
    centers <- sort(unique(centers))
    dist <- function(a, b) abs(a - b)
    d <- outer(x, centers, dist)
    m <- apply(d, 1, which.min)
    centers[m]
  }
}

#' @rdname binningfunctions
#' @export
bin_by_ntile <- function(nbins) {
  nbins <- .check_nbins(nbins)
  function(x, ...) {
    nbins <- .resolve_nbins(nbins, ...)
    
    # Mimic the function from dplyr
    len <- sum(!is.na(x))
    r <- rank(x, ties.method="first", na.last="keep")
    as.integer(floor(nbins*(r - 1)/len + 1))
  }
}

#' @rdname binningfunctions
#' @export
bin_by_eqcut <- function(nbins) {
  nbins <- .check_nbins(nbins)
  function(x, ..., quantile.type=7) {
    nbins <- .resolve_nbins(nbins, ...)
    
    # Mimic the function from table1
    breaks <- quantile(x, probs=seq.int(nbins - 1)/nbins, na.rm=T, type=quantile.type)
    cut_at(breaks)(x)
  }
}

#' @rdname binningfunctions
#' @export
bin_by_pam <- function(nbins) {
  has_cluster <- requireNamespace("cluster", quietly=TRUE)
  if (!has_cluster) {
    stop("Package 'cluster' is required to use the binning method. Please install it.")
  }
  nbins <- .check_nbins(nbins)
  function(x, ...) {
    nbins <- .resolve_nbins(nbins, ...)
    
    centers <- sort(cluster::pam(x, nbins)$medoids)
    nearest(centers)(x)
  }
}

#' @rdname binningfunctions
#' @export
bin_by_classInt <- function(style, nbins=NULL) {
  has_classInt <- requireNamespace("classInt", quietly=TRUE)
  if (!has_classInt) {
    stop("Package 'classInt' is required to use the binning method. Please install it.")
  }
  style <- style
  if (!is.null(nbins)) {
    nbins <- .check_nbins(nbins)
  }
  function(x, ...) {
    args <- list(var=x, style=style)
    if (!is.null(nbins)) {
      nbins <- .resolve_nbins(nbins, ...)
      args$n <- nbins
    }
    args <- c(args, list(...))
    if (style %in% c("kmeans", "hclust", "dpih")) {
      # These don't accept '...' arguments
      args1 <- args[intersect(names(args), methods::formalArgs(classInt::classIntervals))]
      args2 <- if (style == "kmeans") {
        args[intersect(names(args), methods::formalArgs(stats::kmeans))]
      } else if (style == "hclust") {
        args[intersect(names(args), methods::formalArgs(stats::hclust))]
      } else if (style == "dpih") {
        has_KernSmooth <- requireNamespace("KernSmooth", quietly=TRUE)
        if (!has_KernSmooth) {
          stop("Package 'KernSmooth' is required to use the binning method. Please install it.")
        }
        args[intersect(names(args), methods::formalArgs(KernSmooth::dpih))]
      } else {
        list()
      }
      args <- c(args1, args2)
    }
    args <- args[!duplicated(args)]
    breaks <- do.call(classInt::classIntervals, args)$brks
    cut_at(breaks)(x)
  }
}

#' Perform a consistency check on observed and simulated data.
#' 
#' This function performs a simple consistency check on an observed and
#' simulated dataset to make sure they are consistent with respect to ordering
#' as required by the other functions used in the VPC calculation.
#'
#' The consistency check is performed by comparing a combination of unique
#' subject identifier (ID) and time. Both `data.frame`s must be given with
#' those in positions 1 and 2 repectively.
#'
#' @param obs,sim A `data.frame` with 2 columns (see Details).
#' @param tol A tolerance for comparing time values.
#' @return The number of replicates contained in `sim`.
#' @seealso \code{\link{observed}}, \code{\link{simulated}}.
#' @examples
#'
#' \dontrun{
#' library(vpc)
#'
#' exampleobs <- as.data.table(vpc::simple_data$obs)[MDV == 0]
#' examplesim <- as.data.table(vpc::simple_data$sim)[MDV == 0]
#'
#' check_order(exampleobs[, .(ID, TIME)], examplesim[, .(ID, TIME)])
#' }
#' @export
check_order <- function(obs, sim, tol=1e-5) {
  if (nrow(sim) %% nrow(obs) != 0) {
    stop("Rows in sim is not a multiple of rows in obs")
  }
  if (is.numeric(obs[[1]])) obs[[1]] <- as.numeric(obs[[1]])
  if (is.numeric(sim[[1]])) sim[[1]] <- as.numeric(sim[[1]])
  if (is.factor(obs[[1]])) obs[[1]] <- as.character(obs[[1]])
  if (is.factor(sim[[1]])) sim[[1]] <- as.character(sim[[1]])
  if (!identical(rep(obs[[1]], len=nrow(sim)), sim[[1]])) {
    stop("ID columns are not identical")
  } 
  if (!all(abs(rep(obs[[2]], len=nrow(sim)) - sim[[2]]) < tol)) {
    stop("Time columns are not equal")
  } 
  nrow(sim) / nrow(obs)
}


# Internal function
.check_centers <- function(centers) {
  if (is.data.frame(centers)) {
    centers <- as.data.table(centers)
    if (is.null(centers$centers)) {
      stop("centers data.frame must contain column centers")
    }
    if (any(is.na(centers$centers))) {
      stop("centers cannot contain missing values")
    }
    keycols <- setdiff(names(centers), "centers")
    setkeyv(centers, keycols)
  } else if (is.numeric(centers)) {
    if (any(is.na(centers))) {
      stop("centers cannot contain missing values")
    }
  } else {
    stop("centers must be a numeric vector or data.frame")
  }
  centers
}

# Internal function
.resolve_centers <- function(centers, ...) {
  if (is.data.table(centers)) {
    keycols <- key(centers)
    key <- as.data.table(list(...)[keycols])
    centers <- unique(centers[key]$centers)
  }
  if (is.null(centers) || !is.numeric(centers) || any(is.na(centers))) {
    stop("invalid centers")
  }
  centers
}

# Internal function
.check_breaks <- function(breaks) {
  if (is.data.frame(breaks)) {
    breaks <- as.data.table(breaks)
    if (is.null(breaks$breaks)) {
      stop("breaks data.frame must contain column breaks")
    }
    if (any(is.na(breaks$breaks))) {
      stop("breaks cannot contain missing values")
    }
    keycols <- setdiff(names(breaks), "breaks")
    setkeyv(breaks, keycols)
  } else if (is.numeric(breaks)) {
    if (any(is.na(breaks))) {
      stop("breaks cannot contain missing values")
    }
  } else {
    stop("breaks must be a numeric vector or data.frame")
  }
  breaks
}

# Internal function
.resolve_breaks <- function(breaks, ...) {
  if (is.data.table(breaks)) {
    keycols <- key(breaks)
    key <- as.data.table(list(...)[keycols])
    breaks <- breaks[key]$breaks
  }
  if (is.null(breaks) || !is.numeric(breaks) || any(is.na(breaks))) {
    stop("invalid breaks")
  }
  breaks
}

# Internal function
.check_nbins <- function(nbins) {
  if (is.data.frame(nbins)) {
    nbins <- as.data.table(nbins)
    if (is.null(nbins$nbins)) {
      stop("nbins data.frame must contain column nbins")
    }
    if (any(is.na(nbins$nbins))) {
      stop("nbins cannot contain missing values")
    }
    keycols <- setdiff(names(nbins), "nbins")
    setkeyv(nbins, keycols)
  } else if (is.numeric(nbins) && length(nbins) == 1) {
    if (any(is.na(nbins))) {
      stop("nbins cannot contain missing values")
    }
  } else {
    stop("nbins must be a numeric vector of length 1 or data.frame")
  }
  nbins
}

# Internal function
.resolve_nbins <- function(nbins, ...) {
  if (is.data.table(nbins)) {
    keycols <- key(nbins)
    key <- as.data.table(list(...)[keycols])
    nbins <- unique(nbins[key]$nbins)
  }
  if (is.null(nbins) || !(is.numeric(nbins) && length(nbins) == 1 && !is.na(nbins))) {
    stop("nbins must be uniquely determined")
  }
  nbins
}


binlessaugment <- function(o, qpred = c(0.05, 0.50, 0.95), interval = c(0,7), loess.ypc = FALSE, ...) { 
  l.ypc <- strat.split <- y <- NULL
  
  qpred <- sort(qpred)
  obs <- o$obs
  log <- o$predcor.log
  
  environment(.autoloess) <- environment()
  
  if (loess.ypc) {  #Split data on strata to optimize loess
    if (!is.null(o$strat)) {
      pred <- o$pred
      obs <- cbind(obs, pred)
      strat <- o$strat
      strat.split <- split(obs, strat)
      loess.mod.strat <- vector("list", length(strat.split))
      names(loess.mod.strat) <- names(strat.split)
      if(isTRUE(o$predcor.log)) {
        for (i in seq_along(strat.split)) {
          loess.mod.strat[[i]] <-  .autoloess(loess(pred ~ x, span = .5, data = strat.split[[i]]))
          strat.split[[i]][, lpred := fitted(loess(pred ~ x, span = loess.mod.strat[[i]]$span, na.action = na.exclude))]
          strat.split[[i]][, l.ypc := (lpred - pred) + y]
        }
      } else {
        for (i in seq_along(strat.split)) {
          loess.mod.strat[[i]] <-  .autoloess(loess(pred ~ x, span = .5, data = strat.split[[i]]))
          strat.split[[i]][, lpred := fitted(loess(pred ~ x, span = loess.mod.strat[[i]]$span, na.action = na.exclude))]
          strat.split[[i]][, l.ypc := (lpred/pred) * y]
        }
      }
      span <- .getspan(loess.mod.strat)
    } else {
      pred <- o$pred
      obs <- cbind(obs, pred)
      loess.mod <-  .autoloess(loess(pred ~ x, span = .5, data = obs))
      lpred <- fitted(loess.mod$fit)
      span <- loess.mod$span
      if (isTRUE(o$predcor.log)) {
        obs[, l.ypc := (lpred - pred) + y]
      } else {
        obs[, l.ypc := (lpred/pred) * y]
      }
    }
  }
  
  if(!loess.ypc && !is.null(o$strat)) {
    strat <- o$strat
    strat.split <- split(obs, strat)
  }
  
  # Internal Function
  .sic.strat.ypc <- function(llam, quant) {
    a <- AIC(
      rqss(
        l.ypc ~
          qss(x, lambda=exp(llam)),
        tau=quant, na.action=na.exclude, data = strat.split[[i]]
      ),
      k=-1
    )
  }
  .sic.strat <- function(llam, quant){
    a <- AIC(
      rqss(
        y ~
          qss(x, lambda=exp(llam)),
        tau=quant, na.action=na.exclude, data = strat.split[[i]]
      ),
      k=-1
    )
  }
  
  .sic.ypc <- function(llam, quant){
    a <- AIC(
      rqss(
        l.ypc ~
          qss(x, lambda=exp(llam)),
        tau=quant, na.action=na.exclude, data = obs
      ),
      k=-1
    )
  }
  
  .sic <- function(llam, quant){
    a <- AIC(
      rqss(
        y ~
          qss(x, lambda=exp(llam)),
        tau=quant, na.action=na.exclude, data = obs
      ),
      k=-1
    )
  }
  
  if(loess.ypc) {
    if(!is.null(o$strat)){
      llamoptimize <- .sic.strat.ypc
    } else {
      llamoptimize  <- .sic.ypc
    }
  }
  
  if(!loess.ypc) {
    span <- NULL
    if(!is.null(o$strat)) {
      llamoptimize <- .sic.strat
    } else {
      llamoptimize <- .sic
    }
  }
  
  environment(llamoptimize) <- environment()
  . <- list
  if(!is.null(o$strat.split)) {
    llam.strat.lo  <- vector("list", length(strat.split))
    llam.strat.med <- vector("list", length(strat.split))
    llam.strat.hi  <- vector("list", length(strat.split))
    
    for (i in seq_along(strat.split)) {
      llam.strat.lo[[i]]    <- strat.split[[i]][, .(llam.lo  = optimize(llamoptimize, quant = qpred[1], interval = interval)$min)][,.(lo = unlist(llam.lo))]
      names(llam.strat.lo)  <- names(strat.split)
      setnames(llam.strat.lo[[i]], paste0("q", qpred[1]))
      llam.strat.med[[i]]   <- strat.split[[i]][, .(llam.med = optimize(llamoptimize, quant = qpred[2], interval = interval)$min)][,.(med = unlist(llam.med))]
      names(llam.strat.med) <- names(strat.split)
      setnames(llam.strat.med[[i]], paste0("q", qpred[2]))
      llam.strat.hi[[i]]    <- strat.split[[i]][, .(llam.hi  = optimize(llamoptimize, quant = qpred[3], interval = interval)$min)][,.(hi = unlist(llam.hi))]
      names(llam.strat.hi)  <- names(strat.split)
      setnames(llam.strat.hi[[i]], paste0("q", qpred[3]))
    }
    
    llam.qpred <- cbind(list(llam.strat.lo, llam.strat.med, llam.strat.hi))
    names(llam.qpred) <- paste0("q", qpred)
  } else {
    llam.lo  <- obs[, .(llam.lo = optimize(llamoptimize, quant = qpred[1], interval = interval)$min)]
    llam.med <- obs[, .(llam.med = optimize(llamoptimize, quant = qpred[2], interval = interval)$min)]
    llam.hi  <- obs[, .(llam.hi = optimize(llamoptimize, quant = qpred[3], interval = interval)$min)]
    
    llam.qpred <- c(llam.lo, llam.med, llam.hi)
    names(llam.qpred) <- paste0("q", qpred)
    llam.qpred <- unlist(llam.qpred)
  }
  
  update(o, llam.qpred = llam.qpred, span = span, qpred = qpred, loess.ypc = loess.ypc)
}


binlessfit <- function(o, conf.level = .95, llam.quant = NULL, span = NULL, ...){
  y <- l.ypc <- repl <- NULL  
  . <- list
  
  qpred <- o$qpred
  qnames <- paste0("q", as.character(qpred))
  qconf <- c(0, 0.5, 1) + c(1, 0, -1)*(1 - conf.level)/2
  
  obs <- o$obs
  sim <- o$sim
  
  if(isTRUE(o$loess.ypc)) {
    pred <- o$pred
    obs <- cbind(obs, pred)
    sim[, pred := rep(pred, len=.N), by = .(repl)]
    if(is.null(span)) {
      span <- o$span  
    }
  }
  getllam <- function(qnames, userllam, stratlev) {
    userllam <- as.data.frame(userllam)
    userllam <- userllam[, order(names(userllam))]
    llam.list <- vector("list", length(qnames))
    names(llam.list) <- qnames
    if(stratlev == 2) {
      for (i in seq_along(llam.list)) {
        llam.list[[i]] <- list(lambda = userllam[i, 1], lambda = userllam[i,2])
      }
    }
    if(stratlev == 3) {
      for (i in seq_along(llam.list)) {
        llam.list[[i]] <- list(lambda = userllam[i, 1], lambda = userllam[i,2], lambda = userllam[i,3])
      }
    }
    if(stratlev == 4) {
      for (i in seq_along(llam.list)) {
        llam.list[[i]] <- list(lambda = userllam[i, 1], lambda = userllam[i,2], lambda = userllam[i,3], lambda = userllam[i,4])
      }
    }
    if(stratlev == 5) {
      for (i in seq_along(llam.list)) {
        llam.list[[i]] <- list(lambda = userllam[i, 1], lambda = userllam[i,2], lambda = userllam[i,3], lambda = userllam[i,4], lambda = userllam[i,5])
      }
    }
    names(llam.list[[1]]) <- names(o$strat.split)
    names(llam.list[[2]]) <- names(o$strat.split)
    names(llam.list[[3]]) <- names(o$strat.split)
    return(llam.list)
  }
  
  if(is.null(llam.quant)) {
    if(is.null(o$llam.qpred)) {
      stop("Must specify llambda for binlessfit. Include binlessaugment() before running binlessfit() for optimized llambda values using AIC.")
    } else {
      llam.qpred <- o$llam.qpred
    }
  } else if(!is.null(llam.quant) && !is.null(o$strat)) {
    stratlev <- lapply(o$strat, unique) 
    stratlev <- length(stratlev[[1]])
    #environment(.getllam) <- environment()
    llam.qpred <- getllam(qnames, llam.quant, stratlev)
  } else { 
    llam.qpred <- llam.quant
  } 
  
  if(is.null(span)) {
    if(!is.null(o$span) && isTRUE(o$loess.ypc)) {
      span <- o$span
    } else {
      span <- o$span
    }
  } 
  
  if(!is.null(o$strat)) {
    strat <- o$strat
    strat.split <- split(obs, strat)
    x.strat <- c("x", names(strat))
    sim.strat <- sim[, c(names(strat)) := rep(strat, len = .N), by = .(repl)]
    strat.split.sim <- split(sim, strat)
  }
  
  if(isTRUE(o$loess.ypc) && !is.null(o$strat)) {
    if(isTRUE(o$predcor.log)) {
      for(i in seq_along(strat.split)) {
        strat.split[[i]][, l.ypc := y +  (fitted(loess(pred ~ x, span = span[[i]], na.action = na.exclude, .SD)) - pred)]
      }
    } else {
      for(i in seq_along(strat.split)) {
        strat.split[[i]][, l.ypc := y * (fitted(loess(pred ~ x, span = span[[i]], na.action = na.exclude, .SD)) / pred)]
      } 
    }
    obs <- rbindlist(strat.split)
    o <- update(o, obs = obs)
  }
  
  if(isTRUE(o$loess.ypc) && is.null(o$strat)) {
    if(isTRUE(o$predcor.log)) {
      obs[, l.ypc := y + (fitted(loess(pred ~ x, span = span, na.action = na.exclude, .SD)) - pred)]
    } else {
      obs[, l.ypc := y * (fitted(loess(pred ~ x, span = span, na.action = na.exclude, .SD)) / pred)]
    }
    o <- update(o, obs = obs)
  }
  
  if (is.null(o$strat)) {
    if (isTRUE(o$loess.ypc)) {
      rqss.obs.fits <- .fitobs(obs, llam.qpred, qpred, l.ypc = TRUE)
      if(isTRUE(o$predcor.log)) {
        rqss.sim.fits <- .fitsim(sim, llam.qpred, span, qpred, qconf, l.ypc = TRUE, log = TRUE)
      } else {
        rqss.sim.fits <- .fitsim(sim, llam.qpred, span, qpred, qconf, l.ypc = TRUE)
      }
    } else {
      rqss.obs.fits <- .fitobs(obs, llam.qpred, qpred)
      rqss.sim.fits <- .fitsim(sim, llam.qpred, qpred = qpred, qconf= qconf)
    }
  } 
  
  
  if(!is.null(o$strat)){
    if(isTRUE(o$loess.ypc)){
      if(isTRUE(o$predcor.log)) {
        rqss.obs.fits <- .fitobs.strat(strat.split = strat.split, x.strat = x.strat, llam.qpred = llam.qpred, qpred = qpred, l.ypc = TRUE)
        rqss.sim.fits <- .fitsim.strat(strat.split.sim = strat.split.sim, x.strat = x.strat, llam.qpred = llam.qpred, span = span, qpred = qpred, qconf = qconf, l.ypc = TRUE, log = TRUE)
      } else {
        rqss.obs.fits <- .fitobs.strat(strat.split = strat.split, x.strat = x.strat, llam.qpred = llam.qpred, qpred = qpred, l.ypc = TRUE)
        rqss.sim.fits <- .fitsim.strat(strat.split.sim = strat.split.sim, x.strat = x.strat, llam.qpred = llam.qpred, span = span, qpred = qpred, qconf = qconf, l.ypc = TRUE)
      }
    } else {
      rqss.obs.fits <- .fitobs.strat(strat.split = strat.split, x.strat = x.strat, llam.qpred = llam.qpred, qpred = qpred)
      rqss.sim.fits <- .fitsim.strat(strat.split.sim = strat.split.sim, x.strat = x.strat, llam.qpred = llam.qpred, qpred = qpred, qconf = qconf)
    }
  }
  
  update(o, rqss.obs.fits = rqss.obs.fits, rqss.sim.fits = rqss.sim.fits, llam.qpred = llam.qpred, span = span, conf.level = conf.level)
  
}

# Internal Function
#Below function used for fitting rqss
.fitobs <- function(obs, llam.qpred, qpred, l.ypc = FALSE) {
  rqsslo <- rqssmed <- rqsshi <- NULL
  qnames <- paste0("q", as.character(qpred))
  
  if(l.ypc) {
    obs[, rqsslo := fitted(rqss(l.ypc ~ qss(x, lambda = exp(llam.qpred[[1]])), tau = qpred[1], na.action = na.exclude))]
    obs[, rqssmed := fitted(rqss(l.ypc ~ qss(x, lambda = exp(llam.qpred[[2]])), tau = qpred[2], na.action = na.exclude))]
    obs[, rqsshi := fitted(rqss(l.ypc ~ qss(x, lambda = exp(llam.qpred[[3]])), tau = qpred[3], na.action = na.exclude))]
    setnames(obs, c("rqsslo", "rqssmed", "rqsshi"), qnames)
    obs.fits <- melt(obs, id.vars = "x", measure.vars = qnames)
    obs.fits <- setnames(obs.fits, c("variable", "value"), c("qname", "fit"))
  } else {
    obs[, rqsslo := fitted(rqss(y ~ qss(x, lambda = exp(llam.qpred[[1]])), tau = qpred[1], na.action = na.exclude))]
    obs[, rqssmed := fitted(rqss(y ~ qss(x, lambda = exp(llam.qpred[[2]])), tau = qpred[2], na.action = na.exclude))]
    obs[, rqsshi := fitted(rqss(y ~ qss(x, lambda = exp(llam.qpred[[3]])), tau = qpred[3], na.action = na.exclude))]
    #obs.fits <- data.table(cbind(x,fitted(lo),fitted(med),fitted(hi)))
    setnames(obs, c("rqsslo", "rqssmed", "rqsshi"), qnames)
    obs.fits <- melt(obs, id.vars = "x", measure.vars = qnames)
    obs.fits <- setnames(obs.fits, c("variable", "value"), c("qname", "fit"))
  }
  return(obs.fits)
}

# Internal Function
.fitobs.strat <- function(strat.split, x.strat, llam.qpred, qpred, l.ypc = FALSE) {
  rqsslo <- rqssmed <- rqsshi <- NULL
  qnames <- paste0("q", as.character(qpred))
  
  if(l.ypc) {
    for (i in seq_along(strat.split)) {
      strat.split[[i]][, rqsslo  := fitted(rqss(l.ypc ~ qss(x, lambda = exp(llam.qpred[[1]][[i]][[1]])),tau= qpred[1], na.action = na.exclude))]
      strat.split[[i]][, rqssmed := fitted(rqss(l.ypc ~ qss(x, lambda = exp(llam.qpred[[2]][[i]][[1]])),tau= qpred[2], na.action = na.exclude))]
      strat.split[[i]][, rqsshi  := fitted(rqss(l.ypc ~ qss(x, lambda = exp(llam.qpred[[3]][[i]][[1]])),tau= qpred[3], na.action = na.exclude))]
    }
  } else {
    for (i in seq_along(strat.split)) {
      strat.split[[i]][, rqsslo  := fitted(rqss(y ~ qss(x, lambda = exp(llam.qpred[[1]][[i]][[1]])),tau= qpred[1], na.action = na.exclude))]
      strat.split[[i]][, rqssmed := fitted(rqss(y ~ qss(x, lambda = exp(llam.qpred[[2]][[i]][[1]])),tau= qpred[2], na.action = na.exclude))]
      strat.split[[i]][, rqsshi  := fitted(rqss(y ~ qss(x, lambda = exp(llam.qpred[[3]][[i]][[1]])),tau= qpred[3], na.action = na.exclude))]
    }
  }
  strat.obs.fits <- rbindlist(strat.split)
  strat.obs.fits <- setnames(strat.obs.fits, c("rqsslo", "rqssmed", "rqsshi"), qnames)
  strat.obs.fits <- melt(strat.obs.fits, id.vars = x.strat, measure.vars = qnames)
  strat.obs.fits <- setnames(strat.obs.fits, c("variable", "value"), c("qname", "fit"))
}

# Internal Function
.fitsim <- function(sim, llam.qpred, span = NULL, qpred, qconf, l.ypc = FALSE, log = FALSE) {
  . <- list
  rqsslo <- rqssmed <- rqsshi <- y <- pred <- repl <- NULL
  qnames <- paste0("q", as.character(qpred))
  
  if(l.ypc) {
    if(log) {
      sim[, l.ypc := y + (fitted(loess(pred ~ x, span = span, na.action = na.exclude, .SD)) - pred), by = .(repl)]
    } else {
      sim[, l.ypc := y * (fitted(loess(pred ~ x, span = span, na.action = na.exclude, .SD)) / pred), by = .(repl)]
    }
    suppressWarnings(sim[, rqsslo := fitted(rqss(l.ypc ~ qss(x, lambda = exp(llam.qpred[[1]])),
                                                 tau  = qpred[1], na.action = na.exclude, .SD)), by = .(repl)]) 
    suppressWarnings(sim[, rqssmed := fitted(rqss(l.ypc ~ qss(x, lambda = exp(llam.qpred[[2]])),
                                                  tau = qpred[2], na.action = na.exclude, .SD)), by = .(repl)]) 
    suppressWarnings(sim[, rqsshi := fitted(rqss(l.ypc ~ qss(x, lambda = exp(llam.qpred[[3]])),
                                                 tau  = qpred[3], na.action = na.exclude, .SD)), by = .(repl)])
    setnames(sim, c("rqsslo", "rqssmed", "rqsshi"), qnames)
  } else {
    suppressWarnings(sim[, rqsslo := fitted(rqss(y ~ qss(x, lambda = exp(llam.qpred[[1]])),
                                                 tau  = qpred[1], na.action = na.exclude, .SD)), by = .(repl)]) 
    suppressWarnings(sim[, rqssmed := fitted(rqss(y ~ qss(x, lambda = exp(llam.qpred[[2]])),
                                                  tau = qpred[2], na.action = na.exclude, .SD)), by = .(repl)]) 
    suppressWarnings(sim[, rqsshi := fitted(rqss(y ~ qss(x, lambda = exp(llam.qpred[[3]])),
                                                 tau  = qpred[3], na.action = na.exclude, .SD)), by = .(repl)])
    setnames(sim, c("rqsslo", "rqssmed", "rqsshi"), qnames)
  }
  
  sim.lb <- sim[, lapply(.SD, quantile, probs = qconf[[1]]), by = "x"] #CI lb
  sim.lb <- setnames(melt(sim.lb, id.vars = "x", measure.vars = qnames), c("x", "qname", "fit.lb"))  #Wide to long
  
  sim.ub <- sim[, lapply(.SD, quantile, probs = qconf[[3]]), by = "x"] #CI ub
  sim.ub <- setnames(melt(sim.ub, id.vars = "x", measure.vars = qnames), c("x", "qname", "fit.ub"))   #Wide to long
  
  sim <-  sim[, lapply(.SD, median, na.rm = TRUE), by = "x"] #Med fits
  sim <- setnames(melt(sim, id.vars = "x", measure.vars = qnames), c("x", "qname", "fit")) #wide to long
  
  sim <- sim[sim.lb, on = c("x", "qname")]
  sim <- sim[sim.ub, on = c("x", "qname")]
}

# Internal Function
.fitsim.strat <- function(strat.split.sim, x.strat, llam.qpred, span = NULL, qpred, qconf, l.ypc = FALSE, log = FALSE) {
  . <- list
  rqsslo <- rqssmed <- rqsshi <- y <- pred <- repl <- NULL
  qnames <- paste0("q", as.character(qpred))
  
  if(l.ypc) {
    if(log) {
      for (i in seq_along(strat.split.sim)) {
        strat.split.sim[[i]][, l.ypc := y + (fitted(loess(pred ~ x, span = span[[i]], na.action = na.exclude, .SD)) - pred), by = .(repl)]
        
        suppressWarnings(strat.split.sim[[i]][, rqsslo := fitted(rqss(l.ypc ~ qss(x, lambda = exp(llam.qpred[[1]][[i]][[1]])),
                                                                      tau  = qpred[1], na.action = na.exclude, .SD)), by = .(repl)][,.(fit.lo = unlist(rqsslo))])
        
        suppressWarnings(strat.split.sim[[i]][, rqssmed := fitted(rqss(l.ypc ~ qss(x, lambda = exp(llam.qpred[[2]][[i]][[1]])),
                                                                       tau = qpred[2], na.action = na.exclude, .SD)), by = .(repl)][,.(fit.med = unlist(rqssmed))])
        
        suppressWarnings(strat.split.sim[[i]][, rqsshi := fitted(rqss(l.ypc ~ qss(x, lambda = exp(llam.qpred[[3]][[i]][[1]])),
                                                                      tau  = qpred[3], na.action = na.exclude, .SD)), by = .(repl)][,.(fit.hi = unlist(rqsshi))])
      }
    } else {
      for (i in seq_along(strat.split.sim)) {
        strat.split.sim[[i]][, l.ypc := y * (fitted(loess(pred ~ x, span = span[[i]], na.action = na.exclude, .SD)) / pred), by = .(repl)]
        
        suppressWarnings(strat.split.sim[[i]][, rqsslo := fitted(rqss(l.ypc ~ qss(x, lambda = exp(llam.qpred[[1]][[i]][[1]])),
                                                                      tau  = qpred[1], na.action = na.exclude, .SD)), by = .(repl)][,.(fit.lo = unlist(rqsslo))])
        
        suppressWarnings(strat.split.sim[[i]][, rqssmed := fitted(rqss(l.ypc ~ qss(x, lambda = exp(llam.qpred[[2]][[i]][[1]])),
                                                                       tau = qpred[2], na.action = na.exclude, .SD)), by = .(repl)][,.(fit.med = unlist(rqssmed))])
        
        suppressWarnings(strat.split.sim[[i]][, rqsshi := fitted(rqss(l.ypc ~ qss(x, lambda = exp(llam.qpred[[3]][[i]][[1]])),
                                                                      tau  = qpred[3], na.action = na.exclude, .SD)), by = .(repl)][,.(fit.hi = unlist(rqsshi))])
      }
    }
  } else {  
    for (i in seq_along(strat.split.sim)) {
      suppressWarnings(strat.split.sim[[i]][, rqsslo := fitted(rqss(y ~ qss(x, lambda = exp(llam.qpred[[1]][[i]][[1]])),
                                                                    tau  = qpred[1], na.action = na.exclude, .SD)), by = .(repl)][,.(fit.lo = unlist(rqsslo))])
      
      suppressWarnings(strat.split.sim[[i]][, rqssmed := fitted(rqss(y ~ qss(x, lambda = exp(llam.qpred[[2]][[i]][[1]])),
                                                                     tau = qpred[2], na.action = na.exclude, .SD)), by = .(repl)][,.(fit.med = unlist(rqssmed))])
      
      suppressWarnings(strat.split.sim[[i]][, rqsshi := fitted(rqss(y ~ qss(x, lambda = exp(llam.qpred[[3]][[i]][[1]])),
                                                                    tau  = qpred[3], na.action = na.exclude, .SD)), by = .(repl)][,.(fit.hi = unlist(rqsshi))])
    }
  }
  
  sim <- rbindlist(strat.split.sim)
  setnames(sim, c("rqsslo", "rqssmed", "rqsshi"), qnames)
  
  sim.lb <- sim[, lapply(.SD, quantile, probs = qconf[[1]]), by = x.strat, .SDcols = qnames]
  sim.lb <- setnames(melt(sim.lb, id.vars = x.strat, measure.vars = qnames), c("variable", "value"), c("qname", "fit.lb"))
  
  sim.ub <- sim[, lapply(.SD, quantile, probs = qconf[[3]]), by = x.strat, .SDcols = qnames]
  sim.ub <- setnames(melt(sim.ub, id.vars = x.strat, measure.vars = qnames), c("variable", "value"), c("qname", "fit.ub"))
  
  sim <- sim[, lapply(.SD, median, na.rm = TRUE), by = x.strat, .SDcols = qnames]
  sim <- setnames(melt(sim, id.vars = x.strat, measure.vars = qnames), c("variable", "value"), c("qname", "fit"))
  
  sim <- sim[sim.lb, on = c(x.strat, "qname")]
  sim <- sim[sim.ub, on = c(x.strat, "qname")]
}

# Internal function for optimizing loess fit
.aicc.loess <- function(fit){
  #compute AIC_C for a LOESS fit, from:
  #Hurvich, C.M., Simonoff, J.S., and Tsai, C. L. 1998. Smoothing
  #parameter selection in nonparametric regression using an improved
  #Akaike Information Criterion. Journal of the Royal Statistical
  #Society B 60: 271–293.
  
  stopifnot(inherits(fit, 'loess'))
  n <- fit$n
  trace <- fit$trace.hat
  sigma2 <- sum(resid(fit) ^ 2) / (n - 1)
  return(log(sigma2) + 1 + (2 * (trace + 1)) / (n - trace - 2))
}

# Internal function for optimizing loess fit
.autoloess <- function(fit, span=c(.1, .9), ...){
  #compute loess fit which has span minimizes AIC_C
  #fit = loess fit; span parameter value doesn't matter
  #span = a two-value vector representing the minimum and maximum span values
  #Returns LOESS fit with span minimizing the AIC_C function
  
  stopifnot(inherits(fit, 'loess')) #, length(span) == 2)
  
  #loss function in form to be used by optimize
  f <- function(span) .aicc.loess(update(fit, span=span))
  
  #find best loess according to loss function
  opt.fit  <- update(fit, span=optimize(f, span)$minimum)
  opt.span <- optimize(f, span)$minimum
  return(list(fit=opt.fit, span=opt.span))
  
}

# Internal function for returning optimized span by strata
.getspan <- function(x) {
  span <- vector("list", length(x))
  for (i in seq_along(x)) {
    span[[i]] <- x[[i]]$span
  }
  names(span) <- names(x)
  return(span)
}

# Internal function for returning l.ypc fits by strata
.getfitted <- function(x) {
  fits <- vector("list", length(x))
  for (i in seq_along(x)) {
    fits[[i]] <- as.data.table(x[[i]]$loess[[1]]$fit$fitted)
    fits[[i]] <- setnames(fits[[i]], "fitted")
  }
  names(fits) <- names(x)
  return(fits)
}


.getfit <- function(x) {
  fit <- vector("list", length(x))
  for (i in seq_along(x)) {
    fit[[i]] <- x[[i]]$loess[[1]]$fit
  }
  names(fit) <- names(x)
  return(fit)
}

