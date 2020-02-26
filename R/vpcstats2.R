#' Perform a Visual Predictive Check (VPC) computation
#'
#' These functions work together to claculate the statistics that are plotted
#' in a VPC. They would typically be chained together using the "pipe" operator
#' (see Examples).
#'
#' @param o An object.
#' @param ... Additional arguments.
#' @examples
#'
#' \dontrun{
#' library(data.table)
#' library(magrittr)
#' library(vpc)
#'
#' exampleobs <- as.data.table(vpc::simple_data$obs)[MDV == 0]
#' examplesim <- as.data.table(vpc::simple_data$sim)[MDV == 0]
#' exampleobs$PRED <- examplesim[REP == 1, PRED]
#' exampleobs$SEX <- rep(c("F", "M"), len=nrow(exampleobs))
#' 
#' vpc <- observed(exampleobs, x=TIME, y=DV) %>%
#'     simulated(examplesim, y=DV) %>%
#'     stratify(~ ISM) %>%
#'     binning(TIME) %>%
#'     predcorrect(pred=PRED) %>%
#'     vpcstats()
#'
#' plot(vpc)
#'
#' # Example with 2-way stratification
#'
#' vpc <- vpc %>%
#'     stratify(SEX ~ ISM) %>%
#'     binning(TIME) %>%
#'     predcorrect(pred=PRED) %>%
#'     vpcstats()
#'
#' plot(vpc)
#' }
#'
#' @import data.table
#' @import magrittr
#' @importFrom stats median model.frame quantile setNames update
#' @name generics
NULL

#' @rdname generics
#' @export
observed <- function(o, ...) UseMethod("observed")

#' @rdname generics
#' @export
simulated <- function(o, ...) UseMethod("simulated")

#' @rdname generics
#' @export
censoring <- function(o, ...) UseMethod("censoring")

#' @rdname generics
#' @export
stratify <- function(o, ...) UseMethod("stratify")

#' @rdname generics
#' @export
binning <- function(o, ...) UseMethod("binning")

#' @rdname generics
#' @export
predcorrect <- function(o, ...) UseMethod("predcorrect")

#' @rdname generics
#' @export
nopredcorrect <- function(o, ...) UseMethod("nopredcorrect")

#' @rdname generics
#' @export
vpcstats <- function(o, ...) UseMethod("vpcstats")

#' @export
update.vpcstatsobj <- function(object, ...) {
    args <- list(...)
    for (i in names(args)) {
        object[[i]] <- args[[i, exact=TRUE]]
    }
    object
}

#' @export
observed.data.frame <- function(o, x, yobs, pred=NULL, blq, lloq=-Inf, uloq=Inf, ...) {
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
        blq  <- rlang::eval_tidy(rlang::enquo(blq),  data)
    }
    blq  <- as.logical(blq)

    if (missing(alq)) {
        alq <- (yobs > uloq)
    } else {
        alq  <- rlang::eval_tidy(rlang::enquo(alq),  data)
    }
    alq  <- as.logical(alq)

    obs <- data.table(x, y=yobs, blq, lloq, alq, uloq)

    o <- structure(list(data=data), class="vpcstatsobj")
    update(o, obs=obs, pred=pred)
}

#' @export
simulated.vpcstatsobj <- function(o, data, ysim, ...) {
    ysim <- rlang::eval_tidy(rlang::enquo(ysim), data)

    obs  <- o$obs
    x    <- obs$x
    nrep <- length(ysim)/nrow(obs)
    repl <- rep(1:nrep, each=nrow(obs))

    sim <- data.table(x, y=ysim, repl)
    update(o, sim=sim)
}

#' @export
censoring.vpcstatsobj <- function(o, blq, lloq, alq, uloq, data=o$data, ...) {
    if (missing(blq)) {
        blq <- o$blq
    } else {
        blq <- rlang::eval_tidy(rlang::enquo(blq), data)
    }
    if (missing(lloq)) {
        lloq <- o$lloq
    } else {
        lloq <- rlang::eval_tidy(rlang::enquo(lloq), data)
    }
    if (missing(alq)) {
        alq <- o$alq
    } else {
        alq <- rlang::eval_tidy(rlang::enquo(alq), data)
    }
    if (missing(uloq)) {
        uloq <- o$uloq
    } else {
        uloq <- rlang::eval_tidy(rlang::enquo(uloq), data)
    }

    if (is.null(blq)) {
        stop("No blq specified")
    }
    if (is.null(alq)) {
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

#' @export
stratify.vpcstatsobj <- function(o, formula, data=o$data, ...) {
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
    update(o, strat=strat, strat.formula=formula)
}

#' @export
binning.vpcstatsobj <- function(o, bin, data=o$data, ..., xbin="xmedian", centers, breaks, nbins, altx, stratum=NULL, by.strata=T) {

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

#' @export
predcorrect.vpcstatsobj <- function(o, pred, data=o$data, ..., log=FALSE) {

    ypc <- y <- NULL

    if (missing(pred)) {
        pred <- o$pred
    } else {
        pred <- rlang::eval_tidy(rlang::enquo(pred), data)
    }
    if (is.null(pred)) {
        stop("No pred specified")
    }

    stratbin <- o$.stratbin
    if (is.null(stratbin)) {
        stop("Need to specify binning before pred correction")
    }

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

    update(o, predcor=TRUE, pred=pred)
}

#' @export
nopredcorrect.vpcstatsobj <- function(o, ...) {
    update(o, predcor=FALSE)
}

#' @export
vpcstats.vpcstatsobj <- function(o, qpred=c(0.05, 0.5, 0.95), ..., conf.level=0.95, quantile.type=7) {

    repl <- ypc <- y <- blq <- lloq <- alq <- uloq <- NULL
    . <- list

    obs      <- o$obs
    sim      <- o$sim
    predcor  <- o$predcor
    stratbin <- o$.stratbin
    xbin     <- o$xbin

    if (is.null(stratbin)) {
        stop("Need to specify binning before calling vpcstats")
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

#' @describeIn bininfo Method for \code{vpcstatsobj}.
#' @param by.strata Should the calculations be done by strata? Defaults to what
#' was specified when the binning was done.
#' @export
bininfo.vpcstatsobj <- function(o, by.strata=o$bin.by.strata, ...) {

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

#' Print a \code{vpcstatsobj}.
#' @param x An object.
#' @param ... Further arguments can be specified but are ignored.
#' @return Returns \code{x} invisibly.
#' @export
print.vpcstatsobj <- function(x, ...) {
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

#' Plot a \code{vpcstatsobj}.
#' @param x An object.
#' @param show.points Should the observed data points be plotted?
#' @param show.boundaries Should the bin boundary be displayed?
#' @param show.stats Should the VPC stats be displayed?
#' @param show.binning Should the binning be displayed by coloring the observed data points by bin?
#' @param xlab A character label for the x-axis.
#' @param ylab A character label for the y-axis.
#' @param color A character vector of colors for the percentiles, from low to high.
#' @param linetype A character vector of linetyps for the percentiles, from low to high.
#' @param legend.position A character string specifying the position of the legend.
#' @param facet.scales A character string specifying the `scales` argument to use for facetting.
#' @param ... Further arguments can be specified but are ignored.
#' @return A `ggplot` object.
#' @seealso
#' \code{ggplot}
#' @export
plot.vpcstatsobj <- function(x, ..., show.points=TRUE, show.boundaries=TRUE, show.stats=!is.null(x$stats), show.binning=isFALSE(show.stats), xlab=NULL, ylab=NULL, color=c("red", "blue", "red"), linetype=c("dotted", "solid", "dashed"), legend.position="top", facet.scales="free") {

    xbin <- lo <- hi <- qname <- md <- y <- xleft <- xright <- ypc <- NULL
    . <- list

    vpc <- x

    qlvls <- levels(vpc$stats$qname)
    qlbls <- paste0(100*as.numeric(sub("^q", "", qlvls)), "%")

    if (isTRUE(vpc$predcor)) {
        ylab <- paste0(ylab, "\nPrediction Corrected")
    }

    has_ggplot2 <- requireNamespace("ggplot2", quietly=TRUE)
    if (!has_ggplot2) {
        stop("Package 'ggplot2' is required for plotting. Please install it to use this method.")
    }
    if (show.stats) {
        g <- ggplot2::ggplot(vpc$stats, ggplot2::aes(x=xbin)) +
            ggplot2::geom_ribbon(ggplot2::aes(ymin=lo, ymax=hi, fill=qname, col=qname, group=qname), alpha=0.1, col=NA) +
            ggplot2::geom_line(ggplot2::aes(y=md, col=qname, group=qname)) +
            ggplot2::geom_line(ggplot2::aes(y=y, linetype=qname), size=1) +
            ggplot2::scale_colour_manual(
                name=sprintf("Simulated Percentiles\nMedian (lines) %s%% CI (areas)", 100*vpc$conf.level),
                values=color,
                breaks=qlvls,
                labels=qlbls) +
            ggplot2::scale_fill_manual(
                name=sprintf("Simulated Percentiles\nMedian (lines) %s%% CI (areas)", 100*vpc$conf.level),
                values=color,
                breaks=qlvls,
                labels=qlbls) +
            ggplot2::scale_linetype_manual(
                name="Observed Percentiles\n(black lines)",
                values=linetype,
                breaks=qlvls,
                labels=qlbls) +
            ggplot2::guides(
                fill=ggplot2::guide_legend(order=2),
                colour=ggplot2::guide_legend(order=2),
                linetype=ggplot2::guide_legend(order=1))
    } else {
        g <- ggplot2::ggplot(vpc$strat)
    }

    g <- g + ggplot2::theme_bw() +
        ggplot2::theme(
            legend.key.width=ggplot2::unit(2, "lines"),
            legend.position=legend.position) +
        ggplot2::labs(x=xlab, y=ylab)

    if (show.points) {
        points.dat <- copy(vpc$obs)
        if (isTRUE(vpc$predcor)) {
            points.dat[, y := ypc]
        }
        if (show.binning) {
            reorder2 <- function(y, x) {
                y <- stats::reorder(y, x)
                (1:nlevels(y))[y]
            }
            points.dat[, color := reorder2(factor(bin), x), by=vpc$strat]
            points.dat[, color := factor(color)]
            g <- g + ggplot2::geom_point(data=points.dat, ggplot2::aes(x=x, y=y, color=color), size=1, alpha=0.4, show.legend=F) +
                scale_color_brewer(palette="Set1")
        } else {
            g <- g + ggplot2::geom_point(data=points.dat, ggplot2::aes(x=x, y=y), size=1, alpha=0.4)
        }
    }

    if (show.boundaries) {
        if (!is.null(vpc$strat)) {
            boundaries <- bininfo(vpc)[, .(x=sort(unique(c(xleft, xright)))), by=names(vpc$strat)]
        } else {
            boundaries <- bininfo(vpc)[, .(x=sort(unique(c(xleft, xright))))]
        }
        if (show.binning) {
            g <- g + ggplot2::geom_vline(data=boundaries, ggplot2::aes(xintercept=x), size=rel(0.5), col="gray80") + 
                ggplot2::theme(panel.grid=ggplot2::element_blank())
        }
        g <- g + ggplot2::geom_rug(data=boundaries, ggplot2::aes(x=x), sides="t", size=1)
    }

    if (!is.null(vpc$strat)) {
        if (length(as.list(vpc$strat.formula)) == 3) {
            g <- g + ggplot2::facet_grid(vpc$strat.formula, scales=facet.scales)
        } else {
            g <- g + ggplot2::facet_wrap(names(vpc$strat), scales=facet.scales)
        }
    }

    g
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

