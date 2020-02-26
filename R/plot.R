#' Plot a \code{tidyvpcobj}.
#' 
#' Use ggplot2 graphics to plot and customize the appearance of VPC
#' 
#' @title plot
#' @param x A tidyvpcobj object.
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
#' @param custom.theme A Character string specifying theme from ggplot2 package
#' @param ... Further arguments can be specified but are ignored.
#' @return A `ggplot` object.
#' @seealso
#' \code{ggplot}
#' @rdname plot
#' @export
plot.tidyvpcobj <- function(x, ..., show.points=TRUE, show.boundaries=TRUE, show.stats=!is.null(x$stats), show.binning=isFALSE(show.stats), xlab=NULL, ylab=NULL, color=c("red", "blue", "red"), linetype=c("dotted", "solid", "dashed"), legend.position="top", facet.scales="free", custom.theme = "theme_bw") {
  
  xbin <- lo <- hi <- qname <- md <- y <- xleft <- xright <- ypc <- l.ypc <- bin <- NULL
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
    if (!is.null(vpc$rqss.obs.fits)) {
      g <- ggplot2::ggplot(vpc$stats, ggplot2::aes(x = x)) +
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
      g <- ggplot2::ggplot(vpc$stats, ggplot2::aes(x = xbin)) +
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
    }
  } else {
    g <- ggplot2::ggplot(vpc$strat)
  }
  
  
  g <- g + eval(parse(text = paste0(custom.theme, "()"))) +
    ggplot2::theme(
      legend.key.width=ggplot2::unit(2, "lines"),
      legend.position=legend.position) +
    ggplot2::labs(x=xlab, y=ylab)
  
  if (show.points) {
    points.dat <- copy(vpc$obs)
    if (isTRUE(vpc$predcor)) {
      if(isTRUE(vpc$loess.ypc)) {
        points.dat[, y := l.ypc]
      } else {
        points.dat[, y := ypc]
      }
    }
    if (show.binning) {
      reorder2 <- function(y, x) {
        y <- stats::reorder(y, x)
        (1:nlevels(y))[y]
      }
      points.dat[, color := reorder2(factor(bin), x), by=vpc$strat]
      points.dat[, color := factor(color)]
      g <- g + ggplot2::geom_point(data=points.dat, ggplot2::aes(x=x, y=y, color=color), size=1, alpha=0.4, show.legend=F) +
        ggplot2::scale_color_brewer(palette="Set1")
    } else {
      g <- g + ggplot2::geom_point(data=points.dat, ggplot2::aes(x=x, y=y), size=1, alpha=0.4)
    }
  }
  
  if (show.boundaries) {
    if(is.null(vpc$rqss.obs.fits)) {
      if (!is.null(vpc$strat)) {
        boundaries <- bininfo(vpc)[, .(x=sort(unique(c(xleft, xright)))), by=names(vpc$strat)]
      } else {
        boundaries <- bininfo(vpc)[, .(x=sort(unique(c(xleft, xright))))]
      }
      if (show.binning) {
        g <- g + ggplot2::geom_vline(data=boundaries, ggplot2::aes(xintercept=x), size=ggplot2::rel(0.5), col="gray80") + 
          ggplot2::theme(panel.grid=ggplot2::element_blank())
      }
      g <- g + ggplot2::geom_rug(data=boundaries, ggplot2::aes(x=x), sides="t", size=1)
    }
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