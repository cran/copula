## Copyright (C) 2012 Marius Hofert and Martin Maechler
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.


### Pairs plot for a cu.u object ###############################################

##' Pairs plot for a cu.u object
##'
##' @title Pairs plot for a cu.u object
##' @param gcu.u (n,d,d)-array of pairwise Rosenblatt-transformed u's
##'	   as returned by \code{pairwiseCcop()}
##' @param panel pairs() argument
##' @param bgColList list of pairwise background colors and information as returned
##'	   by pairsColList()
##' @param labels pairs() argument; can be missing (in which case a suitable
##'	   default is chosen or can be "none" [or something else])
##' @param ... pairs() argument (can contain font.main, cex.main, line.main,
##'	   adj.main for title adjustments)
##' @param text.panel pairs() argument
##' @param label.pos pairs() argument
##' @param cex.labels pairs() argument
##' @param font.labels pairs() argument
##' @param gap pairs() argument
##' @param axes logical indicating whether axes are drawn (if g1 is not a
##'	   function, they are omitted automatically)
##' @param panel.border logical indicating whether a border is drawn around the
##'	   pairs (or not; to mimic the behavior of image())
##' @param key logical indicating whether a color key is drawn
##' @param key.space white space in height of characters in inch to specify the
##'	   the distance of the key to the pairs plot
##' @param key.width key width in height of characters in inch
##' @param key.axis logical indicating whether an axis for the color key is
##'	   drawn
##' @param key.title key title
##' @param key.line key placement (horizontal distance from color key in lines)
##' @param main title
##' @param main.centered logical indicating if the title should be centered or not;
##'	   the default (FALSE) centers it according to the pairs plot, not the whole
##'	   plotting region
##' @param main.line title placement (vertical distance from pairs plot in lines)
##' @param sub sub-title
##' @param sub.centered logical indicating if the sub-title should be centered or not;
##'	   the default (FALSE) centers it according to the pairs plot, not the whole
##'	   plotting region
##' @param sub.line sub-title placement (vertical distance from pairs plot in lines)
##' @return invisible
##' @author Marius Hofert
##' Note: based on pairs.default() and filled.contour() from R-2.14.1
pairs2 <- function(gcu.u,
		   panel=points, bgColList,
		   labels, ..., text.panel=textPanel,
		   label.pos=0.5, cex.labels=NULL, font.labels=1, gap=0,
		   axes=TRUE, panel.border=TRUE,
		   key=TRUE, key.space=2.5, key.width=1.5, key.axis=TRUE,
		   key.rug.at = numeric(), key.title=NULL, key.line=5, main=NULL,
		   main.centered=FALSE, main.line=2, sub=NULL, sub.centered=FALSE,
		   sub.line=4)
{
    ## checking
    stopifnot(is.array(gcu.u), is.numeric(gcu.u), length(dc <- dim(gcu.u)) == 3,
	      dc == c((n <- dc[1]),(d <- dc[2]), d), d>=2)

    panel <- match.fun(panel)
    if(d * d > 500 - 3) stop("Currently, layout() only allows 500 different plot areas\n This limit was now reached. It might be extended in the future.")

    ## preliminaries ###########################################################

    ## for plotting text on the diagonal panels
    textPanel <-
	function(x = 0.5, y = 0.5, txt, cex, font)
	    text(x, y, txt, cex = cex, font = font)

    ## function for drawing the (local) axes
    localAxis <- function(side, x, y, xpd, bg, col=NULL, main, sub, oma, ...) {
      ## Explicitly ignore any color argument passed in as
      ## it was most likely meant for the data points and
      ## not for the axis.
	if(side %%2 == 1) Axis(x, side=side, xpd=NA, ...)
	else Axis(y, side=side, xpd=NA, ...)
    }

    ## determine the labels on the diagonal
    has.labs <- TRUE
    if(missing(labels)){ ## if no labels are giving, choose a reasonable default:
	## the jth column in the pairs plot corresponds to the plot of (u_j, C(u_i|u_j)) => choose  .|u_j
	labels <- as.expression( lapply(1:d, function(j)
	    substitute(phantom()%.%"|"* italic(u)[ind], list(ind=j))) )
    }
    else { # labels are given (possibly == "none")
	labels <- match.arg(labels, "none")
	stopifnot(length(labels) == d || !(has.labs <- !(labels=="none")))
    }

    xls <- vapply(1:d, function(j)range(gcu.u[,j, j]), numeric(2))# those in the diag.
    yls <- vapply(1:d, function(i)range(gcu.u[,i,-i]), numeric(2))# those not in diag

    ## get ... and setup plot options
    op <- par(no.readonly = TRUE) # the whole list of settable par's.
    on.exit(par(op)) # on exit of this function, restore previous settings
    dots <- list(...) # list of ... arguments
    nmdots <- names(dots)
    oma <- if("oma" %in% nmdots) dots$oma else c(if(!is.null(sub)) 6 else 4, # below
						 4, # left
						 if(!is.null(main)) 6 else 4, # top
						 if(key) 6 else 4) # right
    par. <- par(mar=rep.int(gap/2,4), oma=oma, las=1) # set mar, oma, las and return previous settings of these parameters
    dev.hold()
    ## on.exit(dev.flush(), add=TRUE)
    ## care about *order* of exit calls:
    on.exit({dev.flush(); par(op)})
    ## cat("op:\n"); str(op)
    ## cat("on.exit():\n"); ex <- sys.on.exit(); str(ex)

    ## determine layout
    pc <- d*d # plot counter
    layout.mat <- matrix(1:pc, ncol=d) # first d^2 many plots for the pairs plot
    if(key){
	pc <- pc+1
	layout.mat <- cbind(layout.mat, rep(0, d), rep(pc, d)) # d^2 + 1 plot for color key; bind column of white space and the plot region for the color key to the layout
    }
    nLineAdded <- 0 # for substraction later on (in heights.)
    if(!is.null(main)){ # we have a title
	## Idea: If main.centered, create a full-width bar at the top of the plotting
	##	 region. If !main.centered, create such a bar only on top of the pairs
	##	 plot. Note that the height + distance do not play a role since we
	##	 turn off the clipping later anyway.
	pc <- pc+1
	layout.mat <- rbind(c(rep(pc, d), rep(if(main.centered)
					      pc else 0, 2)), layout.mat) # d^2 + 2 plot
	nLineAdded <- nLineAdded+1
    }
    if(!is.null(sub)){ # as for main
	pc <- pc+1
	layout.mat <- rbind(layout.mat, c(rep(pc, d), rep(if(sub.centered)
							  pc else 0, 2))) # d^2 + 3 plot
	nLineAdded <- nLineAdded+1
    }

    ## set the layout
    unit <- par("csi")*2.54
    widths. <- c(rep(1,d), if(key) lcm(c(key.space, key.width)*unit))
    heights. <- rep(1, nrow(layout.mat)-nLineAdded)
    ## add space for main (size does not matter since we will turn off clipping anyway)
    if(!is.null(main)) heights. <- c(lcm(1*unit), heights.)
    if(!is.null(sub)) heights. <- c(heights., lcm(1*unit)) # similarly for sub
    layout(layout.mat, respect=TRUE, widths=widths., heights=heights.)
    ## for debugging purposes:
    ## layout.show(pc)
    ## stop("stopped after showing layout")

    ## Pairs plot ##############################################################

    ## go over all panels (from top to bottom, left to right)
    has.bgColList <- !missing(bgColList)
    for(j in 1:d) { # column
	x <- gcu.u[,j,j]
	for(i in 1:d) { # row

	    ## start a new plot
	    plot.new()
	    ## determine whether we are on the diagonal
	    mfg <- par("mfg") # for checking afterwards
	    ok <- TRUE
	    if(i == j) {
		if(has.labs) {
		    ## The following line is commented out since it changes par("usr")
		    ## on the diagonal. I don't see why this is useful. Indeed, it
		    ## messes up the drawing of the second axes as "global box" if
		    ## panel.border==TRUE.
		    ## par(usr=c(0, 1, 0, 1))

		    plot.window(xlim= 0:1, ylim= 0:1)
		    ## determine the label sizes
		    if(is.null(cex.labels)){ # default
			l.wid <- strwidth(labels, "user")
			cex.labels <- max(0.8, min(2, .9/max(l.wid)))
		    }

		    ## draw the labels
		    text.panel(0.5, label.pos, labels[i], cex=cex.labels,
			       font=font.labels)
		}
	    } else { # i != j
		if(any(is.na(xls[,j]), is.na(yls[,i])))
		    ok <- FALSE
		if(ok) {
		    plot.window(xlim=xls[,j], ylim=yls[,i])
		    y <- gcu.u[,i,j]
		}
		## draw the colored background (rectangle)
		bg <- if(has.bgColList) bgColList$colMat[i,j] else "transparent"
		ll <- par("usr")
		rect(ll[1], ll[3], ll[2], ll[4], col=bg, if(!panel.border) border=NA) # set background color; xleft, ybottom, xright, ytop

		if(ok)## actual (panel) plot
		    panel(x, y, ...)
		##  =====
	    }
	    if(any(par("mfg")!=mfg)) stop("the 'panel' function made a new plot")

	    ## Now draw the axes (if required); reason for doing it here instead of
	    ## (as pairs()) above: since we have colored backgrounds, the colors would
	    ## overwrite parts of the axes if panel.border=TRUE
	    if(axes && ok) {
		if(i==1 && !(j %% 2)) localAxis(1 + 2, x, y, ...)
		if(i==d &&   j %% 2)  localAxis(3 - 2, x, y, ...)
		if(j==1 && !(i %% 2)) localAxis(2,     x, y, ...)
		if(j==d &&   i %% 2)  localAxis(4,     x, y, ...)
	    }

	    ## determine whether we draw boxes around the panels (also leading to a
	    ## global box) or just one global box
	    if(panel.border){ # a border is drawn (easy case)
		box() # draw it by simply drawing boxes around each panel
	    } else { # no border is drawn around the panels, but we still want to have a "global" box
		ll <- par("usr") # note: this changes if par(usr=c(0, 1, 0, 1)) is used above
		## note: par("usr") = c(x1, x2, y1, y2) for lower left point (x1, y1)
		##	 and upper right point (x2, y2)

		## rough idea in the following: use axis() in order to make a box
		## nicely aligned with (possibly) other axes drawn above; also tried
		## to work with segments but line thickness etc. is not the same.

		## upper border for 1st row
		if(i==1) axis(x, side=3, lwd.ticks=0, labels=FALSE,
		   at=c(ll[1], ll[2]), ...)
		## left border for 1st column
		if(j==1) axis(x, side=2, lwd.ticks=0, labels=FALSE,
		   at=c(ll[3], ll[4]), ...)
		## right border for last column
		if(j==d) axis(x, side=4, lwd.ticks=0, labels=FALSE,
		   at=c(ll[3], ll[4]), ...)
		## lower border for last row
		if(i==d) axis(x, side=1, lwd.ticks=0, labels=FALSE,
		   at=c(ll[1], ll[2]), ...)
	    }
	}
    }

    ## color key ###############################################################

    ## Note: font.main, cex.main are also needed below (if key == TRUE)
    font.main <- if("font.main" %in% nmdots) dots$font.main else par("font.main")
    cex.main <- if("cex.main" %in% nmdots) dots$cex.main else par("cex.main")
    if(key){

	## pick out color info
	cols <- bgColList$bucketCols # colors of the colorkey
	levels <- bgColList$pvalueBuckets # levels of the color key
	## checks (also that levels > 0 due to log-scale (e)axis)
	stopifnot((nl <- length(levels)) >= 2, length(cols) == nl-1)
	if(min(levels) <= 0)
	    stop("pvalueBuckets has to be > 0; most likely you want to specify pmin0 (> 0)") # longer error message

	## go to the space reserved for plotting the key
	plot.new()
	plot.window(xlim=c(0,1), ylim=range(levels), log="y", xaxs="i", yaxs="i")

	## plot rectangles that build the color key
	rect(0, levels[-nl], 1, levels[-1L], col=cols)

	## plot color key axis
	if(key.axis) sfsmisc::eaxis(4)

	## draw "data" ticks at left side of key
	if(length(key.rug.at) > 0)
	    axis(2, at=key.rug.at, labels=FALSE, tcl = -0.25)

	## plot color key title
	if(!is.null(key.title)) {
	    par(las=0) # label axis parallel to the axis
	    mtext(key.title, side=4, line=key.line, outer=FALSE, cex=cex.main,
		  font=font.main)
	}
    }

    ## title ###################################################################

    adj <- if("adj" %in% nmdots) dots$adj else par("adj") # also for sub
    if(!is.null(main)){
	plot.new()
	plot.window(xlim=c(0,1), ylim=c(0,1))
	## turn clipping off so that the title is visible if broader than the plot region
	par(xpd=NA)
	if(!is.list(main)){
	    mtext(main, side=3, line=main.line, adj=adj, cex=cex.main, font=font.main) # title
	} else { # for multi-line titles we accept a list
	    stopifnot(length(main) == length(main.line))
	    for(i in seq_along(main))
		mtext(main[[i]], side=3, line=main.line[i], adj=adj,
		      cex=cex.main, font=font.main)
	}
	par(xpd=FALSE) # restore default
    }

    ## sub-title ###############################################################

    if(!is.null(sub)){
	plot.new()
	plot.window(xlim=c(0,1), ylim=c(0,1))
	## turn clipping off so that the sub-title is visible if broader than the plot region
	par(xpd=NA)
	mtext(sub, side=3, line=-sub.line, adj=adj, cex=0.75*cex.main,
	      font=0.75*font.main) # title
	par(xpd=FALSE) # restore default
    }

    ## return invisibly
    invisible()
} # end{pairs2}


### colors #####################################################################

##' Create a list containing information about colors for a given matrix of p-values
##'
##' @title Create a list containing information about colors for a given matrix of
##'	   p-values
##' @param P matrix of p-values
##' @param pdiv vector of strictly increasing p-values in (0,1) that determine the
##'	   "buckets" for the background colors of pairs2()
##' @param signif.P significance level (must be an element of pdiv)
##' @param pmin0 a numeric indicating the lower endpoint of pvalueBuckets if pmin=0.
##'	   If set to 0, the lowest pvalueBucket will also be 0 (as pmin). If you want
##'	   to use pairsColList() for pairs2(), pmin0 should be in (0, min(pdiv))
##' @param colsBelow vector of (dark) colors for the coloring below signif.P.
##'	   Needs to contain at least two different colors as long as there are at least
##'	   two buckets below signif.P
##' @param colsAbove vector of (bright) colors for the coloring above signif.P.
##'	   Needs to contain at least two different colors as long as there are at least
##'	   two buckets above signif.P
##' @param colDiag color for the diagonal of P
##' @param space see ?colorRampPalette
##' @param ... additional arguments passed to colorRampPalette()
##' @return a list containing
##'	    1) a matrix colMat (colors corresponding to P)
##'	    2) bucketCols (a vector containing the colors corresponding to
##'	       pvalueBuckets)
##'	    3) a vector pvalueBuckets (containing the endpoints of the p-value buckets
##'	       determining the colors)
##' @author Marius Hofert
pairsColList <- function(P, pdiv=c(1e-4, 1e-3, 1e-2, 0.05, 0.1, 0.5), signif.P=0.05,
			 pmin0=0,
			 colsBelow=c("blue", "red"), colsAbove=c("orange", "white"),
			 colDiag="transparent", space="Lab", ...)
{
    stopifnot(is.na(P) | (0 <= P & P <= 1), is.matrix(P), (d <- ncol(P))==nrow(P),
	      (lp <- length(pdiv)) >= 1, 0 < pdiv, pdiv < 1, diff(pdiv) > 0,
	      length(signif.P)==1, any(is.Pv <- (signif.P == pdiv)),
	      length(pmin0)==1, 0 <= pmin0, pmin0 <= 1,
	      length(colDiag)==1, (ncolsBelow <- length(colsBelow)) >= 1,
	      (ncolsAbove <- length(colsAbove)) >= 1)

    ## 1) Determine the p-value buckets
    rP <- range(P, na.rm=TRUE) # range of p-values
    pmin <- rP[1] # minimal attained p-value
    pmax <- rP[2] # maximal attained p-value
    pvalueBuckets <- pdiv # p-value buckets endpoints
    if(pmin < pvalueBuckets[1]){ # if pmin < min(pvalueBuckets), extend pvalueBuckets to include pmin
	## Note: pmin could be 0 (due to simulation-based p-values, for example).
	##	 In this case, if pmin.adjust = TRUE, we adjust pmin to
	##	 (min(P, na.rm=TRUE) + min(pdiv))/2. This is useful when using
	##	 pairsColList for pairs2() and plotting the color key in log-scale;
	##	 there we have to avoid 0
	pvalueBuckets <- if(pmin > 0){
	    c(pmin, pvalueBuckets)
	} else { # pmin=0; if we are working in log-scale, we can't allow pmin=0
	    if(pmin0 >= pdiv[1]) stop("pmin0 must be smaller than min(pdiv)")
	    c(pmin0, pvalueBuckets) # note: this contains indeed 0 iff pmin0=0
	}
	lp <- lp+1
    }
    if(pmax > pvalueBuckets[lp]){ # if pmax > max(pvalueBuckets), extend pvalueBuckets to include pmax
	pvalueBuckets <- c(pvalueBuckets, pmax)
	lp <- lp+1
    }
    ## => pvalueBuckets now is typically in (0,1] and contains p_min and p_max.
    ##	  Note that 0 is included iff pmin0=0.

    ## Deal with the degenerate case that all p-values (entries of P) are equal to
    ## a single given pdiv. In this case, use (pdiv, 1) as single p-value bucket.
    ## This should rarely happen but it could (if d=1 and p-values are simulation
    ## based [=> grid], for example).
    if(lp == 1){
	pvalueBuckets <- c(pdiv, 1)
	lp <- lp+1
    }

    ## From here on it is guaranteed that length(pvalueBuckets) >= 2, which means
    ## we can use this "span" (determined by the (at least) two values) as proper
    ## endpoints of a bucket.

    ## 2) Determine which of the values of pvalueBuckets is signif.P
    is.Pv <- (signif.P == pvalueBuckets)
    if(!any(is.Pv)) stop("signif.P is not in pvalueBuckets (= extended pdiv)")
    num.signif.P <- which(is.Pv) # => num.signif.P can be all of {1,..,lp}

    ## 3) Determine number of different colors
    ncolsBelow. <- num.signif.P - 1 # number of different colors below signif.P
    nb <- lp-1 # number of buckets; lp = length(pvalueBuckets)
    ncolsAbove. <- nb - ncolsBelow.  # number of different colors above signif.P
    ## => Can both be 0, but not simultaneously:
    ##	  Example: if pdiv <- 0.05 = signif.P and all p-values (entries of P)
    ##		   are > 0.05 => no colors below signif.P.

    ## 4) Determine colors
    ##	  Note: As for the default, the larger the index, the brighter the color
    ##		vector should be so that the plot makes small p-values visible by
    ##		dark colors.

    ## colsBelow.
    colsBelow. <- if(ncolsBelow.==0){ # there is no bucket below signif.P
	character()
    } else if(ncolsBelow.==1){ # there is precisely one bucket below signif.P
	## in this case, we only need one color of colsBelow; if there are more given,
	## let's pick out the first one (typically the darkest)
	colsBelow[1]
    } else { # there are at least two buckets below signif.P
	## The problem is that colorRampPalette() fails if there is only one color
	## => we have to check that
	if(ncolsBelow==1) stop("colorRampPalette() needs at least two colors specified by 'colsBelow' to interpolate between in this case")
	colorRampPalette(colsBelow, space=space, ...)(ncolsBelow.)
    }

    ## colsAbove.
    colsAbove. <- if(ncolsAbove.==0){ # there is no bucket above signif.P
	character()
    } else if(ncolsAbove.==1){ # there is precisely one bucket above signif.P
	## in this case, we only need one color of colsAbove; if there are more given,
	## let's pick out the first one (more consistent than picking out the last
	## one which is typically the brightest)
	colsAbove[1]
    } else { # there are at least two buckets above signif.P
	## The problem is that colorRampPalette() fails if there is only one color
	## => we have to check that
	if(ncolsAbove==1) stop("colorRampPalette() needs at least two colors specified by 'colsAbove' to interpolate between in this case")
	colorRampPalette(colsAbove, space=space, ...)(ncolsAbove.)
    }

    ## colors for all buckets
    bucketCols <- c(colsBelow., colsAbove.)

    ## 5) find colors according to p-values
    ind <- findInterval(P, pvalueBuckets, all.inside=TRUE)
    ## Note:
    ## - findInterval can only return 0 if there is an entry in P which is smaller
    ##	 than the smalles value of pvalueBuckets. Due to the way pvalueBuckets is
    ##	 created, this can only happen if pmin=0 and pmin0>0. But in this case,
    ##	 all.inside=TRUE forces findInterval to return 1 instead of 0 which is what
    ##	 we need for an index
    ## - ind should now all be in {1,..,nb}
    ## - if P[i,j] == pvalueBuckets[k] for some k, then the bucket with
    ##	 upper endpoint pvalueBuckets[k] is returned
    ## - still contains NA (for the diagonal elements (order preserved))
    cols <- rep(colDiag, d*d) # vector of colors (filled with colDiag)
    nNA <- !is.na(ind)
    cols[nNA] <- bucketCols[ind[nNA]]
    colMat <- matrix(cols, nrow=nrow(P), ncol=ncol(P))

    ## 6) return a list containing the matrix of colors, the p-value "buckets"
    ##	  and corresponding colors
    list(colMat=colMat, bucketCols=bucketCols, pvalueBuckets=pvalueBuckets)
}


### Pairwise Rosenblatt / copula QQ-plot #######################################

##' Compute pairs plot based on pairwise Rosenblatt-transformed data
##'
##' @title Compute pairs plot based on pairwise Rosenblatt-transformed data
##' @param cu.u (n,d,d)-array of pairwise Rosenblatt-transformed u's
##'	   as returned by \code{pairwiseCcop()}
##' @param pvalueMat (d,d)-matrix of p-values
##' @param method string indicating the plot method
##' @param g1 function [0,1]^n -> [0,1]^n (n = length(u)); "x" for plotting
##'	   in one panel
##' @param g2 function [0,1]^{n x 2} -> [0,1]^n; "y" for plotting in one panel
##' @param bgColList list of pairwise background colors and information as returned
##'	   by pairsColList()
##' @param main title
##' @param sub sub title with a smart default
##' @param ... additional arguments passed to pairs2()
##' @return invisible
##' @author Marius Hofert
pairsRosenblatt <- function(cu.u, pvalueMat=pviTest(pairwiseIndepTest(cu.u)),
			    method = c("scatter", "QQchisq", "QQgamma",
			    "PPchisq", "PPgamma", "none"),
			    g1, g2, bgColList=pairsColList(pvalueMat, pmin0=1e-5),
			    main=NULL, sub = gpviString(pvalueMat),
			    key.title="p-value", key.rug=TRUE, ...)
{
    stopifnot(is.array(cu.u), length(dc <- dim(cu.u)) == 3,
	      dc == c((n <- dc[1]),(d <- dc[2]), d),
	      is.matrix(pvalueMat), dim(pvalueMat) == c(d,d))
    if(!missing(method)) {		# method given
	method <- match.arg(method)
	if(!missing(g1)) warning("'g1' will be ignored")
	if(!missing(g2)) warning("'g2' will be ignored")
	switch(method,
	       "scatter" = {
		   g1 <- function(u) u
		   g2 <- function(u, v) v
	       },
	       "QQchisq" = {
		   g1 <- function(u) qchisq(ppoints(length(u)), df=2)
		   g2 <- function(u, v) sort(qnorm(u)^2 + qnorm(v)^2)
	       },
	       "QQgamma" = {
		   g1 <- function(u) qgamma(ppoints(length(u)), shape=2)
		   g2 <- function(u, v) sort(-log(u) -log(v))
	       },
	       "PPchisq" = {
		   g1 <- function(u) ppoints(length(u))
		   g2 <- function(u, v) pchisq(sort(qnorm(u)^2 + qnorm(v)^2), df=2)
	       },
	       "PPgamma" = {
		   g1 <- function(u) ppoints(length(u))
		   g2 <- function(u, v) pgamma(sort(-log(u) -log(v)), shape=2)
	       },
	       "none" = {
		   g1 <- g2 <- NA
	       },
	       stop("unsupported method: ", method))
    } else {				# no method given
	if(!missing(g1) && !missing(g1)) {
	    ## check if they are functions and "look" reasonable:
	    stopifnot(is.function(g1), is.function(g2), length(formals(g2)) >= 2)
	} else { ## no method, no g1, no g2
	    if(!missing(g1) || !missing(g2)) stop("specify both 'g1' and 'g2' or none")
	    method <- "scatter"
	    g1 <- function(u) u
	    g2 <- function(u, v) v
	}
    }

    ## cu.u : (n,d,d)-array cu with cu[,i,j] containing C(u[,i]|u[,j]) for i!=j
    ##	       and u[,j] for i=j
    ## gcu.u: (n,d,d)-array gcu with gcu[,i,j] containing g2(u[,j], C(u[,i]|u[,j])) for i!=j
    ##	       and g1(u[,j]) for i=j

    gcu.u <- array(NA_real_, dim=c(n,d,d))
    if(is.function(g1)) {
	for(j in 1:d) {			# column
	    cuj <- cu.u[,,j]
	    uj <- cuj[,j]
	    gcu.u[,j,j] <- g1(uj)
	    for(i in (1:d)[-j])		# row
		gcu.u[,i,j] <- g2(uj, cuj[,i])
	}
    }
    ## plot
    pairs2(gcu.u, bgColList=bgColList, main=main, sub=sub,
	   key.title=key.title, key.rug.at = if(key.rug) pvalueMat, ...)
}


