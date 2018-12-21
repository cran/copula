## Copyright (C) 2012 Marius Hofert, Ivan Kojadinovic, Martin Maechler, and Jun Yan
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


### 0 Auxiliary functions ######################################################

##' @title Checker for 'FUN'
##' @param FUN A function such as pCopula, dcopula _or_ F.n()
##' @return A logical being TRUE if "like pCopula" _or_ F.n()
##' @author Martin Maechler
chkFun <- function(FUN)
{
    if(!is.function(FUN)) stop("the 'FUN' argument is not even a function")
    isObj <- function(nm) any(nm == c("copula", "mvdc")) ## [pdq][Cc]opula
    nf <- names(formals(FUN))
    if(isObj(nf[2]) || all(nf[1:2] == c("x","X"))) ## || is F.n()
        TRUE
    else if(isObj(nf[1])) FALSE
    else NA # and the caller will produce an error eventually
}


### 1 Elementary graphical tools for detecting dependence ######################

## NOTE: chiPlot() and Kplot() are currently not exported

## See Genest and Farve (2007, Journal of Hydrologic Engineering)
## Originally proposed by Fisher and Switzer (1985 Biometrika, 2001 Am. Stat.)
chiPlot <- function(x, plot = TRUE, pval = 0.95, ...)
{
    if (!is.matrix(x)) x <- as.matrix(x) # (n, 2)-matrix
    if(ncol(x) != 2) stop("must be a matrix of two columns")
    n <- nrow(x)
    hfg <- t(sapply(1:n, function(i) {
        tF <- (x[,1L] <= x[i,1L])
        tG <- (x[,2L] <= x[i,2L])
        c(H = sum(tF & tG), F = sum(tF), G = sum(tG))
    }) - 1) / (n - 1)
    H <- hfg[,"H"]
    F <- hfg[,"F"]
    G <- hfg[,"G"]
    chi <-(H - F * G) / sqrt(F * (1 - F) * G * (1 - G))
    lambda <- 4 * sign( (F - 0.5) * (G - 0.5) ) * pmax( (F - 0.5)^2, (G - 0.5)^2 )
    cp <- c(1.54, 1.78, 2.18) # critical points
    stopifnot(is.numeric(pval), length(pval) == 1L)
    idx <- which(abs(pval - c(0.9, 0.95, 0.99)) < 1e-6)
    if(length(idx) != 1L)
        stop("pval must be one of 0.9, 0.95, 0.99.")
    cp <- cp[idx]
    ymax <- max(abs(na.omit(chi)), cp / sqrt(n))
    if(plot) {
        plot(lambda, chi, xlim=c(-1, 1), ylim=c(-ymax, ymax), ...)
        abline(0, 0, lty = 3, col="gray")
        abline(cp / sqrt(n), 0, lty = 3, col="blue")
        abline(- cp / sqrt(n), 0, lty = 3, col="blue")
        lines(c(0, 0), c(-2, 2), lty = 3, col="gray")
    }
    invisible(cbind(hfg, chi, lambda))
}

## Genest and Boies (2003, American Statistician)
Kplot <- function(x, plot = TRUE, xlim = c(0,1), ylim = c(0,1), ...) {
    if (!is.matrix(x)) x <- as.matrix(x)
    if(ncol(x) != 2) stop("must be a matrix of two columns")
    n <- nrow(x)
### FIXME: should define another 'offset' (or offset of length 2!) such that we could use
### ----   H <- .Fn(x, offset = c(-1,-1))
    H <- (sapply(1:n,
                 function(i) (sum(x[,1] <= x[i,1] & x[,2] <= x[i,2]))) - 1) / (n - 1)
    H <- sort(H)
    K0 <- function(x) x - x * log(x)
    k0 <- function(x) - log(x)
    integrand <- function(w, i) w * k0(w) * K0(w)^(i - 1) * (1 - K0(w))^(n - i)
    W <- sapply(1:n,
                function(i) integrate(integrand, 0, 1, i = i,
                                      rel.tol=.Machine$double.eps^0.25)$value)
    W <- n * choose(n - 1, 1:n - 1) * W
    if (plot) {
        plot(W, H, xlim=xlim, ylim=ylim, ...)
        curve(K0(x), add=TRUE, col="blue")
        abline(0, 1, col="gray")
    }
    invisible(cbind(H, W))
}

## ==> Example for help file (and test):
## --> ../tests/misc.R
##     ^^^^^^^^^^^^^^^

### 2 Base graphics ############################################################

### 2.1 Legacy persp() and contour() methods (still improved, though) ##########

##' @title Perspective Plot Method for Class "Copula"
##' @param x An object of class "Copula"
##' @param FUN A function like dCopula or pCopula
##' @param n.grid The (vector of) number(s) of grid points in each dimension
##' @param delta Distance from the boundary of [0,1]^2
##' @param xlab The x-axis label
##' @param ylab The y-axis label
##' @param zlab The z-axis label
##' @param theta See ?persp
##' @param phi See ?persp
##' @param expand See ?persp
##' @param ticktype See ?persp
##' @param ... Additional arguments passed to persp()
##' @return invisible()
##' @author Marius Hofert (based on Ivan's/Jun's code)
##' @note *ARGHHHHHHHH...* persp() doesn't allow for plotmath x-/y-labels!
perspCopula <- function(x, FUN, n.grid = 26, delta = 0,
                        xlab = "u1", ylab = "u2",
                        zlab = deparse(substitute(FUN))[1], zlim = NULL,
                        theta = -30, phi = 30, expand = 0.618, ticktype = "detail", ...)
{
    stopifnot(n.grid >= 2)
    if(length(n.grid) == 1) n.grid <- rep(n.grid, 2)
    stopifnot(length(n.grid) == 2, 0 <= delta, delta < 1/2)
    xName <- deparse(substitute(x))
    fName <- deparse(substitute(FUN))
    x. <- seq(0 + delta, 1 - delta, length.out = n.grid[1])
    y. <- seq(0 + delta, 1 - delta, length.out = n.grid[2])
    grid <- as.matrix(expand.grid(x = x., y = y., KEEP.OUT.ATTRS = FALSE))
    z.mat <- matrix(if(chkFun(FUN)) FUN(grid, x) else FUN(x, grid), ncol = n.grid[2])
    if(is.null(zlim)) {
	zlim <- range(z.mat, na.rm = TRUE)
	if(diff(zlim) == 0) # persp.default would stop() --> better error message here :
	    stop(gettextf(
		"The non-NA values of '%s' are all equal to %g.  Not allowed for persp()",
		paste0(fName,"(xy, ",xName,")"), zlim[1]),
		call. = FALSE, domain = NA)
    }
    T <- persp(x = x., y = y., z = z.mat, zlim = zlim,
               xlab = xlab, ylab = ylab, zlab = zlab,
               theta = theta, phi = phi, expand = expand, ticktype = ticktype, ...)
    invisible(list(x = x., y = y., z = z.mat, persp = T))
}

##' @title Perspective Plot Method for Class "mvdc"
##' @param x An object of class "mvdc"
##' @param FUN A function like dCopula or pCopula
##' @param xlim The x-axis limits
##' @param ylim The y-axis limits
##' @param n.grid The (vector of) number(s) of grid points in each dimension
##' @param xlab The x-axis label
##' @param ylab The y-axis label
##' @param zlab The z-axis label
##' @param theta See ?persp
##' @param phi See ?persp
##' @param expand See ?persp
##' @param ticktype See ?persp
##' @param ... Additional arguments passed to persp()
##' @return invisible()
##' @author Marius Hofert (based on Ivan's/Jun's code)
##' @note *ARGHHHHHHHH...* persp() doesn't allow for plotmath x-/y-labels!
perspMvdc <- function(x, FUN, xlim, ylim, n.grid = 26,
                      xlab = "x1", ylab = "x2", zlab = deparse(substitute(FUN))[1],
                      theta = -30, phi = 30, expand = 0.618, ticktype = "detail", ...)
{
    stopifnot(n.grid >= 2)
    if(length(n.grid) == 1) n.grid <- rep(n.grid, 2)
    stopifnot(length(n.grid) == 2)
    x. <- seq(xlim[1], xlim[2], length.out = n.grid[1])
    y. <- seq(ylim[1], ylim[2], length.out = n.grid[2])
    grid <- as.matrix(expand.grid(x = x., y = y., KEEP.OUT.ATTRS = FALSE))
    z.mat <- matrix(if(chkFun(FUN)) FUN(grid, x) else FUN(x, grid), ncol = n.grid[2])
    T <- persp(x = x., y = y., z = z.mat,
               xlab = xlab, ylab = ylab, zlab = zlab,
               theta = theta, phi = phi, expand = expand, ticktype = ticktype, ...)
    invisible(list(x = x., y = y., z = z.mat, persp = T))
}

##' @title Contour Plot Method for Class "Copula"
##' @param x An object of class "Copula"
##' @param FUN A function like dCopula or pCopula
##' @param n.grid The (vector of) number(s) of grid points in each dimension
##' @param delta Distance from the boundary of [0,1]^2
##' @param xlab The x-axis label
##' @param ylab The y-axis label
##' @param box01 A logical indicating whether a box is drawn to on the boundary of [0,1]^2
##' @param ... Additional arguments passed to contour()
##' @return invisible()
##' @author Marius Hofert (based on Ivan's/Jun's code)
contourCopula <- function(x, FUN, n.grid = 26, delta = 0,
                          xlab = quote(u[1]), ylab = quote(u[2]),
                          box01 = TRUE, ...)
{
    stopifnot(n.grid >= 2)
    if(length(n.grid) == 1) n.grid <- rep(n.grid, 2)
    stopifnot(length(n.grid) == 2, 0 <= delta, delta < 1/2)
    x. <- seq(0 + delta, 1 - delta, length.out = n.grid[1])
    y. <- seq(0 + delta, 1 - delta, length.out = n.grid[2])
    grid <- as.matrix(expand.grid(x = x., y = y., KEEP.OUT.ATTRS = FALSE))
    z.mat <- matrix(if(chkFun(FUN)) FUN(grid, x) else FUN(x, grid), ncol = n.grid[2])
    contour(x = x., y = y., z = z.mat, xlab = xlab, ylab = ylab, ...)
    if(box01) rect(0, 0, 1, 1, border = "gray40", lty = 2)
    invisible(list(x = x., y = y., z = z.mat))
}

##' @title Contour Plot Method for Class "mvdc"
##' @param x An object of class "mvdc"
##' @param FUN A function like dCopula or pCopula
##' @param xlim The x-axis limits
##' @param ylim The y-axis limits
##' @param n.grid The (vector of) number(s) of grid points in each dimension
##' @param xlab The x-axis label
##' @param ylab The y-axis label
##' @param box01 A logical indicating whether a box is drawn to on the boundary of [0,1]^2
##' @param ... Additional arguments passed to contour()
##' @return invisible()
##' @author Marius Hofert
contourMvdc <- function(x, FUN, xlim, ylim, n.grid = 26,
                        xlab = quote(x[1]), ylab = quote(x[2]),
                        box01 = FALSE, ...)
{
    stopifnot(n.grid >= 2)
    if(length(n.grid) == 1) n.grid <- rep(n.grid, 2)
    else stopifnot(length(n.grid) == 2)
    x. <- seq(xlim[1], xlim[2], length.out = n.grid[1])
    y. <- seq(ylim[1], ylim[2], length.out = n.grid[2])
    grid <- as.matrix(expand.grid(x = x., y = y., KEEP.OUT.ATTRS = FALSE))
    z.mat <- matrix(if(chkFun(FUN)) FUN(grid, x) else FUN(x, grid), ncol = n.grid[2])
    contour(x = x., y = y., z = z.mat, xlab = xlab, ylab = ylab, ...)
    if(box01) rect(0, 0, 1, 1, border = "gray40", lty = 2)
    invisible(list(x = x., y = y., z = z.mat))
}

## Define persp() and contour() methods for objects of class "Copula" and "mvdc"
## Note: - No generic method (via setGeneric()) needed here as already provided by R
##       - The objects 'x' have to be named the same in perspCopula() and perspMvdc()
setMethod("persp",   signature = signature("Copula"), definition = perspCopula)
setMethod("persp",   signature = signature("mvdc"),   definition = perspMvdc)
setMethod("contour", signature = signature("Copula"), definition = contourCopula)
setMethod("contour", signature = signature("mvdc"),   definition = contourMvdc)

## "F.n(), C.n()" -- once we have empirical
## setMethod("persp", signature("mvFn"),
##           function(x, ...) {
##               perspCopula(x, F.n, ...)
##               warning("persp(<mvFn>, ..)  implementation unfinished;
##  contact maintainer(\"copula\")") ## FIXME : use trans3d() to add data; *or*
##               ## ensure seeing vertical jumps, by using {x_i-delta, x_i+delta} or ..
##           })


### 2.2 Enhanced Q-Q plot with rugs and confidence intervals ###################

##' @title Q-Q plot with rugs and pointwise asymptotic confidence intervals
##' @param x data (n-vector)
##' @param qF theoretical quantile function
##' @param log character string indicating whether log-scale should be used; see
##'        ?plot.default
##' @param qqline.args argument list passed to qqline(); use NULL to omit Q-Q line
##' @param rug.args argument list passed to rug(); use NULL to omit rugs
##' @param alpha significance level
##' @param CI.args argument list passed to lines() for drawing confidence intervals;
##'        use NULL to omit confidence intervals
##' @param CI.mtext argument list for information about confidence intervals; use
##'        NULL to omit information about confidence intervals
##' @param main title (can be an expression; use "" for no title)
##' @param main.args argument list passed to mtext() for drawing title
##' @param xlab x axis label
##' @param ylab y axis label
##' @param file file name (with extension .pdf) or "" (no pdf)
##' @param width width parameter of pdf()
##' @param height height parameter of pdf()
##' @param ... additional arguments passed to plot()
##' @return Q-Q plot
##' @author Marius Hofert
##' Note: - used in Genest, Hofert, Neslehova (2013)
##'       - better than pointwise asymptotic CIs would be (non-parametric)
##'         bootstrapped ones
qqplot2 <- function(x, qF, log="", qqline.args=if(log=="x" || log=="y") list(untf=TRUE) else list(),
                    rug.args=list(tcl=-0.6*par("tcl")),
                    alpha=0.05, CI.args = list(col="gray40"),
                    CI.mtext=list(text=paste0("Pointwise asymptotic ", 100*(1-alpha),
                                  "% confidence intervals"), side=4,
                                  cex=0.6*par("cex.main"), adj=0, col="gray40"),
                    main = quote(bold(italic(F) ~~ "Q-Q plot")),
                    main.args=list(text=main, side=3, line=1.1,
                                   cex=par("cex.main"), font=par("font.main"),
                                   adj=par("adj"), xpd=NA),
                    xlab="Theoretical quantiles", ylab="Sample quantiles",
                    file="", width=6, height=6, ...)
{
    x. <- sort(x) # drops NA
    n <- length(x.)
    p <- ppoints(n)
    q <- qF(p)

    ## Plot points
    doPDF <- nchar(file)
    if(doPDF) pdf(file=file, width=width, height=height)
    plot(q, x., xlab=xlab, ylab=ylab, log=log, ...) # empirical vs. theoretical quantiles
    if(length(main) && !identical(main, ""))
        do.call(mtext, main.args, quote = any(vapply(main.args, is.language, NA)))

    ## Plot the line (overplots points, but that's good for the eye!)
    if(!is.null(qqline.args))
        if(nchar(log)==1 && (is.null(untf <- qqline.args$untf) || !untf))
            warning("for a Q-Q line in x-log- or y-log-scale, specify 'untf = TRUE' in qqline.args")
    ## Draw the line (through the first and third quartile; see ?qqline)
    ## Note: abline(..., untf=TRUE) displays a curve (proper line in log-space)
    ##       *unless* both axes are in log-scale
        else {
            do.call(qqline,
                    args = c(if(log=="xy")
                                 list(y = log10(x.),
                                      distribution = function(p) log10(qF(p)))
                             else
                                 list(y = x., distribution = qF),
                             qqline.args),
                    quote = any(vapply(qqline.args, is.language, NA)))
        }
    ## Rugs
    if(!is.null(rug.args)) {
        doQ <- any(vapply(rug.args, is.language, NA))
        do.call(rug, c(list(q,  side=1), rug.args), quote=doQ)
        do.call(rug, c(list(x., side=2), rug.args), quote=doQ)
    }
    ## Confidence intervals
    if(!is.null(CI.args)) {
        ## Pointwise approximate (CLT) two-sided 1-alpha confidence intervals
        ## With up/low = p +/- qnorm(1-a/2)*sqrt(p(1-p)/n) it follows that
        ##   IP(F^{-1}(p) in [F_n^{-1}(low), F_n^{-1}(up)]) (p ~> F(x))
        ## = IP(x in [F_n^{-1}(low), F_n^{-1}(up)])
        ## = IP(F_n(x) in [low, up]) ~= 1-a since (CLT)
        ##   sqrt(n)*(F_n(x)-p)/sqrt(p*(1-p)) ~= N(0,1) since
        ## X_i~F => F_n(x)=bar{Z}_n with Z_i=I_{X_i<=x}~B(1, p=F(x)) => mean=p, var=p(1-p)
        SE <- sqrt(p*(1-p)/n) # standard error
        qa2 <- qnorm(1-alpha/2)
        ## Lower
        low <- p - qa2 * SE
        ilow <- 0 < low & low < 1 # in (0,1)
        low.y <- quantile(x., probs=low[ilow]) # F_n^{-1}(0<low<1)
        low.x <- q[ilow] # corresponding theoretical quantile at each ppoint
        ## Upper
        up <- p + qa2 * SE
        iup <- 0 < up & up < 1 # in (0,1)
        up.y <- quantile(x., probs=up[iup]) # F_n^{-1}(0<up<1)
        up.x <- q[iup] # corresponding theoretical quantile at each ppoint
        ## Draw lines
        doQ <- any(vapply(CI.args, is.language, NA))
        do.call(lines, c(list(low.x, low.y), CI.args), quote=doQ)
        do.call(lines, c(list( up.x,  up.y), CI.args), quote=doQ)
        ## Info
        if(!is.null(CI.mtext))
            do.call(mtext, CI.mtext, quote=any(vapply(CI.mtext, is.language, NA)))
    }
    if(doPDF) dev.off()
    invisible()
}


### 2.3 plot() methods #########################################################

##' @title Scatter Plot Method for Class "Copula"
##' @param x An object of class "Copula" (bivariate)
##' @param n The sample size
##' @param ... Additional arguments passed to plot()
##' @return invisible()
##' @author Marius Hofert
plotCopula <- function(x, n, xlim = 0:1, ylim = 0:1,
                       xlab = quote(U[1]), ylab = quote(U[2]), main=NULL, ...)
{
    stopifnot(n >= 1)
    if(dim(x) != 2)
        stop("The copula needs to be bivariate.")
    U <- rCopula(n, copula = x)
    xnam <- deparse(substitute(x, parent.frame())) # parent: as called from plot()
    if(isTRUE(main) || (is.null(main) && grepl("[Cc]opula\\b", xnam)))
	main <- xnam
    plot(U, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, main=main, ...)
}

##' @title Scatter Plot Method for Class "mvdc"
##' @param x An object of class "mvdc" (bivariate)
##' @param n The sample size
##' @param ... Additional arguments passed to plot()
##' @return invisible()
##' @author Marius Hofert
plotMvdc <- function(x, n, xlim = NULL, ylim = NULL,
                     xlab = quote(X[1]), ylab = quote(X[2]), ...)
{
    stopifnot(n >= 1)
    if(dim(x) != 2)
        stop("The multivariate distribution needs to be bivariate.")
    X <- rMvdc(n, mvdc = x)
    if(is.null(xlim)) xlim <- range(X[,1], finite = TRUE)
    if(is.null(ylim)) ylim <- range(X[,2], finite = TRUE)
    plot(X, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, ...)
}

## Define plot() methods for objects of various classes
## Note: 'x' is a "copula" or "mvdc" object
setMethod("plot", signature(x = "Copula"), plotCopula)
setMethod("plot", signature(x = "mvdc"),   plotMvdc)


### 2.4 Enhanced pairs() #######################################################

##' @title A Pairs Plot with Nice Defaults
##' @param x A numeric matrix or as.matrix(.)able
##' @param labels The labels, typically unspecified
##' @param labels.null.lab A character string to determine 'labels'
##'        in case 'labels' is NULL and 'x' does not have all column names given
##' @param ... Further arguments passed to pairs()
##' @return invisible()
##' @author Marius Hofert
pairs2 <- function(x, labels = NULL, labels.null.lab = "U", ...)
{
   stopifnot(is.numeric(x <- as.matrix(x)), (d <- ncol(x)) >= 1)
   if(is.null(labels)) {
       stopifnot(length(labels.null.lab) == 1, is.character(labels.null.lab))
       colnms <- colnames(x)
       labels <-
           if(sum(nzchar(colnms)) != d)
	       as.expression(lapply(1:d, function(i)
		   substitute(v[I], list(v = as.name(labels.null.lab), I = 0+i))))
           else # 'x' has column names => parse them
               parse(text = colnms)
   }
   pairs(x, gap = 0, labels = labels, ...)
}


### 3 Lattice graphics #########################################################

### 3.1 contourplot() methods ##################################################

##---> ../man/contourplot2-methods.Rd
##     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

##' @title Contour Plot Method for Classes "matrix" and "data.frame"
##' @param x A numeric matrix or as.matrix(.)able
##' @param aspect The aspect ratio
##' @param xlim, ylim  The x- and y-axis limits
##' @param xlab, ylab  The x-/y-axis label
##' @param cuts The number of levels
##' @param labels A logical indicating whether the contours should be labeled.
##'        Use labels = list(cex = 0.7) to have nicer small labels
##' @param scales See ?contourplot
##' @param region A logical indicating whether regions should be colored
##' @param ... Further arguments passed to contourplot()
##' @param col.regions The colors used for the colored regions
##' @return A contourplot() object
##' @author Marius Hofert
contourplot2Matrix <- function(x, aspect = 1,
                               xlim = range(x[,1], finite = TRUE),
                               ylim = range(x[,2], finite = TRUE),
                               xlab = NULL, ylab = NULL,
                               cuts = 16, labels = !region, pretty = !labels,
                               scales = list(alternating = c(1,1), tck = c(1,0)),
                               region = TRUE, ...,
                               ## Note that '...' must be before 'col.regions', otherwise
                               ## col = "red" is interpreted as 'col.regions = "red"'
                               col.regions = gray(seq(0.4, 1, length.out = max(100,4*cuts))))
{
    ## Checking
    if(!is.matrix(x)) x <- as.matrix(x)
    if(ncol(x) != 3) stop("'x' should have three columns")

    ## Labels
    colnms <- colnames(x)
    if(is.null(xlab)) xlab <- if(is.null(colnms)) "" else parse(text = colnms[1])
    if(is.null(ylab)) ylab <- if(is.null(colnms)) "" else parse(text = colnms[2])

    ## Contour plot
    contourplot(x[,3] ~ x[,1] * x[,2], aspect = aspect,
                xlim = extendrange(xlim, f = 0.025), # exactly the right balance to show all labels but no whitespace
                ylim = extendrange(ylim, f = 0.025),
                xlab = xlab, ylab = ylab, labels = labels, region = region,
                col.regions = col.regions, cuts = cuts, scales = scales, ...)
}

##' @title Contourplot Plot Method for Class "Copula"
##' @param x An object of class "Copula"
##' @param FUN A function like dCopula or pCopula
##' @param n.grid The (vector of) number(s) of grid points in each dimension
##' @param xlim The x-axis limits
##' @param ylim The y-axis limits
##' @param xlab The x-axis label
##' @param ylab The y-axis label
##' @param ... Additional arguments passed to contourplot2Matrix()
##' @return A contourplot() object
##' @author Marius Hofert
contourplot2Copula <- function(x, FUN, n.grid = 26, delta = 0,
                               xlim = 0:1, ylim = 0:1,
                               xlab = quote(u[1]), ylab = quote(u[2]),
                               ...)
{
    stopifnot(dim(x) == 2, n.grid >= 2)
    if(length(n.grid) == 1) n.grid <- rep(n.grid, 2)
    stopifnot(length(n.grid) == 2, 0 <= delta, delta < 1/2)
    x. <- seq(xlim[1] + delta, xlim[2] - delta, length.out = n.grid[1])
    y. <- seq(ylim[1] + delta, ylim[2] - delta, length.out = n.grid[2])
    grid <- as.matrix(expand.grid(x = x., y = y., KEEP.OUT.ATTRS = FALSE))
    z <- if(chkFun(FUN)) FUN(grid, x) else FUN(x, grid)
    val <- cbind(grid, z = z)
    contourplot2Matrix(val, xlim = xlim, ylim = ylim,
                       xlab = xlab, ylab = ylab, ...)
}

##' @title Contourplot Plot Method for Class "mvdc"
##' @param x An object of class "mvdc"
##' @param FUN A function like dCopula or pCopula
##' @param n.grid The (vector of) number(s) of grid points in each dimension
##' @param xlim The x-axis limits
##' @param ylim The y-axis limits
##' @param xlab The x-axis label
##' @param ylab The y-axis label
##' @param ... Additional arguments passed to contourplot2Matrix()
##' @return A contourplot() object
##' @author Marius Hofert
contourplot2Mvdc <- function(x, FUN, n.grid = 26, xlim, ylim,
                             xlab = quote(x[1]), ylab = quote(x[2]),
                             ...)
{
    stopifnot(dim(x) == 2, n.grid >= 2)
    if(length(n.grid) == 1) n.grid <- rep(n.grid, 2)
    stopifnot(length(n.grid) == 2)
    x. <- seq(xlim[1], xlim[2], length.out = n.grid[1])
    y. <- seq(ylim[1], ylim[2], length.out = n.grid[2])
    grid <- as.matrix(expand.grid(x = x., y = y., KEEP.OUT.ATTRS = FALSE))
    z <- if(chkFun(FUN)) FUN(grid, x) else FUN(x, grid)
    val <- cbind(grid, z = z)
    contourplot2Matrix(val, xlim = xlim, ylim = ylim,
                       xlab = xlab, ylab = ylab, ...)
}

## Define contourplot2() methods for objects of various classes
## Note: 'x' is a "matrix", "data.frame", "Copula" or "mvdc" object
setGeneric("contourplot2", function(x, ...) standardGeneric("contourplot2"))
setMethod("contourplot2", signature(x = "matrix"),     contourplot2Matrix)
setMethod("contourplot2", signature(x = "data.frame"), contourplot2Matrix)
setMethod("contourplot2", signature(x = "Copula"),     contourplot2Copula)
setMethod("contourplot2", signature(x = "mvdc"),       contourplot2Mvdc)


### 3.2 wireframe() methods ####################################################

##' @title Wireframe Plot Method for Classes "matrix" (and "data.frame")
##' @param x A numeric matrix or as.matrix(.)able
##' @param xlim, ylim, zlim  The x-, y-, and z-axis limits
##' @param xlab, ylab, zlab  The x-, y-, and z-axis label
##' @param alpha.regions See ?wireframe
##' @param scales See ?wireframe
##' @param par.settings Additional arguments passed to 'par.settings' (some are set)
##' @param draw.4.pCoplines logical indicating whether the four boundary
##'        lines are displayed
##' @param ... Further arguments passed to wireframe()
##' @return A wireframe() object
##' @author Marius Hofert
##' @note - axis.line makes the outer box disappear
##'       - 'col = "black"' in scales is required to make the ticks visible again
##'       - 'clip' is set to off to avoid axis labels being clipped
wireframe2Matrix <- function(x,
                             xlim = range(x[,1], finite = TRUE),
                             ylim = range(x[,2], finite = TRUE),
                             zlim = range(x[,3], finite = TRUE),
                             xlab = NULL, ylab = NULL, zlab = NULL,
                             alpha.regions = 0.5,
                             scales = list(arrows = FALSE, col = "black"),
                             par.settings = standard.theme(color = FALSE),
			     draw.4.pCoplines = FALSE, ...)
{
    ## Checking
    if(!is.matrix(x)) x <- as.matrix(x)
    if(ncol(x) != 3) stop("'x' should have three columns")

    ## Labels
    colnms <- colnames(x)
    if(is.null(xlab)) xlab <- if(is.null(colnms)) "" else parse(text = colnms[1])
    else ## workaround buglet in lattice [e-mailed Deepayan, 2016-11-16]:
        if(is.language(xlab) && !is.expression(xlab)) xlab <- as.expression(xlab)
    if(is.null(ylab)) ylab <- if(is.null(colnms)) "" else parse(text = colnms[2])
    else if(is.language(ylab) && !is.expression(ylab)) ylab <- as.expression(ylab)
    if(is.null(zlab)) zlab <- list(
                          if(is.null(colnms)) "" else parse(text = colnms[3]),
                          rot = 90)
    else if(is.language(zlab) && !is.expression(zlab)) zlab <- as.expression(zlab)
    par.settings <- modifyList(par.settings,
                               list(axis.line = list(col = "transparent"),
                                    clip = list(panel = "off")))
    ## Return wireframe() object :
    if(is.function(list(...)$panel.3d.wireframe)) ## use it from ... => do *not* set it:
        wireframe(x[,3] ~ x[,1] * x[,2], xlim = xlim, ylim = ylim, zlim = zlim,
                  xlab = xlab, ylab = ylab, zlab = zlab,
                  alpha.regions = alpha.regions, scales = scales,
                  par.settings = par.settings, ...)
    else
        wireframe(x[,3] ~ x[,1] * x[,2], xlim = xlim, ylim = ylim, zlim = zlim,
                  xlab = xlab, ylab = ylab, zlab = zlab,
                  alpha.regions = alpha.regions, scales = scales,
                  panel.3d.wireframe = if(draw.4.pCoplines) panel.3dwire.4 else panel.3dwire,
                  par.settings = par.settings, ...)
}

##' wireframe panel additionally drawing "the" 4-line-segments of every pCopula()
panel.3dwire.4 <- function(x, y, z, rot.mat, distance,
                           xlim.scaled, ylim.scaled, zlim.scaled,
                           ## NB: These are *documented* in ../man/wireframe2-methods.Rd
                           col.4 = "#668b5580", lwd.4 = 5, lty.4 = "82", ...)
{
    panel.3dwire(x, y, z, rot.mat, distance,
                 xlim.scaled=xlim.scaled, ylim.scaled=ylim.scaled, zlim.scaled=zlim.scaled,
                 ...)
    X <- xlim.scaled; x0 <- X[1]; x1 <- X[2]
    Y <- ylim.scaled; y0 <- Y[1]; y1 <- Y[2]
    Z <- zlim.scaled; z0 <- Z[1]# z1 <- Z[2]
    segs <- list(rbind(X,  y0, z0), # X-axis
                 rbind(x0,  Y, z0), # Y-axis
                 rbind(x1,  Y,  Z), # Y-Z diagonal (at x1)
                 rbind( X, y1,  Z)) # X-Z diagonal (at y1)
    for(seg in segs) {
        m <- ltransform3dto3d(seg, rot.mat, distance)
        lsegments(x0= m[1,1], y0= m[2,1], col = col.4,
                  x1= m[1,2], y1= m[2,2], lwd = lwd.4, lty = lty.4)
    }
}

##' @title Wireframe Plot Method for Class "Copula" ---> ../man/wireframe2-methods.Rd
##' @param FUN A function like dCopula or pCopula
wireframe2Copula <- function(x, FUN, n.grid = 26, delta = 0,
                             xlim = 0:1, ylim = 0:1, zlim = NULL,
                             xlab = quote(u[1]), ylab = quote(u[2]),
                             zlab = list(deparse(substitute(FUN))[1], rot = 90),
                             draw.4.pCoplines = identical(FUN, pCopula), ...)
{
    stopifnot(dim(x) == 2, n.grid >= 2)
    if(length(n.grid) == 1) n.grid <- rep(n.grid, 2)
    stopifnot(length(n.grid) == 2, 0 <= delta, delta < 1/2)
    x. <- seq(xlim[1] + delta, xlim[2] - delta, length.out = n.grid[1])
    y. <- seq(ylim[1] + delta, ylim[2] - delta, length.out = n.grid[2])
    grid <- as.matrix(expand.grid(x = x., y = y., KEEP.OUT.ATTRS = FALSE))
    z <- if(chkFun(FUN)) FUN(grid, x) else FUN(x, grid)
    if(is.null(zlim)) zlim <- range(z, finite = TRUE)
    wireframe2Matrix(cbind(grid, z = z), xlim = xlim, ylim = ylim, zlim = zlim,
                     xlab = xlab, ylab = ylab, zlab = zlab,
                     draw.4.pCoplines = draw.4.pCoplines, ...)
}

##' @title Wireframe Plot Method for Class "mvdc"   --->  ../man/wireframe2-methods.Rd
##' @param FUN A function like dMvdc or pMvdc
wireframe2Mvdc <- function(x, FUN, n.grid = 26,
                           xlim, ylim, zlim = NULL,
                           xlab = quote(x[1]), ylab = quote(x[2]),
                           zlab = list(deparse(substitute(FUN))[1], rot = 90), ...)
{
    stopifnot(dim(x) == 2, n.grid >= 2)
    if(length(n.grid) == 1) n.grid <- rep(n.grid, 2)
    stopifnot(length(n.grid) == 2)
    x. <- seq(xlim[1], xlim[2], length.out = n.grid[1])
    y. <- seq(ylim[1], ylim[2], length.out = n.grid[2])
    grid <- as.matrix(expand.grid(x = x., y = y., KEEP.OUT.ATTRS = FALSE))
    z <- if(chkFun(FUN)) FUN(grid, x) else FUN(x, grid)
    if(is.null(zlim)) zlim <- range(z, finite = TRUE)
    val <- cbind(grid, z = z)
    wireframe2Matrix(val, xlim = xlim, ylim = ylim, zlim = zlim,
                     xlab = xlab, ylab = ylab, zlab = zlab, ...)
}

## Define wireframe2() methods for objects of various classes
## Note: 'x' is a "matrix", "data.frame", "Copula" or "mvdc" object
setGeneric("wireframe2", function(x, ...) standardGeneric("wireframe2"))
setMethod("wireframe2", signature(x = "matrix"),     wireframe2Matrix)
setMethod("wireframe2", signature(x = "data.frame"), wireframe2Matrix)
setMethod("wireframe2", signature(x = "Copula"),     wireframe2Copula)
setMethod("wireframe2", signature(x = "mvdc"),       wireframe2Mvdc)


### 3.3 cloud() methods ########################################################

## --> ../man/cloud2-methods.Rd
##            =================
##' @title Cloud Plot Method for Classes "matrix" and "data.frame"
##' @param x A numeric matrix or as.matrix(.)able
##' @param xlim The x-axis limits
##' @param ylim The y-axis limits
##' @param zlim The z-axis limits
##' @param xlab The x-axis label
##' @param ylab The y-axis label
##' @param zlab The z-axis label
##' @param scales See ?cloud
##' @param par.settings Additional arguments passed to 'par.settings' (some are set)
##' @param ... Further arguments passed to cloud()
##' @return A "trellis" object, as returned from cloud()
##' @author Marius Hofert
cloud2Matrix <- function(x,
                         xlim = range(x[,1], finite = TRUE),
                         ylim = range(x[,2], finite = TRUE),
                         zlim = range(x[,3], finite = TRUE),
                         xlab = NULL, ylab = NULL, zlab = NULL,
                         scales = list(arrows = FALSE, col = "black"),
                         par.settings = standard.theme(color = FALSE), ...)
{
    ## Checking
    if(!is.matrix(x)) x <- as.matrix(x)
    if(ncol(x) != 3) stop("'x' should have three columns")

    ## Labels
    colnms <- colnames(x)
    if(is.null(xlab)) xlab <- if(is.null(colnms)) "" else parse(text=colnms[1])
    else ## workaround buglet in lattice [e-mailed Deepayan, 2016-11-16]:
        if(is.language(xlab) && !is.expression(xlab)) xlab <- as.expression(xlab)
    if(is.null(ylab)) ylab <- if(is.null(colnms)) "" else parse(text=colnms[2])
    else if(is.language(ylab) && !is.expression(ylab)) ylab <- as.expression(ylab)
    if(is.null(zlab)) zlab <- if(is.null(colnms)) "" else parse(text=colnms[3])
    else if(is.language(zlab) && !is.expression(zlab)) zlab <- as.expression(zlab)

    ## Cloud plot
    cloud(x[,3] ~ x[,1] * x[,2],
          xlim = extendrange(xlim, f = 0.04),
          ylim = extendrange(ylim, f = 0.04),
          zlim = extendrange(zlim, f = 0.04),
          xlab = xlab, ylab = ylab, zlab = zlab, scales = scales,
          par.settings = modifyList(par.settings,
                                    list(axis.line = list(col = "transparent"),
                                         clip = list(panel = "off"))),
          ...)
}

##' @title Cloud Plot Method for Class "Copula"
##' @param x An object of class "Copula"
##' @param n The sample size
##' @param xlim The x-axis limits
##' @param ylim The y-axis limits
##' @param zlim The z-axis limits
##' @param xlab The x-axis label
##' @param ylab The y-axis label
##' @param zlab The z-axis label
##' @param ... Additional arguments passed to cloud2Matrix()
##' @return A cloud() object
##' @author Marius Hofert
cloud2Copula <- function(x, n,
                         xlim = 0:1, ylim = 0:1, zlim = 0:1,
			 xlab = quote(U[1]), ylab = quote(U[2]), zlab = quote(U[3]),
			 ...)
{
    stopifnot(n >= 1)
    if(dim(x) != 3)
        stop("The copula needs to be of dimension 3.")
    U <- rCopula(n, copula = x)
    cloud2Matrix(U, xlim = xlim, ylim = ylim, zlim = zlim,
                 xlab = xlab, ylab = ylab, zlab = zlab, ...)
}

##' @title Cloud Plot Method for Class "mvdc"
##' @param x An object of class "mvdc"
##' @param n The sample size
##' @param xlim The x-axis limits
##' @param ylim The y-axis limits
##' @param zlim The z-axis limits
##' @param xlab The x-axis label
##' @param ylab The y-axis label
##' @param zlab The z-axis label
##' @param ... Additional arguments passed to cloud2Matrix()
##' @return A cloud() object
##' @author Marius Hofert
cloud2Mvdc <- function(x, n,
                       xlim = NULL, ylim = NULL, zlim = NULL,
		       xlab = quote(X[1]), ylab = quote(X[2]), zlab = quote(X[3]),
		       ...)
{
    stopifnot(n >= 1)
    if(dim(x) != 3)
        stop("The multivariate distribution needs to be of dimension 3.")
    X <- rMvdc(n, mvdc = x)
    if(is.null(xlim)) xlim <- range(X[,1], finite = TRUE)
    if(is.null(ylim)) ylim <- range(X[,2], finite = TRUE)
    if(is.null(zlim)) zlim <- range(X[,3], finite = TRUE)
    cloud2Matrix(X, xlim = xlim, ylim = ylim, zlim = zlim,
                 xlab = xlab, ylab = ylab, zlab = zlab, ...)
}

## Define cloud2() methods for objects of various classes
## Note: 'x' is a "matrix", "data.frame", "Copula" or "mvdc" object
setGeneric("cloud2", function(x, ...) standardGeneric("cloud2"))
setMethod("cloud2", signature(x = "matrix"),     cloud2Matrix)
setMethod("cloud2", signature(x = "data.frame"), cloud2Matrix)
setMethod("cloud2", signature(x = "Copula"),     cloud2Copula)
setMethod("cloud2", signature(x = "mvdc"),       cloud2Mvdc)


### 3.4 splom() methods ########################################################

##' @title A Scatter-plot Matrix (SPLOM) with Nice Defaults
##' @param x A numeric matrix or as.matrix(.)able
##' @param varnames The variable names, typically unspecified
##' @param varnames.null.lab A character string to determine 'varnames'
##'        in case 'varnames' is NULL and 'x' does not have all column names given
##' @param xlab The x-axis label
##' @param col.mat A matrix of colors
##' @param bg.col.mat A matrix of background colors
##' @param ... Further arguments passed to splom()
##' @return An splom() object
##' @author Martin Maechler and Marius Hofert
splom2Matrix <- function(x, varnames = NULL,
                         varnames.null.lab = "U", xlab = "",
                         col.mat = NULL, bg.col.mat = NULL, ...)
{
    if(!is.matrix(x)) x <- as.matrix(x)
    stopifnot((d <- ncol(x)) >= 1)
    if(is.null(varnames)) {
        stopifnot(length(varnames.null.lab) == 1, is.character(varnames.null.lab))
        colnms <- colnames(x)
	varnames <-
	    if(sum(nzchar(colnms)) != d) {
		vNm <- as.name(varnames.null.lab)
		as.expression(lapply(1:d, function(j)
		    substitute(v[I], list(v = vNm, I = j))))
	    } else # 'x' has column names => use them
		parse(text = colnms)
    }
    n <- nrow(x)
    if(is.null(col.mat))
        col.mat <- matrix(trellis.par.get("plot.symbol")$col, n, d)
    else if(length(col.mat) == 1) # for full matrix
        col.mat <- matrix(col.mat, n, d)
    if(is.null(bg.col.mat))
        bg.col.mat <- matrix(trellis.par.get("background")$col, n, d)
    else if(length(bg.col.mat) == 1)
        bg.col.mat <- matrix(bg.col.mat, n, d)

    splom( ~ x, varnames = varnames, xlab = xlab,
          panel = function(x, y, i, j, ...) {
              panel.fill(bg.col.mat[i,j])
              panel.splom(x, y, col = col.mat[i,j], ...)
          }, ...)
}

##' @title Scatter-Plot Matrix Method for Class "Copula"
##' @param x An object of class "Copula"
##' @param n The sample size
##' @param ... Additional arguments passed to splom2Matrix()
##' @return An splom() object
##' @author Marius Hofert
splom2Copula <- function(x, n, ...)
{
    stopifnot(n >= 1)
    if(dim(x) <= 2)
        stop("The copula needs to be of dimension >= 3.")
    U <- rCopula(n, copula = x)
    splom2Matrix(U, ...)
}

##' @title Scatter-Plot Matrix Method for Class "mvdc"
##' @param x An object of class "mvdc"
##' @param n The sample size
##' @param ... Additional arguments passed to splom2Matrix()
##' @return An splom() object
##' @author Marius Hofert
splom2Mvdc <- function(x, n, varnames.null.lab = "X", ...)
{
    stopifnot(n >= 1)
    if(dim(x) <= 2)
        stop("The multivariate distribution needs to be of dimension >= 3.")
    X <- rMvdc(n, mvdc = x)
    splom2Matrix(X, varnames.null.lab = varnames.null.lab, ...)
}

## Define splom2() methods for objects of various classes
## Note: 'x' is a "matrix", "data.frame", "Copula" or "mvdc" object
setGeneric("splom2", function(x, ...) standardGeneric("splom2"))
setMethod("splom2", signature(x = "matrix"),     splom2Matrix)
setMethod("splom2", signature(x = "data.frame"), splom2Matrix)
setMethod("splom2", signature(x = "Copula"),     splom2Copula)
setMethod("splom2", signature(x = "mvdc"),       splom2Mvdc)
