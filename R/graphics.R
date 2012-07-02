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

##' @title Check if function 'fun' is like  pcopula or like pCopula
##' @param fun function such as pCopula, dcopula, ..
##' @return logical: TRUE if "like pCopula"
##' @author Martin Maechler
chkFun <- function(fun) {
    stopifnot(is.function(fun))
    isObj <- function(nm) any(nm == c("copula","mvdc"))
    nf <- names(formals(fun))
    if(isObj(nf[2])) TRUE
    else if(isObj(nf[1])) FALSE
    else NA # and the caller will get an error eventually
}



perspCopula <- function(x, fun, n = 51, theta = -30, phi = 30, expand = 0.618, ...) {
  eps <- (.Machine$double.eps)^(1/4)
  eps <- 0 ## FIXME - argument with default 0 ??
  xis <- yis <- seq(0 + eps, 1 - eps, len = n)
  grids <- as.matrix(expand.grid(xis, yis, KEEP.OUT.ATTRS=FALSE))
  zmat <- matrix(if(chkFun(fun)) fun(grids, x) else fun(x, grids), n, n)
  persp(xis, yis, zmat, theta = theta, phi = phi, expand = expand, ...)
  invisible(list(x = xis, y = yis, z = zmat))
}



contourCopula <- function(x, fun, n = 51,...) {
  eps <- (.Machine$double.eps)^(1/4)
  eps <- 0 ## FIXME - argument with default 0 ??
  xis <- yis <- seq(0 + eps, 1 - eps, len = n)
  grids <- as.matrix(expand.grid(xis, yis, KEEP.OUT.ATTRS=FALSE))
  zmat <- matrix(if(chkFun(fun)) fun(grids, x) else fun(x, grids), n, n)
  contour(xis, yis, zmat, ...)
  invisible(list(x = xis, y = yis, z = zmat))
}


perspMvdc <- function(x, fun,
                      xlim, ylim, nx = 51, ny = 51,
                      theta = -30, phi = 30, expand = 0.618, ...) {
  xis <- seq(xlim[1], xlim[2], length = nx)
  yis <- seq(ylim[1], ylim[2], length = ny)
  grids <- as.matrix(expand.grid(xis, yis, KEEP.OUT.ATTRS=FALSE))
  zmat <- matrix(if(chkFun(fun)) fun(grids, x) else fun(x, grids), nx, ny)
  persp(xis, yis, zmat, theta = theta, phi = phi, expand = expand, ...)
  invisible(list(x = xis, y = yis, z = zmat))
}


contourMvdc <- function(x, fun, xlim, ylim, nx = 51, ny = 51, ...)
{
  xis <- seq(xlim[1], xlim[2], length = nx)
  yis <- seq(ylim[1], ylim[2], length = ny)
  grids <- as.matrix(expand.grid(xis, yis, KEEP.OUT.ATTRS=FALSE))
  zmat <- matrix(if(chkFun(fun)) fun(grids, x) else fun(x, grids), nx, ny)
  contour(xis, yis, zmat, ...)
  invisible(list(x = xis, y = yis, z = zmat))
}

## special for independence copula:
setMethod("persp", signature("indepCopula"), perspCopula)
setMethod("contour", signature("indepCopula"), contourCopula)

## all other copulas:
setMethod("persp", signature("copula"), perspCopula)
setMethod("contour", signature("copula"), contourCopula)

setMethod("persp", signature("mvdc"), perspMvdc)
setMethod("contour", signature("mvdc"), contourMvdc)


### Graphical tools for detecting dependence
### Genest and Farve (2007, Journal of Hydrologic Engineering)

ChiPlot <- function(x, plot=TRUE, pval = 0.95, ...) {
#### originally proposed by Fisher and Switzer (1985 Biometrika, 2001 Am. Stat.)
  ## x is a n by 2 matrix
  if (!is.matrix(x)) x <- as.matrix(x)
  n <- nrow(x)
  hfg <- sapply(1:n,
                function(i) {
                  H <- (sum(x[,1] <= x[i,1] & x[,2] <= x[i,2]) - 1) / (n - 1)
                  F <- (sum(x[,1] <= x[i,1]) - 1) / (n - 1)
                  G <- (sum(x[,2] <= x[i,2]) - 1) / (n - 1)
                  c(H, F, G)
                })
  H <- hfg[1,]
  F <- hfg[2,]
  G <- hfg[3,]
  chi <-(H - F * G) / sqrt(F * (1 - F) * G * (1 - G))
  lambda <- 4 * sign( (F - 0.5) * (G - 0.5) ) * pmax( (F - 0.5)^2, (G - 0.5)^2 )
  cp <- c(1.54, 1.78, 2.18)
  idx <- pmatch(pval, c(0.9, 0.95, 0.99))
  if (is.na(idx)) stop("pval must be one of 0.9, 0.95, 0.99.")
  cp <- cp[idx]
  ymax <- max(abs(na.omit(chi)), cp / sqrt(n))
  if (plot) {
    plot(lambda, chi, xlim=c(-1, 1), ylim=c(-ymax, ymax), ...)
    abline(0, 0, lty = 3, col="gray")
    abline(cp / sqrt(n), 0, lty = 3, col="blue")
    abline(- cp / sqrt(n), 0, lty = 3, col="blue")
    lines(c(0, 0), c(-2, 2), lty = 3, col="gray")
  }
  invisible(cbind(H, F, G, chi, lambda))
}

KPlot <- function(x, plot=TRUE, ...) {
#### Genest and Boies (2003, American Statistician)
  if (!is.matrix(x)) x <- as.matrix(x)
  n <- nrow(x)
  H <- sapply(1:n,
              function(i) (sum(x[,1] <= x[i,1] & x[,2] <= x[i,2]) - 1) / (n - 1))
  H <- sort(H)
  K0 <- function(x) x - x * log(x)
  k0 <- function(x) - log(x)
  integrand <- function(w, i) w * k0(w) * K0(w)^(i - 1) * (1 - K0(w))^(n - i)
  W <- sapply(1:n,
              function(i) integrate(integrand, 0, 1, i = i,
                                    rel.tol=.Machine$double.eps^0.25)$value)
  W <- n * choose(n - 1, 1:n - 1) * W

  if (plot) {
    plot(W, H, xlim=c(0, 1), ylim=c(0, 1))
    curve(K0(x), add=TRUE, col="blue")
    abline(0, 1, col="gray")
  }
  invisible(cbind(H, W))
}

## x <- c(-2.224, -1.538, -0.807, 0.024, 0.052, 1.324)
## y <- c(0.431, 1.035, 0.586, 1.465, 1.115, -0.847)
## ChiPlot(cbind(x, y))
## KPlot(cbind(x, y))

################################################################################


##' @title A scatter plot matrix with nice variable names
##' @param data numeric matrix or as.matrix(.)able
##' @param varnames variable names, typically unspecified
##' @param Vname character string to become "root variable name"
##' @param col.mat matrix of colors
##' @param bg.col.mat matrix of background colors
##' @param ... further arguments to splom()
##' @return a splom() object
##' @author Martin Maechler
splom2 <- function(data, varnames=NULL, Vname="U", xlab="",
                   col.mat=NULL, bg.col.mat=NULL, ...)
{
    stopifnot(require(lattice),
	      is.numeric(data <- as.matrix(data)),
	      (d <- ncol(data)) >= 1)
    if(is.null(varnames)) {
	varnames <- do.call(expression,
			    lapply(1:d, function(i)
				   substitute(italic(A[I]), list(A = as.name(Vname), I=0+i))))
    }
    n <- nrow(data)
    if(is.null(col.mat))
        col.mat <- matrix(trellis.par.get("plot.symbol")$col, n,d)
    if(is.null(bg.col.mat))
        bg.col.mat <- matrix(trellis.par.get("background")$col, n,d)
    ## From Deepayan Sarkar, working around missing feature
    ##		(which should be in next release) of lattice
    my.diag.panel <- function(x, varname, ...)
        diag.panel.splom(x, varname=parse(text=varname), ...)
    ## splom
    splom(~data[,1:d], varnames=varnames, diag.panel=my.diag.panel, xlab="",
          panel = function(x, y, i, j, ...) {
              panel.fill(bg.col.mat[i,j])
              panel.splom(x, y, col=col.mat[i,j], ...)
          }, ...)
}


