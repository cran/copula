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


##' do the margin functions "p<nam>", "d<nam>" exist?
mvd.has.marF <- function(margins, prefix = "p")
    vapply(margins, function(M)
	   existsFunction(paste0(prefix, M)), NA)

mvdCheckM <- function(margins, prefix = "p") {
    ex <- mvd.has.marF(margins, prefix)
    if(any(!ex))
	warning("margins correct? Currently, have no function(s) named: ",
		paste(vapply(unique(margins[!ex]), function(M)
			     paste0(prefix, M), ""), collapse=", "))
}

mvdc <- function(copula, margins, paramMargins, marginsIdentical = FALSE,
		 check = TRUE, fixupNames = TRUE, warnMore = TRUE)
{
    if(check) {
	mvdCheckM(margins, "p")
	mvdCheckM(margins, "d")
    }
    ## marginsIdentical check at end ==> can be fast with just length-1 margins here
    if((fixupNames || check) && all(mvd.has.marF(margins, "p"))) {
	for(i in seq_along(margins)) {
	    n.i <- names(p.i <- paramMargins[[i]])
            ## names of formal args:
            nnms <- names(formals(get(cdfN <- paste0("p",margins[[i]])))[-1])
            ## but not the typical "non-parameter" arguments:
            nnms <- nnms[is.na(match(nnms, c("lower.tail", "log.p")))]
            if(!is.list(p.i) && !is.null(n.i))
                 warning(gettextf("paramMargins[[%d]] is not a list() but named", i), domain = NA)
	    if(fixupNames && (is.null(n.i) || !all(nzchar(n.i)))) {
		if(length(nnms) > (lpi <- length(p.i))) {
                    if(warnMore && (length(nnms) - lpi > 1 || !(cdfN == "pgamma" || "ncp" %in% nnms)))
                        warning(gettextf("%s() has more parameter arguments than length(paramMargins[[%d]]) == %d",
                                         cdfN, i, lpi), domain = NA)
                    length(nnms) <- lpi
                }
		if(length(nnms) > 0 &&
		   (is.null(n.i) || length(nnms) == length(n.i))) # careful ..
		   names(paramMargins[[i]]) <- nnms
	    }
	}
    }
    if (marginsIdentical) {
        dim <- dim(copula)
	if(length(margins) == 1)
	    margins <- rep(margins, dim)
        else if(length(margins) != dim) stop("marginsIdentical requires 1 or <d> identical margins")
	if(length(paramMargins) == 1)
	    paramMargins <- rep(paramMargins, dim)
        else if(length(paramMargins) != dim) stop("marginsIdentical requires 1 or <d> identical paramMargins")
    }
    ## validMvdc() called here:
    new("mvdc", copula = copula, margins = margins, paramMargins = paramMargins,
	marginsIdentical = marginsIdentical)
}

## "dim": via "xcopula" method

##' @title Parameter names of the margins of an "mvdc" object
##' @param mv
##' @return character vector of "the correct" length
##' @author Martin Maechler
margpnames <- function(mv) {
    nMar <- lengths(mv@paramMargins) # or vapply(mv@paramMargins, nFree, 1L)
    if(sum(nMar) == 0L) return(character())
    ## else
    p <- dim(mv@copula)
    pnms <- unlist(lapply(mv@paramMargins, names)) # maybe NULL
    if(mv@marginsIdentical) ## all the same ==> names only of *first* margin
	paste(paste("m", pnms[seq_len(nMar[1])], sep="."))
    else
	paste(paste0("m", rep.int(1:p, nMar)), pnms, sep=".")
}

## Function asCall was kindly supplied by
## Martin Maechler <maechler@stat.math.ethz.ch>,
## motivated by an application of nor1mix and copula
## from Lei Liu <liulei@virginia.edu>.
## They fixes the function getExpr in the old
## version, which assumed that the parameters to
## [rdpq]<distrib> were vectors.

asCall <- function(fun, param)
{
    cc <-
	if (length(param) == 0)
	    quote(FUN(x))
	else if(is.list(param)) {
	    as.call(c(quote(FUN), c(quote(x), as.expression(param))))
	} else { ## assume that [dpq]<distrib>(x, param) will work *correctly*
	    as.call(c(quote(FUN), c(quote(x), substitute(param))))
	}
    cc[[1]] <- as.name(fun)
    cc
}

dMvdc <- function(x, mvdc, log=FALSE) {
  dim <- dim(mvdc@copula)
  densmarg <- if(log) 0 else 1
  if (is.vector(x)) x <- matrix(x, nrow = 1)
  u <- x
  ## FIXME for large 'dim'  and  mvdc@marginsIdentical TRUE  can use *same* expr{ession}s for all i
  for (i in 1:dim) {
    cdf.expr <- asCall(paste0("p", mvdc@margins[i]), mvdc@paramMargins[[i]])
    dFnam <- paste0("d", mvdc@margins[[i]])
    pdfE <- asCall(dFnam, mvdc@paramMargins[[i]])
    pdf.expr <-
      if(log) { ## use d<distr>(x, <par>, log = TRUE) if possible  <==>  if(has.log)
        dFargNms <- names(formals(get(dFnam))[-1])
        has.log <- any("log" == dFargNms)
        pdf.expr <- if(has.log)
                        as.call(c(as.list(pdfE), list(log = TRUE)))
                    else substitute(log(PDF), list(PDF = pdfE))
      } else pdfE
    u[,i] <- eval(cdf.expr, list(x = x[,i]))
    densmarg <-
	if(log)
	    densmarg + eval(pdf.expr, list(x = x[,i]))
	else
	    densmarg * eval(pdf.expr, list(x = x[,i]))
  }
  if(log)
      dCopula(u, mvdc@copula, log=TRUE) + densmarg
  else
      dCopula(u, mvdc@copula) * densmarg
}

pMvdc <- function(x, mvdc) {
  dim <- dim(mvdc@copula)
  if (is.vector(x)) x <- matrix(x, nrow = 1)
  u <- x
  for (i in 1:dim) {
    cdf.expr <- asCall(paste0("p", mvdc@margins[i]), mvdc@paramMargins[[i]])
    u[,i] <- eval(cdf.expr, list(x = x[,i]))
  }
  pCopula(u, mvdc@copula)
}

rMvdc <- function(n, mvdc) {
  dim <- dim(mvdc@copula)
  u <- rCopula(n, mvdc@copula)
  x <- u
  for (i in 1:dim) {
    qdf.expr <- asCall(paste0("q", mvdc@margins[i]), mvdc@paramMargins[[i]])
    x[,i] <- eval(qdf.expr, list(x = u[,i]))
  }
  x
}

dmvdc <- function(mvdc, x, log=FALSE) { .Defunct("dMvdc"); dMvdc(x, mvdc, log) }
pmvdc <- function(mvdc, x) { .Defunct("pMvdc"); pMvdc(x, mvdc) }
rmvdc <- function(mvdc, n) { .Defunct("rMvdc"); rMvdc(n, mvdc) }

print.mvdc <- function(x, digits = getOption("digits"), ...)
{
    cat("Multivariate Distribution Copula based (\"mvdc\")\n @ copula:\n")
    print(x@copula, digits=digits, ...)
    cat(" @ margins:\n")
    print(x@margins, ...)
    margid <- x@marginsIdentical
    p <- dim(x)
    cat("   with", p, if(margid) "identical" else "(not identical)",
        " margins;")
    if(margid) {
        cat(" each with parameters\n")
        print(x@paramMargins[[1]], ...)
    } else {
        cat(" with parameters (@ paramMargins) \n")
        str(x@paramMargins, digits.d = digits)
    }
    invisible(x)
}

setMethod("show", signature("mvdc"), function(object) print.mvdc(object))
