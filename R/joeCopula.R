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


joeCopula <- function(param = NA_real_, dim = 2L,
		      use.indepC = c("message", "TRUE", "FALSE"))
{
  stopifnot(length(param) == 1)
  if(!is.na(param) && param == 1) {
    use.indepC <- match.arg(use.indepC)
    if(!identical(use.indepC, "FALSE")) {
      if(identical(use.indepC, "message"))
        message("parameter at boundary ==> returning indepCopula()")
      return( indepCopula(dim=dim) )
    }
  }
  ## cdf expression
  cdfExpr <- function(d) {
    term <- paste0("(1 - (1 - ", paste0("u", 1:d), ")^alpha)", collapse = " * ")
    cdf <- parse(text=paste0("1 - (1 - ", term, ")^(1 / alpha)"))
  }
  cdf <- cdfExpr(dim)
  pdf <- if (dim <= 6) cdfExpr2pdfExpr(cdf, dim)
  exprdist <- structure(c(cdf = cdf, pdf = pdf),
                        cdfalgr = if (!is.null(cdf)) deriv(cdf, "nothing"), # <-- FIXME: far from optimal
                        pdfalgr = if (!is.null(pdf)) deriv(pdf, "nothing")) # <-- FIXME: ditto

  new("joeCopula",
      dimension = as.integer(dim),
      parameters = param[1],
      exprdist = exprdist,
      param.names = "alpha",
      param.lowbnd = 1, # 0.238733989880086 for tau >= -1 -- is NOT valid
      param.upbnd = Inf,
      fullname = "<deprecated slot>")# "Joe copula family; Archimedean copula"
}

setMethod("rCopula", signature("numeric", "joeCopula"),
	  function (n, copula, ...)
	  racopula(n, copJoe, theta=copula@parameters, d = copula@dimension))

setMethod("pCopula", signature("matrix", "joeCopula"),
	  function(u, copula, ...) .pacopula(u, copJoe, theta=copula@parameters))

setMethod("dCopula" , signature("matrix", "joeCopula"),
	  function (u, copula, log = FALSE, ...)
	  copJoe@dacopula(u, theta=copula@parameters, log=log, ...))

## pCopula() and dCopula() *generic* already deal with non-matrix case!
## setMethod("pCopula", signature("numeric", "joeCopula"),
## 	  function(u, copula, ...) pacopula(u, copJoe, theta=copula@parameters))
## setMethod("dCopula", signature("numeric", "joeCopula"),
## 	  function (u, copula, log = FALSE, ...)
## 	  copJoe@dacopula(rbind(u, deparse.level=0L),
## 			  theta=copula@parameters, log=log, ...))

setMethod("psi", signature("joeCopula"),
	  function(copula, s) copJoe@psi(t=s, theta=copula@parameters))
setMethod("iPsi", signature("joeCopula"),
	  function(copula, u) copJoe@iPsi(u, theta=copula@parameters))
setMethod("diPsi", signature("joeCopula"),
	  function(copula, u, degree=1, log=FALSE, ...) {
    s <- if(log || degree %% 2 == 0) 1. else -1.
    s* copJoe@absdiPsi(u, theta=copula@parameters, degree=degree, log=log, ...)
})

setMethod("tau", signature("joeCopula"),
          function(copula) tauJoe(theta=copula@parameters))
setMethod("lambda", signature("joeCopula"),
	  function(copula) c(lower=0,
			     upper=copJoe@lambdaU(theta=copula@parameters)))

setMethod("iTau", signature("joeCopula"), # now that tauJoe() is accurate
	  function(copula, tau, tol = 1e-7) {
    if (any(tau < 0)) warning("For the Joe copula, tau must be >= 0. Replacing negative values by 0.")
    copJoe@iTau(tau, tol=tol)
})

## "TODO"
## setMethod("rho", signature("joeCopula"), ... ? ...)
## setMethod("iRho", signature("joeCopula"), function(copula, rho) ...)

## "TODO"
## setMethod("dRho", signature("joeCopula"), ...)
## setMethod("dTau", signature("joeCopula"), ...)

