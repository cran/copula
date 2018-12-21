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


### Marshall--Olkin copula class ###############################################

## Constructor
moCopula <- function(param = NA_real_, dim = 2L)
{
    if(dim != 2) stop("Marshall-Olkin copulas only available in the bivariate case")
    cdf <- quote(min(u1*u2^(1-theta2), u1^(1-theta1)*u2)) # theta1, theta2 in [0,1]
    if(length(param) == 1 && is.na(param)) param <- rep(param, 2)
    dim <- as.integer(dim)
    new("moCopula",
        dimension = dim,
        exprdist = c(cdf = cdf, pdf = expression()), # used in ./dC-dc.R (?)
        parameters = param,
        param.names = paste0("theta", 1:dim),
        param.lowbnd = rep(0, dim),
        param.upbnd = rep(1, dim),
        fullname = "<deprecated slot>")
}


### Methods ####################################################################

## describe method
setMethod(describeCop, c("moCopula", "character"), function(x, kind, prefix="", ...) {
    kind <- match.arg(kind)
    if(kind == "very short") # e.g. for show() which has more parts
        return(paste0(prefix, "Marshall-Olkin copula"))
    ## else
    d <- dim(x)
    ch <- paste0(prefix, "Marshall-Olkin copula, dim. d = ", d)
    switch(kind <- match.arg(kind),
           short = ch,
           long = ch,
           stop("invalid 'kind': ", kind))
})

## dCopula method
setMethod("dCopula", signature("matrix", "moCopula"),
	  function(u, copula, log = FALSE, ...) {
              stop("Marshall-Olkin copulas do not have densities")## Well, have almost everywhere
})

## pCopula method
setMethod("pCopula", signature("matrix", "moCopula"),
	  function(u, copula, log.p = FALSE, ...) {
    theta <- copula@parameters
    if(log.p) l.u <- log(u)
    pmin(if(log.p) l.u[,1] + (1 - theta[2])*l.u[,2] else u[,1] * u[,2]^(1 - theta[2]),
         if(log.p) l.u[,2] + (1 - theta[1])*l.u[,1] else u[,1]^(1 - theta[1]) * u[,2])
})

## rCopula method
setMethod("rCopula", signature("numeric", "moCopula"),
	  function(n, copula) {
    theta <- copula@parameters
    V <- matrix(runif(n * 3), ncol = 3)
    cbind(pmax(V[,1]^(1/(1 - theta[1])), V[,3]^(1/theta[1])),
          pmax(V[,2]^(1/(1 - theta[2])), V[,3]^(1/theta[2])),
          deparse.level=0L)
})

## Measures of association
setMethod("tau", signature("moCopula"), function(copula) {
    theta <- copula@parameters
    theta[1]*theta[2] / (theta[1] - theta[1]*theta[2] + theta[2])
})
setMethod("rho", signature("moCopula"), function(copula) {
    theta <- copula@parameters
    3*theta[1]*theta[2] / (2*theta[1] - theta[1]*theta[2] + 2*theta[2])
})
setMethod("lambda", signature("moCopula"), function(copula)
    c(lower = 0, upper = min(copula@parameters)))
