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


#### Wrappers and auxiliary functions for dealing with elliptical (Gauss, t_nu)
#### and Archimedean copulas

if(getRversion() < 2.15)
    paste0 <- function (...) paste(..., sep = "")

### copula objects #############################################################

ac.lnames <- if(packageVersion("copula") >= 0.99)
    copula:::c_longNames else nacopula:::c_longNames


##' Determine the copula class for a given copula object
##'
##' @title Copula class for the given copula object
##' @param cop copula object
##' @return "elliptical" or "nArchimedean" depending on the given copula object
##' @author Marius Hofert
copClass <- function(cop)
{
    cls <- class(cop)
    if(is(cop, "copula") && (cls=="normalCopula" || cls=="tCopula")) "elliptical" # note: there could be other "copula" objects which are not elliptical
    else if(cls=="outer_nacopula") "nArchimedean" # can be Archimedean or nested Archimedean
    else stop("not yet supported copula object")
}

##' Determine the copula family for a given copula object
##'
##' @title Copula family for the given copula object
##' @param cop copula object (either elliptical or (nested) Archimedean)
##' @return family string
##' @author Marius Hofert
copFamily <- function(cop)
{
    cls <- class(cop)
    if(is(cop, "copula")){
        if(cls=="normalCopula") "normal"
        else if(cls=="tCopula") "t"
        else stop("unsupported copula family")
    } else if(cls=="outer_nacopula"){
        cop@copula@name # could be nested or not
    } else stop("not yet supported copula object")
}

##' Determine the copula family for a given copula object
##'
##' @title Copula family for the given copula object
##' @param cop copula object (either elliptical or (nested) Archimedean)
##' @return family string
##' @author Marius Hofert
copFamilyClass <- function(family)
{
    if(family=="normal" || family=="t") "elliptical"
    else if(family %in% ac.lnames ||
            family %in% paste0("opower:", ac.lnames)) "nArchimedean"
    else stop("family ", family, " not yet supported")
}

##' Creating elliptical and Archimedean copula objects
##'
##' @title Creating elliptical and Archimedean copula objects
##' @param family specified family (for the elliptical cops, it's "normal" or "t")
##' @param theta parameter (a number or vector)
##' @param d dimension (can be omitted if nacList is given)
##' @param ... additional arguments:
##'     either 'dispstr', 'df', 'df.fixed' for the elliptical copulas
##'        or  'nacList' for the nested Archimedean {use 'theta' above for Archimedean}
##' @return a copula object
##' @author Marius Hofert and Martin Maechler
copCreate <- function(family, theta, d, ...)
{
    switch(copFamilyClass(family),
           "elliptical"={
               stopifnot(is.numeric(theta))
               ellipCopula(family, param=theta, dim=d, ...)
               ## Note: ... args are:
               ## dispstr="ex" (exchangeable), "ar1" (AR1), "toep" (Toeplitz), "un"
               ##         (unstructured)
               ## df (degrees of freedom)
               ## df.fixed=FALSE (TRUE means that d.o.f. will be considered as fixed
               ##          and won't be estimated)
           },
           "nArchimedean"={
               L <- if(hasArg(nacList)){ # nested Archimedean case => we don't need theta
                   list(...)$nacList
               } else {
                   ## Archimedean case => one parameter
                   if(!is.numeric(theta) || length(theta)!=1)
                       stop("must specify 'nacList' (list)  or 'theta' (numeric)")
                   list(theta, 1:d)
               }
               stopifnot(is.list(L))
               onacopulaL(family, L)
               ## Note: ... with an argument "nacList" can be used to construct
               ## a nested Archimedean copula
           },
           stop("family ", family, " not yet supported"))
}


### Kendall's tau ##############################################################

if(FALSE) { ## MM__FIXME__  "old copula" already has kendallsTau() and calibKendallsTau()
showMethods("calibKendallsTau", incl=TRUE)
showMethods("kendallsTau", incl=TRUE)
}## Rather use these by providing "nacopula" methods !!

##' Determine tau from given theta (matricized)
##'
##' @title Determine tau from given theta
##' @param theta parameter or matrix of pairwise parameters
##' @param family copula family
##' @param ... additional args passed to tau-slots
##' @return Kendall's tau (possibly a matrix)
##' @author Marius Hofert
tau <- function(theta, family, ...)
{
    switch(copFamilyClass(family),
           "elliptical"={
               2*asin(theta)/pi
           },
           "nArchimedean"={
               getAcop(family)@tau(theta, ...)
           },
           stop("family ", family, " not yet supported"))
}

##' Determine theta from given tau (matricized)
##'
##' @title Determine theta from given tau
##' @param tau Kendalls tau  or matrix of pairwise Kendall's tau
##' @param family copula family
##' @param ... additional args passed to tau-slots
##' @return parameter theta (possibly a matrix)
##' @author Marius Hofert
tauInv <- function(tau, family, ...)
{
    switch(copFamilyClass(family),
           "elliptical"={
               sin(pi/2*tau)
           },
           "nArchimedean"={
               getAcop(family)@tauInv(tau, ...)
           },
           stop("family ", family, " not yet supported"))
}


### Sampling ###################################################################

## MM: should call rCopula()

##' Sampling elliptical and (nested) Archimedean copulas
##'
##' @title Sampling elliptical and (nested) Archimedean copulas
##' @param n sample size
##' @param cop copula to sample
##' @return an (n x d)-matrix of random variates following the specified copula cop
##' @author Marius Hofert
rcop <- function(n, cop)
{
    switch(copClass(cop),
           "elliptical"={
               rcopula(cop, n)
           },
           "nArchimedean"={
               rnacopula(n, cop)
           },
           stop("not yet supported copula object"))
}


### Densities for elliptical copulas ###########################################

##' Density for elliptical copulas
##'
##' @title Density for elliptical copulas
##' @param u data matrix (in \eqn{[0,1]^d})
##' @param family elliptical family
##' @param P correlation matrix P
##' @param log logical determining if log(density(.)) should be returned
##' @param df degree of freedom parameter (\eqn{\nu}) for t-copulas
##' @param ... additional arguments passed to \code{\link[mvtnorm]{dmvnorm}}
##' or \code{\link[mvtnorm]{dmvt}}.
##' @return density of the specified copula evaluated at \code{u}.
##' @author Marius Hofert (and MMa)
dellip <- function(u, family, P, log=FALSE, df, ...)
{
    ## We assume that this will be part of 'copula' which has dmvt() etc in its NAMESPACE
    if(FALSE)## _OR_ that the caller of this function has executed
        require("mvtnorm")# typically faster than using mvtnorm::* all the time
    val <-
	switch(family,
	       "normal" =
	   {
	       qnu <- qnorm(u)
	       dmvnorm(qnu, sigma=P, log=TRUE) - rowSums(dnorm(qnu, log=TRUE))
	   },
	       "t" =
	   {
	       qtu <- qt(u, df=df)
	       ## Note: for dmvt, log=TRUE is actually the default;
               ##       furthermore, delta=rep(0, length=ncol(u)) is the default
               ##       when delta is missing (although not mentioned on ?dmvt)
	       dmvt(qtu, sigma=P, df=df, log=TRUE) - rowSums(dt(qtu, df=df, log=TRUE))
	   },
	       stop("family ", family, " not yet supported"))
    if(log) val else exp(val)
}
