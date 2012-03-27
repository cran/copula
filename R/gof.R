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


### Goodness-of-fit testing for nested Archimedean copulas

### transformations to univariate quantities ###################################

##' Kendall distribution function
##'
##' @title Kendall distribution function
##' @param u evaluation point(s)
##' @param cop acopula with specified parameter
##' @param d dimension
##' @param n.MC if > 0 a Monte Carlo approach is applied with sample size equal
##'	   to n.MC; otherwise the exact formula is used
##' @return Kendall distribution function at u
##' @author Marius Hofert
K <- function(u, cop, d, n.MC=0)
{
    stopifnot(is(cop, "acopula"))
    th <- cop@theta
    n <- length(u)
    if(n.MC > 0) {
	stopifnot(is.finite(n.MC))
	V <- cop@V0(n.MC,th)
	psiI <- cop@psiInv(u, th)
	unlist(lapply(psiI, function(psInv) mean(ppois(d-1, V* psInv))))
    } else {
	K. <- numeric(n)
	K.[is0 <- u == 0] <- 0
	K.[is1 <- u == 1] <- 1
	if(length(not01 <- seq_len(n)[!(is0 | is1)]))
	    K.[not01] <- if(d == 1) {
		u[not01]
	    } else if(d == 2) {
		u[not01] + cop@psiInv(u[not01], theta=th) / cop@psiInvD1abs(u[not01], th)
	    } else { ## d >= 3 :
		j <- seq_len(d-1)
		lpsiI. <- cop@psiInv(u[not01], theta=th, log=TRUE)
		lpsiDabs <- do.call(rbind,
				    lapply(j, function(j.)
					   cop@psiDabs(exp(lpsiI.),
						       theta=th,
						       degree=j.,
						       log=TRUE))) # (d-1,n)-matrix with n = length(not01)
		lfac.j <- cumsum(log(j)) ## == lfactorial(j)
		r <- colSums(exp(lpsiDabs + j %*% t(lpsiI.) - lfac.j))
		## ensure we are in [0,1] {numerical inaccuracy}
		pmin(1, u[not01] + r)
		## Former code:
		## K2 <- function(psInv) {
		##    lpsiDabs <- unlist(lapply(j, cop@psiDabs,
		##			      u=psInv, theta=th, log=TRUE))
		##    sum(exp(lpsiDabs + j*log(psInv) - lfac.j))
		## }
		## pmin(1, u[not01] + unlist(lapply(psiI[not01], K2)))
		##
		## NB: AMH, Clayton, Frank are numerically not quite monotone near one;
		## --  this does not change that {but maybe slightly *more* accurate}:
		## psiDabs. <- unlist(lapply(j, cop@psiDabs, u = psInv, theta = th,
		##						 log = FALSE))
		##		       sum(psiDabs.*psInv^j/factorial(j))
	    }
	K.
    }
}

##' Transforms supposedly U[0,1]^d distributed vectors of random variates to
##' univariate data (for testing in a one-dimensional setup)
##'
##' @title Transformation to a one-dimensional test setting
##' @param u matrix of random variates to be transformed (typically
##'        pseudo-observations obtained from gtrafo())
##' @param method trafo method. Available are:
##'        "chisq"    : map to a chi-square distribution using the standard normal
##'                     quantile function
##'        "gamma"    : map to an Erlang (= gamma) distribution using the logarithm
##'        "Remillard": the test statistic S_n^{(B)} in Genest, Remillard, Beaudoin (2009)
##'        "Genest"   : the test statistic S_n^{(C)} in Genest, Remillard, Beaudoin (2009)
##' @return the transformed variates
##' @author Marius Hofert and Martin Maechler
## MM: Der Funktionsname passt mir nicht... "ouni" ist hässlich (wieso auch "trafo"? Jede R Funktion ist eine "Trafo" in einem gewissen Sinn..
## Brainstorm   unid21() ("d to 1 dimensional")
gtrafouni <- function(u, method = c("chisq", "gamma", "Remillard", "Genest"))
{
    if(!is.matrix(u)) u <- rbind(u, deparse.level = 0L)
    d <- ncol(u)
    n <- nrow(u)
    method <- match.arg(method)
    switch(method,
	   "chisq" = pchisq(rowSums(qnorm(u)^2), d),
	   "gamma" = pgamma(rowSums(-log(u)), shape=d),
	   "Remillard" = {
	       u <- t(u)
	       sum1 <- exp(lsum(log1p(-u^2)))
	       ## FIXME(MM): should use lsum(log1p(.)) as sum1 ; should work w/o nested apply(.)
	       sum2 <- mean(apply(u, 2, function(u.) { # u. is in R^d
		   sum(apply(u, 2, function(u..) prod(1- pmax(u., u..))))
	       }))
	       n/3^d - sum1/2^(d-1) + sum2
           },
	   "Genest" = {
	       Dn <- apply(u, 1, function(u.){
                   ## u. is one row. We want to know the number of rows of u
                   ## that are (all) componentwise <= u.
                   mean(apply(u <= u., 1, all)) # TRUE <=> the whole row in u is <= u.
               })
               Cperp <- apply(u, 1, prod)
               sum((Dn-Cperp)^2)
           },
 	   stop("gtrafouni: unsupported method ", method))
}


### multivariate transformations ###############################################

##' Transforms vectors of random variates following the given (nested) Archimedean
##' copula (with specified parameters) to U[0,1]^d vectors of random variates
##' via Rosenblatt's transformation
##'
##' @title Rosenblatt transformation for a (nested) Archimedean copula
##' @param u data matrix
##' @param cop an outer_nacopula
##' @param m # order up to which Rosenblatt's transform is computed, i.e.,
##'        C(u_j | u_1,...,u_{j-1}), j=2,..,m
##' @param n.MC parameter n.MC for evaluating the derivatives via Monte Carlo
##' @return matrix of supposedly U[0,1]^d realizations
##' @author Marius Hofert
rtrafo <- function(u, cop, m=d, n.MC=0)
{
    stopifnot(is(cop, "outer_nacopula"), 2 <= m, m <= d)
    if(length(cop@childCops))
        stop("currently, only Archimedean copulas are provided")
    if(!is.matrix(u)) u <- rbind(u, deparse.level = 0L)
    stopifnot(0 <= u, u <=1)
    cop <- cop@copula
    th <- cop@theta
    stopifnot(cop@paraConstr(th))
    dim. <- dim(u)
    n <- dim.[1]
    d <- dim.[2]
    psiI <- cop@psiInv(u, theta=th)
    psiI. <- t(apply(psiI, 1, cumsum))
    ## compute all conditional probabilities
    if(n.MC==0){
	C.j <- function(j){ # computes C(u_j | u_1,...,u_{j-1}) with the same idea as for cacopula
            logD <- cop@psiDabs(c(psiI.[,j], psiI.[,j-1]), theta=th,
                                degree=j-1, n.MC=0, log=TRUE)
            exp(logD[1:n]-logD[(n+1):(2*n)])
        }
    }else{ # n.MC > 0
	## draw random variates
	V <- cop@V0(n.MC, th)
        C.j <- function(j){ # computes C(u_j | u_1,...,u_{j-1}) with the same idea as for cacopula
            ## use same idea as default method of psiDabsMC
            ## only difference: only draw V's once
            arg <- c(psiI.[,j], psiI.[,j-1])
            iInf <- is.infinite(arg)
            logD <- numeric(2*n)
            logD[iInf] <- -Inf
            if(any(!iInf)) logD[!iInf] <- lsum(-V %*% t(arg[!iInf]) +
                                               (j-1) * log(V) - log(n.MC))
            exp(logD[1:n]-logD[(n+1):(2*n)])
        }
    }
    cbind(u[,1],sapply(2:m, C.j))
}

##' Transforms vectors of random variates following the given (nested) Archimedean
##' copula (with specified parameters) to U[0,1]^d (or U[0,1]^(d-1)) vectors of
##' random variates via the transformation of Hering and Hofert (2011)
##'
##' @title Transformation of Hering and Hofert (2011)
##' @param u data matrix
##' @param cop an outer_nacopula
##' @param include.K boolean indicating whether the last component, K, is also
##'        used (include.K = TRUE)
##' @param n.MC parameter n.MC for K
##' @return matrix of supposedly U[0,1]^d realizations
##' @author Marius Hofert and Martin Maechler
htrafo <- function(u, cop, include.K=TRUE, n.MC=0)
{
    stopifnot(is(cop, "outer_nacopula"))
    if(length(cop@childCops))
	stop("currently, only Archimedean copulas are provided")
    if(!is.matrix(u)) u <- rbind(u, deparse.level = 0L)
    stopifnot((d <- ncol(u)) >= 2,
	      0 <= u, u <= 1)
    ## trafo
    th <- cop@copula@theta
    lpsiI <- cop@copula@psiInv(u, th, log=TRUE) # matrix log(psi^{-1}(u))
    lcumsum <- matrix(unlist(lapply(1:d, function(j) lsum(t(lpsiI[,1:j,
                                                                  drop=FALSE])))),
                      ncol=d)
    u. <- matrix(unlist(lapply(1:(d-1), function(k) exp(k*(lcumsum[,k]-
                                                           lcumsum[,k+1])) )),
		 ncol=d-1) # transformed components (uniform under H_0)
    if(include.K) u. <- cbind(u., K(cop@copula@psi(exp(lcumsum[,d]), th),
				    cop=cop@copula, d=d, n.MC=n.MC))
    u.
}


### Gof wrapper ################################################################

##' Conducts a goodness-of-fit test for the given H0 copula cop based on the
##' (copula) data u
##'
##' @title Goodness-of-fit testing for (nested) Archimedean copulas
##' @param u (copula-)data matrix
##' @param cop outer_nacopula with specified H0 parameters
##' @param n.bootstrap number of bootstrap replications
##' @param estimation.method estimation method, see enacopula
##' @param include.K boolean whether K is included in the transformation
##' @param n.MC parameter n.MC for K
##' @param trafo multivariate goodness-of-fit transformation; available are:
##'        "Hering.Hofert" the transformation of Hering, Hofert (2011)
##'        "Rosenblatt" the transformation of Rosenblatt (1952)
##' @param method univariate goodness-of-fit transformation; available are those
##'        of gtrafouni
##' @param verbose if TRUE, the progress of the bootstrap is displayed
##' @param ... additional arguments to enacopula
##' @return htest object
##' @author Marius Hofert and Martin Maechler
gnacopula <- function(u, cop, n.bootstrap,
		      estimation.method= eval(formals(enacopula)$method),
		      include.K= TRUE, n.MC= 0,
		      trafo= c("Hering.Hofert", "Rosenblatt"),
		      method= eval(formals(gtrafouni)$method),
		      verbose= TRUE, ...)
{
    ## setup
    if(!is.matrix(u)) u <- rbind(u, deparse.level = 0L)
    stopifnot(0 <= u, u <= 1, is(cop, "outer_nacopula"), (d <- ncol(u)) >= 2,
              max(cop@comp) == d, n.bootstrap >= 0, n.MC >= 0)
    if(length(cop@childCops))
	stop("currently, only Archimedean copulas are provided")
    if(n.bootstrap == 0)
	stop("Choose a reasonable number of bootstrap replications or
apply the transformations yourself,  see ?gnacopula.")
    u.name <- deparse(substitute(u))

    ## additional warnings for now
    estimation.method <- match.arg(estimation.method)
    if(estimation.method != "mle"){
	if(estimation.method == "smle") warning("'estimation.method = \"smle\"' may be time-consuming!") else
	warning("Consistency for the chosen estimation.method is not clear. Additionally, numerical problems might appear.")
    }

    ## build multivariate transformation
    trafo <- match.arg(trafo)
    method <- match.arg(method)
    gtrafomulti <-
        switch(trafo,
               "Hering.Hofert" = {
                   function(u, cop) htrafo(u, cop=cop, include.K=include.K, n.MC=n.MC)
               },
               "Rosenblatt" = {
                   function(u, cop) rtrafo(u, cop=cop, n.MC=n.MC)
               },
               stop("invalid 'trafo' argument"))

    ## build test statistic function and string describing the method
    string <- "Bootstrapped test of"
    test.stat <-
	switch(method, # define test statistic (and correct string describing the procedure)
	       "chisq" =,
	       "gamma" = {
		   string <- paste0(string, "Anderson and Darling (with trafo = ",
				    trafo, " and method = ", method, ")")
		   function(x) ad.test(x)$statistic
	       },
	       "Remillard" =,
	       "Genest" = {
		   string <- paste0(string, " ", method," (with trafo = ", trafo, ")")
		   function(x) x
	       },
	       stop("wrong 'method' argument"))

    ## main part --- Bootstrapping ------------------

    ## (1) estimate the parameter by the provided estimation method and
    ##	   define the estimated copula
    if(!is.matrix(u)) u <- rbind(u, deparse.level = 0L)
    theta.hat <- enacopula(u, cop, method=estimation.method, ...)
    cop.hat <- onacopulaL(cop@copula@name, list(theta.hat, 1:d)) # copula with theta.hat

    ## (2) transform the data with the copula with estimated parameter
    u.prime <- gtrafomulti(u, cop=cop.hat) # transformed data in the unit hypercube
    Y <- gtrafouni(u.prime, method=method) # transformed data

    ## (3) conduct the Anderson-Darling test or compute the test statistic (depends on method)
    teststat <- test.stat(Y)

    ## (4) conduct the parametric bootstrap
    theta.hat. <- numeric(n.bootstrap) # vector of estimators
    teststat. <- vector("list", n.bootstrap) # vector of test.stat() results
    if(verbose) pb <- txtProgressBar(max=n.bootstrap, style = 3) # setup progress bar
    for(k in 1:n.bootstrap) {

	## (4.1) sample from the copula with estimated parameter and build
	##	     the corresponding pseudo-observations
	u. <- pobs(rnacopula(nrow(u), cop.hat))

	## (4.2) estimate the parameter by the provided method and define
	##	     the estimated copula
	theta.hat.[k] <- enacopula(u., cop, method=estimation.method, ...)
	cop.hat. <- onacopulaL(cop@copula@name, list(theta.hat.[k], 1:d))

	## (4.3) transform the data with the copula with estimated parameter
	u.prime. <- gtrafomulti(u., cop=cop.hat.)
	Y. <- gtrafouni(u.prime., method=method)

	## (4.4) conduct the Anderson-Darling test or compute the test statistic (depends on method)
	teststat.[[k]] <- test.stat(Y.)

        if(verbose) setTxtProgressBar(pb, k) # update progress bar
    }
    if(verbose) close(pb) # close progress bar

    ## (5) build and return results
    structure(class = "htest",
	      list(p.value= mean(unlist(teststat.) > teststat),
                   statistic = teststat, data.name = u.name,
		   method=string, estimator=theta.hat,
		   bootStats = list(estimator=theta.hat., statistic=teststat.)))

}
