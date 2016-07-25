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



##' @title Transformation of Hering and Hofert (2014) or Its Inverse
##' @param u data matrix in [0,1]^(n, d)
##' @param copula an 'outer_nacopula' or 'archmCopula' object
##' @param include.K boolean indicating whether the last component, K, is also
##'        used (include.K = TRUE); ignored when inverse=TRUE since K is crucial
##'        there
##' @param n.MC parameter n.MC for K
##' @param inverse logical indicating whether the inverse of htrafo is computed
##'        (this transformation can also be found in Wu, Valdez, Sherris (2006)).
##' @param method method to compute qK() (see there)
##' @param u.grid argument of qK() (for method "discrete")
##' @param ... additional arguments passed to qK() (see there)
##' @return matrix of transformed realizations
##' @author Marius Hofert and Martin Maechler
htrafo <- function(u, copula, include.K = TRUE, n.MC = 0, inverse = FALSE,
                   method = eval(formals(qK)$method), u.grid, ...)
{
    ## Checks
    if(is(copula, "outer_nacopula")) {
	if(length(copula@childCops))
	    stop("Currently, only Archimedean copulas are supported")
	## outer_nacopula but with no children => an AC => continue
	copula <- copula@copula # class(copula) = "acopula"
	th <- copula@theta
    } else if(is(copula, "archmCopula")) {
	th <- copula@parameters
	d <- copula@dimension
	fam.name <- getAname(copula)
	copula <- onacopulaL(fam.name, nacList = list(th, 1:d))@copula
	## copula is an "acopula" *with* dimension and parameter (as required below)
    } else {
        stop("Not yet implemented for copula class ", class(copula))
    }
    if(!is.matrix(u)) u <- rbind(u, deparse.level = 0L)
    stopifnot((d <- ncol(u)) >= 2, 0 <= u, u <= 1)

    ## Trafos
    if(inverse){ # "simulation trafo" of Wu, Valdez, Sherris (2006)

        ## ingredient 1: log(psi^{-1}(K^{-1}(u_d)))
        KI <- qK(u[,d], copula = copula, d = d, n.MC = n.MC, method = method,
                 u.grid = u.grid, ...) # n-vector K^{-1}(u_d)
        lpsiIKI <- copula@iPsi(KI, th, log = TRUE) # n-vector log(psi^{-1}(K^{-1}(u_d)))
        n <- nrow(u)
        ## ingredient 2: sum_{k=j}^{d-1} log(u_k)/k) for j=1,..,d
        lu. <- log(u[,-d, drop = FALSE]) * (ik <- 1/rep(1:(d-1), each = n))
        cslu. <- apply(lu.[,(d-1):1, drop = FALSE], 1, cumsum) ## note: we apply cumsum to reversed columns,
        ## because we need the "upper partial sums"
        cslu. <- if(d == 2) as.matrix(cslu.) else t(cslu.) # get a result of the right dimensions
        cslu <- cbind(cslu.[,(d-1):1, drop = FALSE], 0) # revert the column order, and bind last column to it
        ## => n x d matrix
        ## ingredient 3:
        l1p <- cbind(0, log1p(-u[,1:(d-1), drop = FALSE]^ ik)) # n x (d-1) matrix + dummy 0's in the first col
        ## finally, compute the transformation
        expo <- rep(lpsiIKI, d) + cslu + l1p
        copula@psi(exp(expo), th)

    } else { # "goodness-of-fit trafo" of Hering and Hofert (2014)

        lpsiI <- copula@iPsi(u, th, log = TRUE) # matrix log(psi^{-1}(u))
        lcumsum <- matrix(unlist(lapply(1:d, function(j)
                                        lsum(t(lpsiI[,1:j, drop = FALSE])))), ncol = d)
        u. <- matrix(unlist(lapply(1:(d-1),
                                   function(k) exp(k*(lcumsum[,k]-
                                                      lcumsum[,k+1])) )),
                     ncol = d-1) # transformed components (uniform under H_0)
	if(include.K) u. <- cbind(u., pK(copula@psi(exp(lcumsum[,d]), th),
                                         copula = copula, d = d, n.MC = n.MC))
        u.

    }
}


### Deprecated stuff ###########################################################

rtrafo <- function(u, copula, indices = 1:dim(copula), inverse = FALSE, log = FALSE)
{
    .Deprecated("cCopula")
    cCopula(u, copula = copula, indices = indices, inverse = inverse, log = log)
}


##' Conducts a goodness-of-fit test for the given H0 copula cop based on the
##' (copula) data u
##'
##' @title Goodness-of-fit testing for (nested) Archimedean copulas
##' @param u (copula-)data matrix
##' @param cop outer_nacopula with specified H0 parameters
##' @param n.bootstrap number of bootstrap replications
##' @param estim.method estimation method, see enacopula
##' @param include.K boolean whether K is included in the transformation
##' @param n.MC parameter n.MC for K
##' @param trafo multivariate goodness-of-fit transformation; available are:
##'        "Hering.Hofert" the transformation of Hering, Hofert (2014)
##'        "Rosenblatt" the transformation of Rosenblatt (1952)
##' @param method test statistic for the test of U[0,1]^d; see gofTstat()
##' @param verbose if TRUE, the progress of the bootstrap is displayed
##' @param ... additional arguments to enacopula
##' @return htest object
##' @author Marius Hofert and Martin Maechler
gnacopula <- function(u, cop, n.bootstrap,
		      estim.method=eval(formals(enacopula)$method),
		      include.K=TRUE, n.MC=0,
		      trafo=c("Hering.Hofert", "Rosenblatt"),
		      method=eval(formals(gofTstat)$method),
		      verbose=TRUE, ...)
{
    .Deprecated("gofCopula")
    gofCopula(cop, x=u, N=n.bootstrap, method=method,
              estim.method=estim.method,
              verbose=verbose, n.MC=n.MC, if(trafo == "htrafo")
              include.K=include.K)

## working (but deprecated)

##     ## setup
##     if(!is.matrix(u)) u <- rbind(u, deparse.level = 0L)
##     stopifnot(0 <= u, u <= 1, is(cop, "outer_nacopula"), (d <- ncol(u)) >= 2,
##               max(cop@comp) == d, n.bootstrap >= 0, n.MC >= 0)
##     if(length(cop@childCops))
## 	stop("currently, only Archimedean copulas are provided")
##     if(n.bootstrap == 0)
## 	stop("Choose a reasonable number of bootstrap replications or
## apply the transformations yourself,  see ?gnacopula.")
##     u.name <- deparse(substitute(u))

##     ## additional warnings for now
##     estim.method <- match.arg(estim.method)
##     if(estim.method != "mle"){
## 	if(estim.method == "smle") warning("'estim.method = \"smle\"' may be time-consuming!") else
## 	warning("Consistency for the chosen estim.method is not clear. Additionally, numerical problems might appear.")
##     }

##     ## build multivariate transformation
##     trafo <- match.arg(trafo)
##     method <- match.arg(method)
##     gtrafomulti <-
##         switch(trafo,
##                "Hering.Hofert" = {
##                    function(u, cop) htrafo(u, copula = cop, include.K=include.K, n.MC=n.MC)
##                },
##                "Rosenblatt" = {
##                    function(u, cop) rtrafo(u, cop=cop, n.MC=n.MC)
##                },
##                stop("invalid 'trafo' argument"))

##     ## build test statistic function and 'meth' string describing the method
##     meth <- paste0("Bootstrapped (B =", n.bootstrap,") test of ")
##     meth2 <- paste0(method,", est.method = ", estim.method)
##     meth <-
## 	switch(method,
##                "Sn" = {
##                    paste0(meth, meth2)
##                },
##                "SnB" =, "SnC" = {
## 		   paste0(meth, meth2," (with trafo = ", trafo, ")")
## 	       },
## 	       "AnChisq" =, "AnGamma" = {
## 		   paste0(meth, "Anderson and Darling (with trafo = ",
##                           trafo, " and method = ", meth2, ")")
## 	       },
## 	       stop("wrong 'method' argument"))

##     ## main part --- Bootstrapping ------------------

##     ## (1) estimate the parameter by the provided estimation method and
##     ##	   define the estimated copula
##     theta.hat <- enacopula(u, cop, method=estim.method, ...)
##     cop.hat <- onacopulaL(cop@copula@name, list(theta.hat, 1:d)) # copula with theta.hat

##     ## (2) transform the data with the copula with estimated parameter
##     u.prime <- gtrafomulti(u, cop=cop.hat) # transformed data in the unit hypercube

##     ## (3) conduct the Anderson-Darling test or compute the test statistic (depends on method)
##     T <- if(method=="Sn") gofTstat(u.prime, method=method, copula=cop.hat)
##          else gofTstat(u.prime, method=method)

##     ## (4) conduct the parametric bootstrap
##     theta.hat. <- numeric(n.bootstrap) # vector of estimators
##     T. <- numeric(n.bootstrap)# vector of gofTstat() results
##     if(verbose) {	     # setup progress bar and close it on exit
## 	pb <- txtProgressBar(max = n.bootstrap, style = if(isatty(stdout())) 3 else 1)
## 	on.exit(close(pb))
##     }
##     for(k in 1:n.bootstrap) {

## 	## (4.1) sample from the copula with estimated parameter and build
## 	##	     the corresponding pseudo-observations
## 	u. <- pobs(rnacopula(nrow(u), cop.hat))

## 	## (4.2) estimate the parameter by the provided method and define
## 	##	     the estimated copula
## 	theta.hat.[k] <- enacopula(u., cop, method=estim.method, ...)
## 	cop.hat. <- onacopulaL(cop@copula@name, list(theta.hat.[k], 1:d))

## 	## (4.3) transform the data with the copula with estimated parameter
## 	u.prime. <- gtrafomulti(u., cop=cop.hat.)

## 	## (4.4) compute the test statistic
##         T.[k] <- if(method=="Sn") gofTstat(u.prime., method=method, copula=cop.hat.)
##         else gofTstat(u.prime., method=method)
##         if(verbose) setTxtProgressBar(pb, k) # update progress bar
##     }

##     ## (5) build and return results
##     structure(class = "htest",
## 	      list(p.value= (sum(T. > T) + 0.5)/(n.bootstrap+1),
##                    statistic = T, data.name = u.name,
## 		   method=meth, estimator=theta.hat,
## 		   bootStats = list(estimator=theta.hat., statistic=T.)))
}
