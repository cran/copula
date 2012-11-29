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

### test statistics ############################################################

##' Test statistics for various goodness-of-fit tests of (supposedly) U[0,1]^d
##' distributed vectors of random variates
##'
##' @title Test statistics for U[0,1]^d
##' @param u matrix of supposedly  U[0,1]^d observations
##' @param method various test statistics. Available are:
##'        "AnChisq": Anderson-Darling test statistic after map to a chi-square
##'                   distribution using the standard normal quantile function
##'        "AnGamma": Anderson-Darling test statistic after map to an Erlang/Gamma
##'                   distribution using the logarithm
##'        "SnB"  : the test statistic S_n^{(B)} in Genest, Remillard, Beaudoin (2009)
##'        "SnC"  : the test statistic S_n^{(C)} in Genest, Remillard, Beaudoin (2009)
##' @return values of the chosen test statistic
##' @author Marius Hofert and Martin Maechler
gofTstat <- function(u, method = c("AnChisq", "AnGamma", "SnB", "SnC"))
{
    if(!is.matrix(u)) u <- rbind(u, deparse.level = 0L)
    d <- ncol(u)
    n <- nrow(u)
    method <- match.arg(method)
    switch(method,
	   "AnChisq" = ad.test( pchisq(rowSums(qnorm(u)^2), d) )$statistic,
	   "AnGamma" = ad.test( pgamma(rowSums(-log(u)), shape=d) )$statistic,
	   "SnB" =
       { ## S_n(B)
           lu2 <- log1p(-u^2) # n x d matrix of log(1-u_{ij}^2)
           ## Note (modulo rowSums/colSums):
           ## Idea: sum1 = sum(prod(1-u^2)) = sum(exp(sum(lu2)))
           ## = exp(log( sum(exp(rowSums(lu2))) )) = exp(lsum(rowSums(lu2)))
           slu2 <- rowSums(lu2) # vector of length n
	   sum1 <- exp(lsum(matrix(slu2, ncol=1))) # lsum() needs a matrix; result: 1 value
           ## 2) The notation here is similar to Genest, Remillard,
           ## Beaudoin (2009) but we interchange k and j (since j always runs
           ## in 1:d). That being said...
	   lu <- t(log1p(-u)) # t(n x d matrix of log(1-u_{ij})) --> accessing columns
           ln <- log(n)
           ## Idea:
           ##   1/n sum_i sum_k prod_j (1-max(u_{ij},u_{kj}))
           ## = 1/n sum_i sum_k exp( sum_j log(1-max{.,.}) )
           ## = 1/n sum_i sum_k exp( sum_j log(min{1-u_{ij},1-u_{kj}}) )
           ## = 1/n sum_i sum_k exp( sum_j min{ log(1-u_{ij}), log(1-u_{kj}) })
           ## = 1/n sum_i sum_k exp( sum(pmin{ lu[i,], lu[k,]}) )
           ## = 1/n sum_i exp( log(sum_k exp( sum(pmin{ lu[i,], lu[k,]}) )) )
           ## = 1/n sum_i exp( lsum( sum(pmin{ lu[i,], lu[k,]}) ) )
           ## = sum_i exp(-log(n) + lsum( sum(pmin{ lu[i,], lu[k,]}) ))
           ## = sum_i exp(-log(n) + lsum_{over k in 1:n}( sum(pmin{ lu[i,], lu[k,]}) ))
           ## => for each fixed i, (l)apply lsum()
	   sum2mands <- unlist(lapply(1:n, function(i){
	       lu.i <- lu[,i] ## now i is fixed
	       sum.k <- vapply(1:n, function(k)# sum over k (n-dim. vector)
			       sum(pmin(lu.i, lu[,k])), 0.)
	       ls.i <- lsum(matrix(sum.k, ncol=1)) # lsum( sum(pmin(...)) ) for fixed i; 1 value
	       exp(-ln + ls.i)
	   }))
	   n/3^d - sum1/2^(d-1) + sum(sum2mands)
       },
	   "SnC" =
       { ## S_n(C)
	   Dn <- apply(u, 1, function(u.){ # Dn is a vector of length n
	       ## u. is one row. We want to know the number of rows of u
	       ## that are (all) componentwise <= u.
	       mean(apply(t(u) <= u., 2, all)) # TRUE <=> the whole row in u is <= u.
	   })
           Cperp <- apply(u, 1, prod)
	   sum((Dn-Cperp)^2)
       },
	   stop("unsupported method ", method))
}


### multivariate transformations ###############################################

##' Transforms vectors of random variates following the given (nested) Archimedean
##' copula (with specified parameters) to U[0,1]^d vectors of random variates
##' via Rosenblatt's transformation
##'
##' @title Rosenblatt transformation for a (nested) Archimedean copula
##' @param u matrix of (pseudo-/copula-)observations
##' @param cop object of class Copula
##' @param j.ind indices j >= 2 for which Rosenblatt's transform is computed, i.e.,
##'        C(u_j | u_1,...,u_{j-1})
##' @param m (deprecated)
##' @param n.MC parameter n.MC for evaluating the derivatives via Monte Carlo
##' @param log
##' @param trafo.only
##' @return matrix U (n x k) of supposedly U[0,1]^k realizations, where k=1+length(j.ind);
##'         where U[,1] == u[,1] in any case.
##' @author Marius Hofert and Martin Maechler
rtrafo <- function(u, cop, j.ind=2:d, m, n.MC=0, log=FALSE, trafo.only=log)
{
    if(!is.matrix(u)) u <- rbind(u, deparse.level = 0L)
    stopifnot(is(cop, "Copula"))
    if(!missing(m)) {
	warning("'m' is deprecated; rather specify 'j.ind = 2:m'")
	stopifnot(m >= 2)
	j.ind <- 2:m
    }
    d. <- dim(u)
    n <- d.[1]
    d <- d.[2]
    stopifnot(0 <= u, u <= 1,
	      2 <= j.ind, j.ind <= d)

    if((nac <- is(cop, "outer_nacopula")) ||
       is(cop, "archmCopula")) {
	if(nac) {
	    if(length(cop@childCops))
		stop("currently, only Archimedean copulas are provided")
	    cop <- cop@copula
	    th <- cop@theta
	} else { # "archmCopula"
	    ## "claytonCopula", "frankCopula", "gumbelCopula", ...
	    th <- cop@parameters
	    ## ccop <- cop
	    cop <- getAcop(cop)
	}
	stopifnot(cop@paraConstr(th))
	psiI <- cop@iPsi(u, theta=th)	   # n x d
	psiI. <- t(apply(psiI, 1, cumsum)) # n x d
	## compute all conditional probabilities
	if(n.MC==0){
	    ## Note: C(u_j | u_1,...,u_{j-1}) = \psi^{(j-1)}(\sum_{k=1}^j \psi^{-1}(u_k)) / \psi^{(j-1)}(\sum_{k=1}^{j-1} \psi^{-1}(u_k))
	    C.j <- function(j) { # computes C(u_j | u_1,...,u_{j-1}) with the
					# same idea but faster than for cacopula()
		logD <- cop@absdPsi(as.vector(psiI.[,c(j,j-1)]), theta=th,
				    degree=j-1, n.MC=0, log=TRUE)
		res <- logD[1:n]-logD[(n+1):(2*n)]
		if(log) res else exp(res)
	    }
	} else {			  # n.MC > 0
	    ## draw random variates
	    V <- cop@V0(n.MC, th)
	    C.j <- function(j){ # computes C(u_j | u_1,...,u_{j-1}) ... see above
		## use same idea as default method of absdPsiMC
		## only difference: only draw V's once
		arg <- c(psiI.[,j], psiI.[,j-1])
		iInf <- is.infinite(arg)
		logD <- numeric(2*n)
		logD[iInf] <- -Inf
		if(any(!iInf)) logD[!iInf] <- lsum(-V %*% t(arg[!iInf]) +
						   (j-1) * log(V) - log(n.MC))
		res <- logD[1:n]-logD[(n+1):(2*n)]
		if(log) res else exp(res)
	    }
	}
    } else if(is(cop, "normalCopula")) {
	sigma <- getSigma(cop)
	n <- nrow(X <- qnorm(u))
	C.j <- function(j) { ## j is the dimension -- really a function of (j, X := F(u))
	    ## log(Determinant) :
	    stopifnot(j >= 2, (dd <- determinant(Sd <- sigma[1:j,1:j]))$sign >= 0)
	    logDet <- dd$modulus
	    ## invsigma is the inverse matrix of sigma[1:j,1:j]
	    invsigma <- solve(Sd)

	    x <- X[,1:j, drop=FALSE]
	    x. <- x[,1:(j-1), drop=FALSE]
	    IS. <- invsigma[-j,-j]
	    sid <- sqrt(invsigma[j,j])

	    ## sum1 = \sum_{i=1}^{j-1}(a_{id}+a_{di})x_i
	    sum1 <- colSums( (invsigma[j,1:(j-1)]+invsigma[1:(j-1),j]) * t(x.) )
	    H <- sid*x[,j] + sum1/(2*sid)

	    ## calculate A
	    ## sum2 = \sum_{j=1}^{j-1}\sum_{i=1}^{j-1}a_{ij} x_i x_j
	    sum2 <- if (j==2) invsigma[1,1]*x.^2 else sapply(1:n,
			function(i) crossprod(x.[i,], IS. %*% x.[i,]))

	    A <- sum2 - sum1^2/(4*invsigma[j,j])

	    ## term= (2\pi)^{(j-1)/2} |\Sigma_j|^{1/2} \sqrt{a_{dd}}
	    term <-  if(log) (j-1)/2*log(2*pi) + logDet/2 + log(sid)
	    else (2*pi)^((j-1)/2)*exp(logDet/2)*sid
	    Id <- if(log) -0.5*A + pnorm(H, log.p=TRUE) - term
	    else exp(-0.5*A)*pnorm(H)/term

	    ## density of Gaussian copula with dimension j-1
	    two <- if (j==2) dnorm(x., log=log) else {
		dmvnorm(x., sigma=sigma[1:(j-1),1:(j-1)], log=log)
	    }
	    if(log) Id - two else (Id/two)

	} ## C.j

    } else if(is(cop, "tCopula")) {
	sigma <- getSigma(cop)
	df <- getdf(cop)
	n <- nrow(X <- qt(u, df=df))
	C.j <- function(j) { ## j is the dimension -- really a function of (j, df, X := F(u))
	    ## log(Determinant) :
	    stopifnot(j >= 2, (dd <- determinant(Sd <- sigma[1:j,1:j]))$sign >= 0)
	    logDet <- dd$modulus
	    ## invsigma is the inverse matrix of sigma[1:j,1:j]
	    invsigma <- solve(Sd)

	    x <- X[,1:j, drop=FALSE]
	    x. <- x[,1:(j-1), drop=FALSE]
	    IS. <- invsigma[-j,-j]
	    sid <- sqrt(invsigma[j,j])

	    ## sum1 = \sum_{i=1}^{j-1}(a_{id}+a_{di}) x_i
	    sum1 <-  colSums( (invsigma[j,1:(j-1)]+invsigma[1:(j-1),j]) * t(x.)	)

	    ## calculate A
	    ## sum2 = \sum_{j=1}^{j-1}\sum_{i=1}^{j-1}a_{ij} x_i x_j
	    sum2 <- if (j==2) invsigma[1,1]*x.^2 else sapply(1:n,
			function(n) crossprod(x.[n,], IS. %*% x.[n,]))

	    A <- sum2-sum1^2/(4*invsigma[j,j])

	    f <- function(z,i) {
		(1 + ((sid*z+sum1[i]/(2*sid))^2 + A[i])/df)^((-df-j)/2)
	    }

	    ## term=\frac{\Gamma{[(\nu+j)/2]}}{\Gamma{(\nu/2)}(\nu\pi)^{j/2}|\Sigma_j|^{1/2}}
	    lterm <- lgamma((df+j)/2)-(lgamma(df/2)+log(pi*df)*(j/2)+0.5*logDet)
	    Id <- sapply(1:n, function(i)
			 integrate(f, lower= -Inf, upper= x[i,j], i=i)$value)
            Id <- if(log) lterm + log(Id) else exp(lterm) * Id

	    ## density of Student-t copula with dimension j-1and degrees of freedom df
	    two <- if(j==2) dt(x.,df=df, log=log) else {
		## require(mnormt)
		## dmt(x. ,df=df, mean = rep(0, (d-1)),S=sigma[1:(d-1),1:(d-1)])
		dmvt(x., df=df, sigma=sigma[1:(j-1),1:(j-1)], log=log)
	    }
	    if(log) pmin(0, Id - two) else pmin(1, Id/two)

	} ## C.j

    } else {
	stop("not yet implemented for copula class ", class(cop))
    }

    if(trafo.only)
	vapply(j.ind, C.j, numeric(n))
    else {
	u[, -1] <- vapply(j.ind, C.j, numeric(n))
	u
    }
}

##' Transforms vectors of random variates following the given (nested) Archimedean
##' copula (with specified parameters) to [0,1]^d (or [0,1]^(d-1)) vectors of
##' random variates via the transformation of Hering and Hofert (2011) or its
##' inverse
##'
##' @title Transformation of Hering and Hofert (2011) or its inverse
##' @param u data matrix in [0,1]^d
##' @param cop an outer_nacopula
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
htrafo <- function(u, cop, include.K=TRUE, n.MC=0, inverse=FALSE,
                   method = formals(qK)$method, u.grid, ...)
{
    ## checks
    stopifnot(is(cop, "outer_nacopula"))
    if(length(cop@childCops))
        stop("currently, only Archimedean copulas are provided")
    if(!is.matrix(u)) u <- rbind(u, deparse.level = 0L)
    stopifnot((d <- ncol(u)) >= 2,
	      0 <= u, u <= 1)
    ## trafos
    th <- cop@copula@theta
    if(inverse){ # "simulation trafo" of Wu, Valdez, Sherris (2006)
        ## ingredient 1: log(psi^{-1}(K^{-1}(u_d)))
        KI <- qK(u[,d], cop=cop@copula, d=d, n.MC=n.MC, method=method,
                 u.grid=u.grid, ...) # n-vector K^{-1}(u_d)
        lpsiIKI <- cop@copula@iPsi(KI, th, log=TRUE) # n-vector log(psi^{-1}(K^{-1}(u_d)))
        n <- nrow(u)
        ## ingredient 2: sum_{k=j}^{d-1} log(u_k)/k) for j=1,..,d
        lu. <- log(u[,-d, drop=FALSE]) * (ik <- 1/rep(1:(d-1), each=n))
        cslu. <- apply(lu.[,(d-1):1, drop=FALSE], 1, cumsum) ## note: we apply cumsum to reversed columns,
        ## because we need the "upper partial sums"
        cslu. <- if(d==2) as.matrix(cslu.) else t(cslu.) # get a result of the right dimensions
        cslu <- cbind(cslu.[,(d-1):1, drop=FALSE], 0) # revert the column order, and bind last column to it
        ## => n x d matrix
        ## ingredient 3:
        l1p <- cbind(0, log1p(-u[,1:(d-1), drop=FALSE]^ ik)) # n x (d-1) matrix + dummy 0's in the first col
        ## finally, compute the transformation
        expo <- rep(lpsiIKI, d) + cslu + l1p
        cop@copula@psi(exp(expo), th)
    } else { # "goodness-of-fit trafo" of Hofert and Hering (2011)
        lpsiI <- cop@copula@iPsi(u, th, log=TRUE) # matrix log(psi^{-1}(u))
        lcumsum <- matrix(unlist(lapply(1:d, function(j)
                                        lsum(t(lpsiI[,1:j, drop=FALSE])))),
                          ncol=d)
        u. <- matrix(unlist(lapply(1:(d-1),
                                   function(k) exp(k*(lcumsum[,k]-
                                                      lcumsum[,k+1])) )),
                     ncol=d-1) # transformed components (uniform under H_0)
	if(include.K) u. <- cbind(u., pK(cop@copula@psi(exp(lcumsum[,d]), th),
                                         cop=cop@copula, d=d, n.MC=n.MC))
        u.
    }
}


### Gof wrapper ################################################################

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
##'        "Hering.Hofert" the transformation of Hering, Hofert (2011)
##'        "Rosenblatt" the transformation of Rosenblatt (1952)
##' @param method test statistic for the test of U[0,1]^d; see gofTstat()
##' @param verbose if TRUE, the progress of the bootstrap is displayed
##' @param ... additional arguments to enacopula
##' @return htest object
##' @author Marius Hofert and Martin Maechler
gnacopula <- function(u, cop, n.bootstrap,
		      estim.method= eval(formals(enacopula)$method),
		      include.K= TRUE, n.MC= 0,
		      trafo= c("Hering.Hofert", "Rosenblatt"),
		      method= eval(formals(gofTstat)$method),
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
    estim.method <- match.arg(estim.method)
    if(estim.method != "mle"){
	if(estim.method == "smle") warning("'estim.method = \"smle\"' may be time-consuming!") else
	warning("Consistency for the chosen estim.method is not clear. Additionally, numerical problems might appear.")
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

    ## build test statistic function and 'meth' string describing the method
    meth <- paste0("Bootstrapped (B =", n.bootstrap,") test of ")
    meth2 <- paste0(method,", est.method = ", estim.method)
    meth <-
	switch(method,
	       "AnChisq" =,  ## A_n
	       "AnGamma" = { ## A_n with gamma distribution
		   paste0(meth, "Anderson and Darling (with trafo = ",
                          trafo, " and method = ", meth2, ")")
	       },
	       "SnB" =, ## S_n(B) and S_n(C)
	       "SnC" = {
		   paste0(meth, meth2," (with trafo = ", trafo, ")")
	       },
	       stop("wrong 'method' argument"))

    ## main part --- Bootstrapping ------------------

    ## (1) estimate the parameter by the provided estimation method and
    ##	   define the estimated copula
    theta.hat <- enacopula(u, cop, method=estim.method, ...)
    cop.hat <- onacopulaL(cop@copula@name, list(theta.hat, 1:d)) # copula with theta.hat

    ## (2) transform the data with the copula with estimated parameter
    u.prime <- gtrafomulti(u, cop=cop.hat) # transformed data in the unit hypercube
    ## (3) conduct the Anderson-Darling test or compute the test statistic (depends on method)
    T <- gofTstat(u.prime, method=method)  # transformed data

    ## (4) conduct the parametric bootstrap
    theta.hat. <- numeric(n.bootstrap) # vector of estimators
    T. <- numeric(n.bootstrap)# vector of gofTstat() results
    if(verbose) {	     # setup progress bar and close it on exit
	pb <- txtProgressBar(max = n.bootstrap, style = if(isatty(stdout())) 3 else 1)
	on.exit(close(pb))
    }
    for(k in 1:n.bootstrap) {

	## (4.1) sample from the copula with estimated parameter and build
	##	     the corresponding pseudo-observations
	u. <- pobs(rnacopula(nrow(u), cop.hat))

	## (4.2) estimate the parameter by the provided method and define
	##	     the estimated copula
	theta.hat.[k] <- enacopula(u., cop, method=estim.method, ...)
	cop.hat. <- onacopulaL(cop@copula@name, list(theta.hat.[k], 1:d))

	## (4.3) transform the data with the copula with estimated parameter
	u.prime. <- gtrafomulti(u., cop=cop.hat.)

	## (4.4) compute the test statistic
        T.[k] <- gofTstat(u.prime., method=method)
        if(verbose) setTxtProgressBar(pb, k) # update progress bar
    }

    ## (5) build and return results
    structure(class = "htest",
	      list(p.value= (sum(T. > T) + 0.5)/(n.bootstrap+1),
                   statistic = T, data.name = u.name,
		   method=meth, estimator=theta.hat,
		   bootStats = list(estimator=theta.hat., statistic=T.)))
}
