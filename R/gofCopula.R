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


### Various goodness-of-fit test statistics and tests ##########################


### Two-sample test statistic of Remillard, Scaillet (2009) ####################

##' @title Test Statistic of Remillard, Scaillet (2009, "Testing for equality
##'        between two copulas")
##' @param u1 (n1, d)-sample of copula observations (in [0,1]^d)
##' @param u2 (n2, d)-sample of copula observations (in [0,1]^d)
##' @param useR logical indicating whether R or C implementations are used
##' @return value of the test statistic
##' @author Marius Hofert
##' @note - See p. 3 in Remillard, Scaillet (2009)
##'       - R version slow for larger n1, n2 or d
gofT2stat <- function(u1, u2, useR = FALSE)
{
    ## Checks, basic variables etc.
    if(!is.matrix(u1)) u1 <- rbind(u1)
    if(!is.matrix(u2)) u2 <- rbind(u2)
    d1 <- ncol(u1)
    d2 <- ncol(u2)
    if(d1 != d2)
        stop("'u1' and 'u2' must be of the same dimension (same number of columns)")
    ## Main
    if(useR)
    {
        n1 <- nrow(u1)
        n2 <- nrow(u2)
        ## Part 1: u1 with u1
        res1 <- sum(vapply(seq_len(n1), function(i) {
            sum(vapply(seq_len(n1), function(k) {
                prod(1-pmax(u1[i,], u1[k,]))
            }, numeric(1))
            )}, numeric(1)))
        ## Part 2: u1 with u2
        res2 <- sum(vapply(seq_len(n1), function(i) {
            sum(vapply(seq_len(n2), function(k) {
                prod(1-pmax(u1[i,], u2[k,]))
            }, numeric(1))
            )}, numeric(1)))
        ## Part 3: u2 with u2
        res3 <- sum(vapply(seq_len(n2), function(i) {
            sum(vapply(seq_len(n2), function(k) {
                prod(1-pmax(u2[i,], u2[k,]))
            }, numeric(1))
            )}, numeric(1)))
        ## Return
        (res1/n1^2 - 2*res2/(n1*n2) + res3/n2^2) / (1/n1 + 1/n2)
    } else {
        .Call(gofT2stat_c, u1, u2)
    }
}


### Test statistics for gofCopula()-based tests ################################

##' @title Test Statistics for Tests of U[0,1]^d
##' @param u (n, d)-matrix of supposedly U[0,1]^d observations
##' @param method various test statistics. Available are:
##'        "Sn"     : the test statistic S_n (Cramer-von Mises) in Genest, Remillard, Beaudoin (2009)
##'        "SnB"    : the test statistic S_n^{(B)} in Genest, Remillard, Beaudoin (2009)
##'        "SnC"    : the test statistic S_n^{(C)} in Genest, Remillard, Beaudoin (2009)
##'        "AnChisq": Anderson-Darling test statistic after map to a chi-square distribution
##'        "AnGamma": Anderson-Darling test statistic after map to an Erlang/Gamma distribution
##' @param useR logical indicating whether R or C implementations are used
##' @param ... additional arguments for the different methods
##' @return n-vector of values of the chosen test statistic
##' @author Marius Hofert and Martin Maechler
gofTstat <- function(u, method = c("Sn", "SnB", "SnC", "AnChisq", "AnGamma"),
		     useR = FALSE, ...)
{
    if(!is.matrix(u)) {
        warning("coercing 'u' to a matrix.")
        stopifnot(is.matrix(u <- as.matrix(u)))
    }
    d <- ncol(u)
    n <- nrow(u)
    method <- match.arg(method)
    switch(method,
           "Sn" =
       { ## S_n
           if(!hasArg(copula)) stop("object 'copula' required for pCopula() call")
           copula <- list(...)$copula
           ## compute \sum_{i=1}^n (\hat{C}_n(\bm{u}_i) - C_{\theta_n}(\bm{u}_i))^2
           ## typically \bm{u}_i ~> pobs \hat{\bm{U}}_i
           C. <- pCopula(u, copula=copula) # C_{\theta_n}(\bm{u}_i), i=1,..,n
           if(useR) {
               C.n <- F.n(u, X=u, ...) # \hat{C}_n(\bm{u}_i), i=1,..,n
               sum((C.n - C.)^2)
           } else {
               .C(cramer_vonMises,
                  as.integer(n),
                  as.integer(d),
                  as.double(u),
                  as.double(C.),
                  stat=double(1))$stat
           }
       },
	   "SnB" =
       { ## S_n(B)
           lu2 <- log1p(-u^2) # n x d matrix of log(1-u_{ij}^2)
           ## Note (modulo rowSums/colSums):
           ## Idea: sum1 = sum(prod(1-u^2)) = sum(exp(sum(lu2)))
           ## = exp(log( sum(exp(rowSums(lu2))) )) = exp(lsum(rowSums(lu2)))
           slu2 <- rowSums(lu2) # vector of length n
	   sum1 <- exp(lsum(cbind(slu2, deparse.level=0L))) # lsum() needs a matrix; result: 1 value
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
			       sum(pmin(lu.i, lu[,k])), NA_real_)
	       ls.i <- lsum(cbind(sum.k, deparse.level=0L)) # lsum( sum(pmin(...)) ) for fixed i; 1 value
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
           "AnChisq" = ad.test( pchisq(rowSums(qnorm(u)^2), df = d) )$statistic,
	   "AnGamma" = ad.test( pgamma(rowSums(-log(u)), shape = d) )$statistic,
	   stop("unsupported method ", method))
}


### The parametric bootstrap (for computing different goodness-of-fit tests) ###

## See gofCopula() for most arguments
gofPB <- function(copula, x, N, method = c("Sn", "SnB", "SnC"),
                  estim.method = c("mpl", "ml", "itau", "irho", "itau.mpl"),
		  trafo.method = if(method == "Sn") "none" else c("cCopula", "htrafo"),
		  trafoArgs = list(), test.method = c("family", "single"),
                  verbose = interactive(), useR = FALSE,
                  ties = NA, ties.method = c("max", "average", "first", "last", "random", "min"),
                  fit.ties.meth = eval(formals(rank)$ties.method), ...)
{
    .Deprecated("gofCopula(*, simulation = \"pb\")")
    .gofPB(copula, x, N, method=method, estim.method=estim.method,
           trafo.method=trafo.method, trafoArgs=trafoArgs, verbose=verbose,
           useR=useR, ties=ties, ties.method=ties.method,
           fit.ties.meth=fit.ties.meth,  ...)
}


## See gofCopula() for most arguments
.gofPB <- function(copula, x, N, method = c("Sn", "SnB", "SnC"),
                   estim.method = c("mpl", "ml", "itau", "irho", "itau.mpl"),
		   trafo.method = if(method == "Sn") "none" else c("cCopula", "htrafo"),
                   trafoArgs = list(), test.method = c("family", "single"),
                   verbose = interactive(), useR = FALSE,
                   ties = NA, ties.method = c("max", "average", "first", "last", "random", "min"),
                   fit.ties.meth = eval(formals(rank)$ties.method), ...)
{
    ## Checks -- NB: let the *generic* fitCopula() check 'copula'
    stopifnot(N >= 1)
    if(!is.matrix(x)) {
        warning("coercing 'x' to a matrix.")
        stopifnot(is.matrix(x <- as.matrix(x)))
    }
    stopifnot((d <- ncol(x)) > 1, (n <- nrow(x)) > 0, dim(copula) == d)
    method <- match.arg(method)
    estim.method <- match.arg(estim.method)
    test.method <- match.arg(test.method)
    if(method != "Sn")
	trafo.method <- match.arg(trafo.method, c("cCopula", "htrafo"))
    if(trafo.method == "htrafo") {
	if(!is(copula, "outer_nacopula"))
	    stop("'trafo.method' = \"htrafo\" only implemented for copula objects of type 'outer_nacopula'")
	if(length(copula@childCops))
	    stop("currently, only Archimedean copulas are supported")
    }

    ## Input checks
    if (method != "Sn" && trafo.method == "none")
        stop(sprintf("'trafo.method' must be \"cCopula\" or \"htrafo\" with 'method'=\"%s\"", method))
    if (method == "Sn" && trafo.method != "none")
        stop(sprintf("'trafo.method' must be \"none\" with 'method'=\"%s\"", method))

    ## Ties: by default, if at least one column has at least one duplicated entry
    if (is.na(ties <- as.logical(ties))) {
	ties <- any(apply(x, 2, anyDuplicated))
        if (ties)
            warning("argument 'ties' set to TRUE")
    }

    ## Progress bar
    if(verbose) {
	pb <- txtProgressBar(max = N+1, style=if(isatty(stdout())) 3 else 1) # setup progress bar
	on.exit(close(pb)) # on exit, close progress bar
    }

    ## 1) Compute the pseudo-observations
    uhat <- pobs(x, ties.method = ties.method)
    uhat.fit <- if (ties == FALSE || ties.method == fit.ties.meth) uhat
                else pobs(x, ties.method = fit.ties.meth)

    ## 2) Fit the copula
    ##    (if test.method = "family", otherwise take the provided copula; this
    ##     is useful for testing random number generators)
    C.th.n <- if(test.method == "family") {
                  fitter <- function(..., test.method)
                      fitCopula(copula, uhat.fit, method = estim.method,
                                estimate.variance = FALSE, ...)@copula
                  fitter(...) # avoids passing on 'test.method' to optim() [can be omitted if 'test.method' is a formal arg of gofCopula()]; see https://stackoverflow.com/questions/7028385/can-i-remove-an-element-in-dot-dot-dot-and-pass-it-on
              } else copula

    ## 3) Compute the realized test statistic
    doTrafo <- (method != "Sn" && trafo.method != "none") # (only) transform if method != "Sn" and trafo.method given
    u <- if(doTrafo) {
	stopifnot(is.list(trafoArgs))
	if(length(names(trafoArgs)) != length(trafoArgs))
	    stop("'trafoArgs' must be a fully named list")
	switch(trafo.method,
	       "cCopula"= do.call(cCopula, c(list(uhat, copula = C.th.n), trafoArgs)),
	       "htrafo" = do.call(htrafo,  c(list(uhat, copula = C.th.n), trafoArgs)),
	       stop("wrong transformation method"))
    } else uhat
    T <- if(method == "Sn") gofTstat(u, method = method, copula = C.th.n, useR = useR)
         else gofTstat(u, method = method)
    if(verbose) setTxtProgressBar(pb, 1) # update progress bar

    ## 4) Simulate the test statistic under H_0

    ## If ties, get tie structure from x  (FIXME?: what if 'x' has no ties, but 'ties = TRUE' ?? )
    ## IK: If 'x' has no ties, but 'ties = TRUE', extra computations for nothing
    if (ties)
        ir <- apply(x, 2, function(y) rank(sort(y)))

    T0 <- vapply(1:N, function(k) {

        ## 4.1) Sample the fitted (if test.method = "family") copula
        U <- rCopula(n, C.th.n)
        if(ties) { ## Sample x may have ties -- Reproduce tie structure of x
            for (i in 1:d) {
                U <- U[order(U[,i]),]
                U[,i] <- U[ir[,i], i]
            }
        }
        Uhat <- pobs(U, ties.method = ties.method)
        Uhat.fit <- if (ties == FALSE || ties.method == fit.ties.meth) Uhat
                    else pobs(U, ties.method = fit.ties.meth)

        ## 4.2) Fit the copula (if test.method = "family"; see Step 2))
        C.th.n. <- if(test.method == "family") {
                       fitter <- function(..., test.method)
                           fitCopula(copula, Uhat.fit, method = estim.method,
                                     estimate.variance = FALSE, ...)@copula
                       fitter(...) # see Step 2)
                   } else copula

        ## 4.3) Compute the test statistic
        u. <- if(doTrafo) { # (no checks needed; all done above)
	    switch(trafo.method,
		   "cCopula"= do.call(cCopula, c(list(Uhat, copula = C.th.n.), trafoArgs)),
		   "htrafo" = do.call(htrafo,  c(list(Uhat, copula = C.th.n.), trafoArgs)))
        } else Uhat
        T0. <- if(method=="Sn") gofTstat(u., method = method, copula = C.th.n., useR = useR)
               else gofTstat(u., method=method)

        if(verbose) setTxtProgressBar(pb, k+1) # update progress bar
        T0. # return
    }, NA_real_)

    ## 5) Return result object
    tr.string <- if (trafo.method == "none") ""
                 else sprintf(", 'trafo.method'=\"%s\"", trafo.method)
    structure(class = "htest",
	      list(method = paste0(.gofTstr("Parametric", copula, test.method),
				   sprintf(", with 'method'=\"%s\", 'estim.method'=\"%s\"%s:",
					   method, estim.method, tr.string)),
                   parameter = c(parameter = getTheta(C.th.n)),
                   statistic = c(statistic = T),
                   p.value = (sum(T0 >= T) + 0.5) / (N + 1), # typical correction => p-values in (0, 1)
                   data.name = deparse(substitute(x))))
}

## Auxiliary function for informative output
.gofTstr <- function(type, copula, test) {
    paste(type,
          "bootstrap-based goodness-of-fit test of",
          if(test == "single") "a single", # single: *the* exception;
          ## strongly recommended default "family": not mentioned (if only for back-compatib.)
          describeCop(copula, kind = "short"))
}


### The multiplier bootstrap (for computing different goodness-of-fit tests) ###

##' @title J-score \hat{J}_{\theta_n} for the multiplier bootstrap
##' @param copula An 'copula' (or 'parCopula' or ..)
##' @param u An (n, d)-matrix of (pseudo-)observations
##' @param method "mpl" or one of "itau", "irho"
##' @return A p by n matrix containing \hat{J}_{\theta_n}
##' @author Marius Hofert (based on ideas of Ivan Kojadinovic)
##' Note: References:
##'       * For estim.method="mpl":
##'         I. Kojadinovic and J. Yan (2011), A goodness-of-fit test for multivariate
##'         multiparameter copulas based on multiplier central limit theorems, Statistics
##'         and Computing 21:1, pages 17-30
##'       * For estim.method="itau" or ="irho" (d=2):
##'         I. Kojadinovic, J. Yan and M. Holmes (2011), Fast large-sample goodness-of-fit
##'         tests for copulas, Statistica Sinica 21:2, pages 841-871
##'       * for the function in general: van der Vaart "Asymptotic Statistics" (2000, p. 179 top)
Jscore <- function(copula, u, method)
{
    ## Checks
    stopifnot(is(copula, "Copula"))# and let methods below check more
    if(!is.matrix(u)) {
        warning("coercing 'u' to a matrix.")
        stopifnot(is.matrix(u <- as.matrix(u)))
    }
    stopifnot((n <- nrow(u)) > 0, (d <- ncol(u)) > 1, dim(copula) == d)

    ## Deal with different methods
    switch(method,
           "mpl"=
       {   ## See page 7 in Kojadinovic and Yan (2011)
           if(has.par.df(copula))
               stop("Jscore() cannot be computed for 'df.fixed = FALSE'")
           ## Integrals computed from n realizations by Monte Carlo
           influ0 <- dlogcdtheta(copula, u) # (n, p)-matrix
           derArg <- dlogcdu    (copula, u) # (n, d)-matrix
           influ <- lapply(1:d, function(i) influ0*derArg[,i])
           p <- nParam(copula, freeOnly=TRUE)
           S <- matrix(0, n, p)
           for(j in 1:d) {
               ij <- order(u[,j], decreasing=TRUE)
               ijb <- vapply(1:n, function(i) sum(t(u[,j]) <= u[i,j]), NA_real_)
	       S <- S +
		   rbind(rep.int(0, p),
			 apply(influ[[j]][ij,,drop=FALSE], 2L, cumsum))[n+1-ijb,,drop=FALSE] / n -
		   matrix(colMeans(influ[[j]]*u[,j]), n, p, byrow=TRUE)
           }
           Sigma.n <- crossprod(influ0) / n # = A^T A / n for A = influ0
           solve(Sigma.n, t(dlogcdtheta(copula, u) - S)) # solve(A, B) solves Ax = B => x = A^{-1}B
       },
           "ml"=
       {
           stop("method='ml' not available")
       },
           "itau"=
       { # See page 849 in Kojadinovic, Yan, and Holmes (2011)
           stopifnot(d == 2)
           (4/dTau(copula)) * ( 2*pCopula(u, copula) - rowSums(u) + (1-tau(copula))/2 )
       },
           "irho"=
       { # See Equation (3.5) on page 847 in Kojadinovic, Yan, and Holmes (2011)
           stopifnot(d == 2)
           i1 <- order(u[,1], decreasing=TRUE)
           i2 <- order(u[,2], decreasing=TRUE)
           n <- nrow(u)
           i1b <- vapply(1:n, function(k) sum(t(u[,1]) <= u[k,1]), NA_real_)
           i2b <- vapply(1:n, function(k) sum(t(u[,2]) <= u[k,2]), NA_real_)
           rmu <- rowMeans(u)
           term <- c(0, cumsum(u[,2][i1]))[n+1-i1b] / n - rmu +
               c(0, cumsum(u[,1][i2]))[n+1-i2b] / n - rmu
           (1/dRho(copula)) * ( 12*(apply(u, 1, prod) + term) - 3 - rho(copula) )
       },
           stop("wrong method"))
}

gofMB <- function(copula, x, N, method = c("Sn", "Rn"),
                  estim.method = c("mpl", "ml", "itau", "irho"),
                  test.method = c("family", "single"),
                  verbose = interactive(), useR = FALSE, m = 1/2,
                  zeta.m = 0, b = 1/sqrt(nrow(x)),
                  ties.method = c("max", "average", "first", "last", "random", "min"),
                  fit.ties.meth = eval(formals(rank)$ties.method), ...)
{
    .Deprecated("gofCopula(*, simulation = \"mult\")")
    .gofMB(copula, x, N, method=method, estim.method=estim.method,
           verbose=verbose, useR=useR, m=m, zeta.m=zeta.m, b=b,
           ties.method=ties.method, fit.ties.meth=fit.ties.meth,...)
}

##' revert to call  gofMB()  once that is gone
.gofMB <- function(copula, x, N, method = c("Sn", "Rn"),
                   estim.method = c("mpl", "ml", "itau", "irho"),
                   test.method = c("family", "single"),
                   verbose = interactive(), useR = FALSE, m = 1/2,
                   zeta.m = 0, b = 1/sqrt(nrow(x)),
                   ties.method = c("max", "average", "first", "last", "random", "min"),
                   fit.ties.meth = eval(formals(rank)$ties.method), ...)
{
    ## Checks -- NB: let the *generic* fitCopula() check 'copula'
    stopifnot(N >= 1)
    if(!is.matrix(x)) {
        warning("coercing 'x' to a matrix.")
        stopifnot(is.matrix(x <- as.matrix(x)))
    }
    stopifnot((d <- ncol(x)) > 1, (n <- nrow(x)) > 0, dim(copula) == d)
    method <- match.arg(method)
    estim.method <- match.arg(estim.method)
    test.method <- match.arg(test.method)
    if(estim.method == "ml") stop("estim.method='ml' not available")
    if(estim.method %in% c("irho", "itau") && d > 2)
        stop("only bivariate case possible for estim.method='irho' or ='itau'")

    ## 1) Compute the pseudo-observations
    u. <- pobs(x, ties.method = ties.method)
    u.fit <- if (ties.method == fit.ties.meth) u.
             else pobs(x, ties.method = fit.ties.meth)

    ## 2) Fit the copula
    ##    ('copula' is the H_0 copula; C.th.n = C_{\theta_n}(.)
    ##     unless 'test.method = "single"' in which case it is the
    ##     provided 'copula'; see Step 2) of .gofPB())
    C.th.n <- if(test.method == "family") {
                  fitter <- function(..., test.method)
                      fitCopula(copula, u.fit, method=estim.method, estimate.variance=FALSE, ...)@copula
                  fitter(...)
              } else copula

    ## 3) Compute the realized test statistic
    C.th.n. <- pCopula(u., C.th.n) # n-vector
    denom <-
        if (method == "Rn") # Rn
            (C.th.n.*(1-C.th.n.) + zeta.m)^m # n-vector
        else ##  "Sn"
            rep(1, n) # Sn
    Cn. <- F.n(u., u.) # n-vector
    T <- sum( ((Cn. - C.th.n.)/denom)^2 ) # test statistic Sn or Rn

    if (useR) { ## R version ################################################

        ## Obtain approximate realizations of the test statistic under H_0

        if (verbose)
            warning("The 'verbose' argument is ignored in 'gofMB' when 'useR' is set to TRUE")

        ## The multipliers
        Z <- matrix(rnorm(N*n), nrow=N, ncol=n) # (N, n)-matrix
        Zbar <- rowMeans(Z) # N-vector
        Zcent <- (Z-Zbar) # (N, n)-matrix

        ## The non-parametric part
        ## use a trick similar to C.n()
        ## note: - u. is an (n, d)-matrix
        ##       - t(u.)<=u.[i,] is an (d, n)-matrix
        ##       - colSums(t(u.)<=u.[i,])==d is an n-vector of logicals with kth entry
        ##         I_{\{\hat{\bm{U}}_k <= \bm{u}_i\}}
        ind <- vapply(1:n, function(i) colSums(t(u.) <= u.[i,]) == d, logical(n)) # (n, n)-matrix
        hCnh. <- Zcent %*% ind # (N, n)-matrix * (n, n)-matrix = (N, n)-matrix
        for(j in 1:d) { # bivariate only here
            ## build a matrix with 1s only, except for the jth column, which contains u.[,j]
            u.j <- matrix(1, nrow=n, ncol=d)
            u.j[,j] <- u.[,j]
            ind.u.j <- vapply(1:n, function(i) colSums(t(u.) <= u.j[i,]) == d, logical(n)) # (n, n)-matrix
            Cnh <- Zcent %*% ind.u.j # (N, n)-matrix * (n, n)-matrix = (N, n)-matrix
            ## finite-differences of the empirical copula at u.
            Cjn. <- dCn(u., U=u., j.ind=j, b=b) # n vector
            ## compute the sum
            Cjn.. <- matrix(rep(Cjn., each=N), nrow=N, ncol=n) # auxiliary (N, n)-matrix
            hCnh. <- hCnh. - Cjn.. * Cnh # (N, n)-matrix
        }

        ## Add the parametric part and compute multiplier replicates of T
        num <- (hCnh. - Z %*% t(dCdtheta(C.th.n, u=u.) %*%
                                Jscore(C.th.n, u=u., method=estim.method)))  # (N, n)-matrix
        T0 <- rowSums( (num %*% diag(1 / denom))^2 ) / n^2 # N vector; test statistic Sn^{(h)} or Rn^{(h)}
    } else {
        ## C version ################################################

        ## .C(cramer_vonMises_grid,
        ##    as.integer(d),
        ##    as.double(u.),
        ##    as.integer(n),
        ##    as.double(u.),
        ##    as.integer(n),
        ##    as.double(pCopula(u., C.th.n)), # C_{\theta_n}(\hat{\bm{U}}_i), i=1,..,n
        ##    stat=double(1))$stat
        ##   gofTstat(u., method="Sn", copula=C.th.n)

        ## 4) Simulate the test statistic under H_0
        T0 <- .C(multiplier,## -> ../src/gof.c
                 p = as.integer(d),
                 U = as.double(u.),
                 n = as.integer(n),
                 G = as.double(u.),
                 g = as.integer(n),
                 b = as.double(b),
                 influ = as.double(dCdtheta(C.th.n, u.) %*%
                                   Jscore(C.th.n, u=u., method=estim.method)),
                 denom = as.double(denom),
                 N = as.integer(N),
                 s0 = double(N),
                 verbose = as.integer(verbose))$s0
    }

    ## Return result object
    structure(class = "htest",
	      list(method = paste0(.gofTstr("Multiplier", copula, test.method),
				   sprintf(", with 'method'=\"%s\", 'estim.method'=\"%s\":",
					  method, estim.method)),
                   parameter = c(parameter = getTheta(C.th.n)),
                   statistic = c(statistic = T),
                   p.value = (sum(T0 >= T) + 0.5) / (N + 1), # typical correction =>  p-values in (0, 1)
                   data.name = deparse(substitute(x))))
}


### Wrapper ####################################################################

##' @title Goodness-of-fit test wrapper function
##' @param copula An object of type 'copula' representing the H_0 copula
##' @param x An (n, d)-matrix containing the data
##' @param N The number of bootstrap (parametric or multiplier) replications
##' @param method A goodness-of-fit test statistic to be used
##' @param estim.method An estimation method for the unknown parameter vector
##' @param simulation The parametric ('pb') or multiplier ('mult') bootstrap
##' @param verbose A logical indicating whether a progress bar is shown
##' @param ties.method passed to pobs (not for fitting)
##' @param fit.ties.meth passed to pobs (fitting only)
##' @param ... Additional arguments passed to the internal auxiliary functions
##'        gofPB() and gofMB()
##' @return An object of class 'htest'
##' @author Ivan Kojadinovic, Marius Hofert
gofCopulaCopula <- function(copula, x, N=1000, method = c("Sn", "SnB", "SnC", "Rn"),
                            estim.method = c("mpl", "ml", "itau", "irho", "itau.mpl"),
                            simulation = c("pb", "mult"), test.method = c("family", "single"),
                            verbose = interactive(),
                            ties = NA, ties.method = c("max", "average", "first", "last", "random", "min"),
                            fit.ties.meth = eval(formals(rank)$ties.method), ...)
{
    ## Checks
    stopifnot(N >= 1)
    if(!is.matrix(x)) {
        warning("coercing 'x' to a matrix.")
        stopifnot(is.matrix(x <- as.matrix(x)))
    }
    stopifnot((d <- ncol(x)) > 1, nrow(x) > 0, dim(copula) == d)
    method <- match.arg(method)
    estim.method <- match.arg(estim.method)
    simulation <- match.arg(simulation)
    test.method <- match.arg(test.method)
    ties.method <- match.arg(ties.method)
    fit.ties.meth <- match.arg(fit.ties.meth)

    if (method == "Sn" && is(copula, "tCopula") && !copula@df.fixed)
	stop("'method'=\"Sn\" not available for t copulas whose df are not fixed as pCopula() cannot be computed for non-integer degrees of freedom yet.")

    ## Deprecation
    ## if (!is.null(print.every)) {
    ##     warning("Argument 'print.every' is deprecated. Please use 'verbose' instead.")
    ##     verbose <- print.every > 0
    ## }
    ## Back-compatibility
    if(missing(estim.method) && !missing(method)) {
        eMeth <- eval(formals()$estim.method)
	if(!is.na(i <- pmatch(method, eMeth))) {
	    warning("old (pre 0.999-*) argument 'method' is now called 'estim.method'")
	    estim.method <- eMeth[i]
	    method <- "Sn"
	}
    }

    ## Distinguish the methods
    switch(simulation,
           "pb" = { ## parametric bootstrap
               .gofPB(copula, x, N = N, method = method, estim.method = estim.method,
                      test.method = test.method, verbose = verbose,
                      ties = ties, ties.method = ties.method,
                      fit.ties.meth = fit.ties.meth, ...)
           },
           "mult" = { ## multiplier bootstrap
               .gofMB(copula, x = x, N = N, method = method, estim.method = estim.method,
                      test.method = test.method, verbose = verbose,
                      ties.method = ties.method, fit.ties.meth = fit.ties.meth, ...)
           },
           stop("Invalid simulation method ", simulation))
}

setMethod("gofCopula", signature("parCopula"), gofCopulaCopula)

