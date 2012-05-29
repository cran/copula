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


## Computes the sum of binomial coefficients
binom.sum <- function(n,k)
{
    bs <- 1
    for (i in 1:k)
        bs <- bs + choose(n,i)
    bs
}

################################################################################

##  Multivariate independence test based on the empirical
##  copula process as proposed by Christian Genest and Bruno
##  Rémillard (2004), Test 13:2, pages 335-369.
##
##  Ivan Kojadinovic, May 2007

# simulate the distribution of TA for A of cardinality 2 to m under independence

indepTestSim <- function(n,p,m=p,N=1000,print.every=100)
{
    if (!is.numeric(n) || (n <- as.integer(n)) < 2)
        stop("n should be an integer greater than 2")
    if (!is.numeric(p) || (p <- as.integer(p)) < 2)
        stop("p should be an integer greater than 2")
    if (!is.numeric(m) || (m <- as.integer(m)) < 2 || m > p)
        stop(paste("m should be an integer greater than 2 and smaller than",p))
    if (!is.numeric(N) || as.integer(N) < 100)
        stop("N should be an integer greater than 100")

    sb <- binom.sum(p,m)

    r. <- .C(simulate_empirical_copula,
             n = as.integer(n),
             N = as.integer(N),
             p = as.integer(p),
             m = as.integer(m),
             TA0 = double(N * (sb - p - 1)),
             G0 = double(N),
             subsets = integer(sb),
             subsets.char = character(sb),
             fisher0 = double(N),
             tippett0 = double(N),
             as.integer(print.every))

    structure(class = "indepTestDist",
              list(sample.size = n,
                   data.dimension = p,
                   max.card.subsets = m,
                   number.repetitions = N,
                   subsets = r.$subsets.char[(p+2):sb],
                   subsets.binary = r.$subsets[(p+2):sb],
                   dist.statistics.independence = matrix(r.$TA0,N,sb-p-1),
                   dist.global.statistic.independence = r.$G0,
                   dist.fisher.independence = r.$fisher0,
                   dist.tippett.independence = r.$tippett0))
}

################################################################################

# test in itself

indepTest <- function(x, d, alpha=0.05)
{
    if (!is.numeric(x <- as.matrix(x)))
        stop("data should be numerical")
    if (any(is.na(x)))
        stop("data cannot contain missing values")
    if (nrow(x) < 2)
        stop("data should contain more than 2 rows")
    if (class(d) != "indepTestDist")
        stop("d should be obtained by means of the function empcopu.simulate")
    if (ncol(x) != d$data.dimension)
      stop("d was not obtained from simulations based on data whose dimension is equal to that of x")
    if (nrow(x) != d$sample.size)
      warning("d was not obtained from simulations based on data whose size is equal to that of x")
    if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1)
        stop("the significance level alpha is not properly set")

    m <- d$max.card.subsets
    p <- ncol(x)
    N <- as.integer(d$number.repetitions)

    ## transform data to ranks
    for (j in 1:p)
        x[,j] <- rank(x[,j])

    sb <- binom.sum(p,m)

    ## perform test
    r. <- .C(empirical_copula_test,
             as.double(x),
             nrow(x),
             as.integer(p),
             as.integer(m),
             as.double(d$dist.statistics.independence),
             as.double(d$dist.global.statistic.independence),
             as.integer(N),
             as.integer(d$subsets.binary),
             TA = double(sb - p - 1),
             G = double(1),
             pval = double(sb - p - 1),
             fisher = double(1),
             tippett = double(1),
             globpval = double(1),
             as.double(d$dist.fisher.independence),
             as.double(d$dist.tippett.independence))

    ## compute critical values at the alpha level
    beta <- (1 - alpha)^(1 / (sb - p -1))
    critical <- numeric(0)
    from <- 1
    for (j in 2:m)
    {
        to <- from + choose(p,j) - 1
        crit <- sort(as.double(d$dist.statistics.independence[,from:to]))[round(beta * N * (to - from + 1))]
        critical <- c(critical,rep(crit,choose(p,j)))
        from <- to + 1
    }

    structure(class = "indepTest",
              list(subsets=d$subsets,statistics=r.$TA, critical.values=critical,
                   pvalues = r.$pval, fisher.pvalue=r.$fisher, tippett.pvalue=r.$tippett, alpha=alpha,
                   beta=beta, global.statistic=r.$G, global.statistic.pvalue=r.$globpval))
}

################################################################################

##  Serial independence test based on the empirical
##  copula process as proposed by Christian Genest and Bruno
##  Rémillard (2004), Test 13:2, pages 335-369.
##
##  Ivan Kojadinovic, May 2007

## simulate the distribution of TAs for A of cardinality 2 to m
## containing 1 under serial independence

serialIndepTestSim <- function(n,lag.max,m=lag.max+1,N=1000,print.every=100)
{
    if (!is.numeric(n) || (n <- as.integer(n)) < 2)
        stop("n should be an integer greater than 2")
    if (!is.numeric(lag.max) || (lag.max <- as.integer(lag.max)) < 1)
        stop("lag.max should be an integer greater than 1")

    p <- lag.max + 1 ## dimension of the virtual data

    if (n-p+1 < 2)
      stop("wrong number of lags with respect to the sample size")

    if (!is.numeric(m) || (m <- as.integer(m)) < 2 || m > p)
        stop(paste("m should be an integer greater than 2 and smaller than",p))
    if (!is.numeric(N) || (N <- as.integer(N)) < 100)
        stop("N should be an integer greater than 100")
    if (!is.numeric(lag.max) || (p <- as.integer(lag.max) + 1) <= 1 || n-p+1 < 2)
      stop("wrong number of lags")

    sb <- binom.sum(p-1,m-1)

    R <- .C(simulate_empirical_copula_serial,
             n = as.integer(n-p+1),
             N = as.integer(N),
             p = as.integer(p),
             m = as.integer(m),
             TA0 = double(N * (sb - 1)),
             G0 = double(N),
             subsets = integer(sb),
             subsets.char = character(sb),
             fisher0 = double(N),
             tippett0 = double(N),
             as.integer(print.every))

    structure(class = "serialIndepTestDist",
	      list(sample.size = n,
		   lag.max = lag.max,
		   max.card.subsets = m,
		   number.repetitions = N,
		   subsets = R $subsets.char[2:sb],
		   subsets.binary = R $subsets[2:sb],
		   dist.statistics.independence = matrix(R $TA0,N,sb-1),
		   dist.global.statistic.independence = R $G0,
		   dist.fisher.independence = R $fisher0,
		   dist.tippett.independence = R $tippett0))
}

################################################################################

# test in itself

serialIndepTest <- function(x, d, alpha=0.05)
{
    if (!is.numeric(x))
      stop("data should be numerical")
    x <- as.numeric(x)
    if (any(is.na(x)))
      stop("data cannot contain missing values")
    if (length(x) < 2)
      stop("data should contain more than 2 observations")
    if (class(d) != "serialIndepTestDist")
      stop("d should be obtained by means of the function empcops.simulate")
    n <- length(x)
    if (n != as.integer(d$sample.size))
      warning("d was not obtained from simulations based on data whose size is equal to that of x")
    if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1)
      stop("the significance level alpha is not properly set")

    ## transform data to pseudo-obs
    x <- rank(x)/n

    m <- d$max.card.subsets
    p <- as.integer(d$lag.max) + 1
    n <- n - p + 1
    N <- as.integer(d$number.repetitions)

    sb <- binom.sum(p-1,m-1)

    ## perform test
    R <- .C(empirical_copula_test_serial,
            as.double(x),
            as.integer(n),
            as.integer(p),
            as.integer(m),
            as.double(d$dist.statistics.independence),
            as.double(d$dist.global.statistic.independence),
            as.integer(N),
            as.integer(d$subsets.binary),
            TA = double(sb - 1),
            G = double(1),
            pval = double(sb - 1),
            fisher = double(1),
            tippett = double(1),
            globpval = double(1),
            as.double(d$dist.fisher.independence),
            as.double(d$dist.tippett.independence))

    ## compute critical values at the alpha level
    beta <- (1 - alpha)^(1 / (sb -1))
    critical <- numeric(0)
    from <- 1
    for (j in 1:(m-1))
    {
        to <- from + choose(p-1,j) - 1
        crit <- sort(as.double(d$dist.statistics.independence[,from:to]))[round(beta * N * (to - from + 1))]
        critical <- c(critical,rep(crit,choose(p-1,j)))
        from <- to + 1
    }

    structure(class = "indepTest",
	      list(subsets=d$subsets,statistics=R $TA, critical.values=critical,
		   pvalues = R $pval, fisher.pvalue=R $fisher, tippett.pvalue=R $tippett, alpha=alpha,
		   beta=beta, global.statistic=R $G, global.statistic.pvalue=R $globpval))
}

################################################################################

##  Independence test among random vectors based on the empirical
##  copula process
##
##  Ivan Kojadinovic, December 2007

multIndepTest <- function(x, d, m=length(d), N=1000, alpha=0.05,
                          print.every=100)
{
    if (!is.numeric(x <- as.matrix(x)))
        stop("data should be numerical")
    if (any(is.na(x)))
        stop("data cannot contain missing values")
    if ((n <- nrow(x)) < 2)
        stop("data should contain more than 2 rows")
    if (!is.numeric(d) || any((d <- as.integer(d)) < 1) || sum(d) != (nc <- ncol(x)))
      stop(paste("wrong vector of dimensions"))

    p <- length(d)

    if (!is.numeric(m) || (m <- as.integer(m)) < 2 || m > p)
        stop(paste("m should be an integer greater than 2 and smaller than",p))
    if (!is.numeric(N) || (N <- as.integer(N)) < 100)
        stop("N should be an integer greater than 100")
    if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1)
        stop("the significance level alpha is not properly set")

    ## transform data to pseudo-observations
    for (j in 1:nc)
        x[,j] <- rank(x[,j])/n

    sb <- binom.sum(p,m)

    ## block indices
    b <- c(0,cumsum(d))

    ## bootstrap
    bootstrap <- .C(bootstrap_MA_I,
                    n = as.integer(n),
                    N = as.integer(N),
                    p = as.integer(p),
                    b = as.integer(b),
                    U = as.double(x),
                    m = as.integer(m),
                    MA0 = double(N * (sb - p - 1)),
                    I0 = double(N),
                    subsets = integer(sb),
                    subsets.char = character(sb),
                    as.integer(print.every))

    subsets <- bootstrap$subsets.char[(p+2):sb]
    subsets.binary <- bootstrap$subsets[(p+2):sb]
    dist.statistics.independence <- matrix(bootstrap$MA0,N,sb-p-1)

    ## perform test
    R <- .C(empirical_copula_test_rv,
             as.double(x),
             as.integer(n),
             as.integer(p),
             as.integer(b),
             as.integer(m),
             as.double(bootstrap$MA0),
             as.double(bootstrap$I0),
             as.integer(N),
             as.integer(subsets.binary),
             MA = double(sb - p - 1),
             I = double(1),
             pval = double(sb - p - 1),
             fisher = double(1),
             tippett = double(1),
             Ipval = double(1)
	     )## FIXME: additional argument  as.integer(print.every)

    ## compute critical values at the alpha level
    beta <- (1 - alpha)^(1 / (sb - p - 1))
    critical <- numeric(sb - p - 1)

    for (k in 1:(sb - p - 1))
      critical[k] <- sort(dist.statistics.independence[,k])[round(beta * N)]

     structure(class = "indepTest",
	       list(subsets=subsets,statistics=R $MA, critical.values=critical,
		    pvalues = R $pval, fisher.pvalue=R $fisher, tippett.pvalue=R $tippett,
		    alpha=alpha, beta=beta,
		    global.statistic=R $I, global.statistic.pvalue=R $Ipval))
}

################################################################################

##  Multivariate serial independence test based on the empirical
##  copula process
##
##  Ivan Kojadinovic, December 2007

multSerialIndepTest <- function(x, lag.max, m=lag.max+1, N=1000, alpha=0.05,
                                print.every=100)
{
    if (!is.numeric(x <- as.matrix(x)))
        stop("data should be numerical")
    if (any(is.na(x)))
        stop("data cannot contain missing values")
    if ((n <- nrow(x)) < 2)
        stop("data should contain more than 2 rows")
    if (!is.numeric(lag.max) || (p <- as.integer(lag.max) + 1) <= 1 || n-p+1 < 2)
      stop("wrong number of lags")
    if (!is.numeric(m) || ((m <- as.integer(m)) < 2) || m > p)
        stop(paste("m should be an integer greater than 2 and smaller than",p))
    if (!is.numeric(N) || (N <- as.integer(N)) < 100)
        stop("N should be an integer greater than 100")
    if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1)
        stop("the significance level alpha is not properly set")

    q <- ncol(x)

    ## transform data to ranks
    for (j in 1:q)
        x[,j] <- rank(x[,j])/n

    n <- n-p+1 ## number of rows of p-dimensional virtual data

    sb <- binom.sum(p-1,m-1)

    ## bootstrap
    bootstrap <- .C(bootstrap_serial,
                    n = as.integer(n),
                    N = as.integer(N),
                    p = as.integer(p),
                    q = as.integer(q),
                    U = as.double(x),
                    m = as.integer(m),
                    MA0 = double(N * (sb - 1)),
                    I0 = double(N),
                    subsets = integer(sb),
                    subsets.char = character(sb),
                    as.integer(print.every))

    subsets <- bootstrap$subsets.char[2:sb]
    subsets.binary <- bootstrap$subsets[2:sb]
    dist.statistics.independence <- matrix(bootstrap$MA0,N,sb-1)

    ## perform test
    R <- .C(empirical_copula_test_rv_serial,
             as.double(x),
             as.integer(n),
             as.integer(p),
             as.integer(q),
             as.integer(m),
             as.double(bootstrap$MA0),
             as.double(bootstrap$I0),
             as.integer(N),
             as.integer(subsets.binary),
             MA = double(sb - 1),
             I = double(1),
             pval = double(sb - 1),
             fisher = double(1),
             tippett = double(1),
             Ipval = double(1))

    ## compute critical values at the alpha level
    beta <- (1 - alpha)^(1 / (sb - 1))
    critical <- numeric(0)
    from <- 1
    for (j in 1:(m-1))
    {
        to <- from + choose(p-1,j) - 1
        crit <- sort(as.double(dist.statistics.independence[,from:to]))[round(beta * N * (to - from + 1))]
        critical <- c(critical,rep(crit,choose(p-1,j)))
        from <- to + 1
    }

    structure(class = "indepTest",
	      list(subsets=subsets,statistics=R $MA, critical.values=critical,
		   pvalues = R $pval, fisher.pvalue=R $fisher, tippett.pvalue=R $tippett, alpha=alpha,
		   beta=beta, global.statistic=R $I, global.statistic.pvalue=R $Ipval))
}

## Dependogram #################################################################

dependogram <- function(test, pvalues=FALSE, print=FALSE)
{
  if (!inherits(test, "indepTest"))
    stop("'test' should be obtained by means of the functions indepTest, multIndepTest, serialIndepTest, or multSerialIndepTest")
  stopifnot(is.logical(pvalues), is.logical(print))

  op <- par(las=3); on.exit(par(op))

  l <- length(test$statistics)
  if (pvalues)
    {
      plot(c(1,1:l),c(0,test$pvalues), type="h", ylab="p-value per subset", xlab="",
           axes=FALSE, main="Dependogram")
      axis(1, at=1:l, labels=test$subsets)
      axis(2)
      points(1:l,rep(1 - test$beta, l),pch=20)
    }
  else ## not pvalues
    {
      plot(c(1,1:l),c(0,test$statistics), type="h", ylab="statistic per subset", xlab="",
           axes=FALSE, main="Dependogram")
      axis(1, at=1:l, labels=test$subsets)
      axis(2)
      points(1:l,test$critical.values,pch=20)
    }
  if (print)
    {
      cat("The subset statistics, p-values and critical values are:\n")
      print(data.frame(subset=test$subsets, statistic=test$statistics,
                       pvalue=test$pvalues, critvalue=test$critical.values))
      cat("The critical values are such that the simultaneous acceptance region \n",
          "has probability 1 -", test$alpha, "under the null.\n")
      cat("The individual rejection probability for any statistic obtained from the Mobius \n",
          "decomposition is 1 -", test$beta, "under the null.\n\n")
    }
}


## summary and show methods for class 'indepTest' ##############################

print.indepTest <- function(x, ...)
{
  cat("\nGlobal Cramer-von Mises statistic:", x$global.statistic,
      "with p-value", x$global.statistic.pvalue, "\n")
  cat("Combined p-values from the Mobius decomposition:\n")
  cat("  ", x$fisher.pvalue, " from Fisher's rule,\n")
  cat("  ", x$tippett.pvalue, " from Tippett's rule.\n")

}

################################################################################

## myEcdf <- function(x) {
##   N <- length(x)
##   fun <- ecdf(x)
##   myfun <- function(x) (0.5 + N * fun(x)) / (N + 1)
##   myfun
## }

## myEcdf.caglad <- function(x) {
##   N <- length(x)
##   fun <- ecdf(-x)
##   myfun <- function(x) (0.5 + N * fun(-x)) / (N + 1)
##   myfun
## }

## combinePvals <- function(MAn) {
##   MAn <- t(MAn)
##   nA <- ncol(MAn)
##   N <- nrow(MAn) - 1
##   ## the observed is in the 1st row
## ##   OneMinusFhat <- apply(MAn, 2,
## ##                         function(x) {
## ##                           myfun <- myEcdf(x[-1])
## ##                           1 - myfun(x)
## ##                         })
## ##   Fisher <- -2 * apply(OneMinusFhat, 1, function(x) sum(log(x)))
## ##   Tippett <-  apply(OneMinusFhat, 1, min)
##   pvalues <- apply(MAn, 2,
##                    function(x) {
##                      myfun <- myEcdf.caglad(x[-1])
##                      myfun(x)
##                    })
##   Fisher <- -2 * apply(pvalues, 1, function(x) sum(log(x)))
##   Tippett <-  apply(pvalues, 1, min)
##   Fisher.pvalue <- (0.5 + sum(Fisher[-1] >= Fisher[1])) / (N + 1)
##   Tippett.pvalue <- (0.5 + sum(Tippett[-1] <= Tippett[1])) / (N + 1)
##   list(Fisher.pvalue=Fisher.pvalue, Tippett.pvalue=Tippett.pvalue)
## }

## myCritical <- function(dist.statistics.independence, beta, p, m, N) {
##   ## compute critical values at the alpha level
##   critical <- numeric(0)
##   from <- 1
##   for (j in 1:(m-1)) {
##     to <- from + choose(p-1,j) - 1
##     crit <- sort(as.double(dist.statistics.independence[,from:to]))[round(beta * N * (to - from + 1))]
##     critical <- c(critical,rep(crit,choose(p-1,j)))
##     from <- to + 1
##   }
##   critical
## }

## empcopsm.test.jy <- function(x, lag.max, m=lag.max + 1, N=1000, alpha=0.05) {
##   if (!is.numeric(x <- as.matrix(x)))
##     stop("data should be numerical")
##   if (any(is.na(x)))
##     stop("data cannot contain missing values")
##   if ((n <- nrow(x)) < 2)
##     stop("data should contain more than 2 rows")
##   if (!is.numeric(lag.max) || (p <- as.integer(lag.max) + 1) <= 1 || n-p+1 < 2)
##     stop("wrong number of lags")

##   q <- ncol(x)
##   n <- n - p + 1 ## number of rows of p-dimensional virtual data

##   if (!is.numeric(m) || ((m <- as.integer(m)) < 2) || m > p)
##     stop(paste("m should be an integer greater than 2 and smaller than",p))
##   if (!is.numeric(N) || (N <- as.integer(N)) < 100)
##     stop("N should be an integer greater than 100")
##   if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1)
##     stop("the significance level alpha is not properly set")

##   ## transform data to ranks
##   for (j in 1:q) x[,j] <- rank(x[,j])

##   sb <- binom.sum(p-1,m-1)

##   ## power set operation
##   powerSet <- 2 ^ as.set(0 : (p - 1))
##   Asize <- unlist(lapply(powerSet, length))
##   hasZero <- unlist(lapply(powerSet, function(x) 0 %e% x))
##   good <- Asize > 1 & Asize <= m & hasZero
##   Asize <- Asize[good]
##   Aidx <- unlist(powerSet[good])
##   nA <- length(Asize)
##   subsets <- powerSet[good]

##   run <- .C("empcopsm", as.integer(x), as.integer(N),
##             as.integer(n), as.integer(p), as.integer(q),
##             as.integer(Aidx), as.integer(Asize), as.integer(nA),
##             In = double(N + 1), MAn = double((N + 1) * nA),
##             package = "copula")
##   MAn <- matrix(run$MAn, nrow=nA)
##   MAn.obs <- MAn[,1]
##   MAn.sim <- MAn[,-1]
##   MAn.pvalues <- apply(MAn, 1, function(x) (sum(x[-1] >= x[1]) + 0.5) / (N + 1) )
##   In.obs <- run$In[1]
##   In.sim <- run$In[-1]
##   In.pvalue <- (sum(In.obs <= In.sim) + 0.5) / (N + 1)

##   MAn.sim <- t(MAn.sim)
##   pcombined <- combinePvals(MAn)

##   beta <- (1 - alpha)^(1 / (sb - 1))
##   critical <- myCritical(MAn.sim, beta, p, m, N)
##   test <- list(subsets = subsets,
##                statistics = MAn.obs,
##                critical.values=critical,
##                pvalues = MAn.pvalues,
##                fisher.pvalue = pcombined$Fisher.pvalue,
##                tippett.pvalue = pcombined$Tippett.pvalue,
##                alpha = alpha, beta = beta,
##                global.statistic = In.obs, global.statistic.pvalue = In.pvalue)
##   class(test) <- "empcop.test"
##   return(test)
## }

## empcopsm.test.ik <- function(x, lag.max, m=lag.max+1, N=1000, alpha=0.05)
## {
##     if (!is.numeric(x <- as.matrix(x)))
##         stop("data should be numerical")
##     if (any(is.na(x)))
##         stop("data cannot contain missing values")
##     if ((n <- nrow(x)) < 2)
##         stop("data should contain more than 2 rows")
##     if (!is.numeric(lag.max) || (p <- as.integer(lag.max) + 1) <= 1 || n-p+1 < 2)
##       stop("wrong number of lags")

##     q <- ncol(x)
##     n <- n-p+1 ## number of rows of p-dimensional virtual data

##     if (!is.numeric(m) || ((m <- as.integer(m)) < 2) || m > p)
##         stop(paste("m should be an integer greater than 2 and smaller than",p))
## ##     if (!is.numeric(N) || (N <- as.integer(N)) < 100)
## ##         stop("N should be an integer greater than 100")
##     if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1)
##         stop("the significance level alpha is not properly set")

##     ## transform data to ranks
##     for (j in 1:q)
##         x[,j] <- rank(x[,j])

##     sb <- binom.sum(p-1,m-1)

##     ## bootstrap
##     bootstrap <- .C("bootstrap_serial",
##                     n = as.integer(n),
##                     N = as.integer(N),
##                     p = as.integer(p),
##                     q = as.integer(q),
##                     R = as.integer(x),
##                     m = as.integer(m),
##                     MA0 = double(N * (sb - 1)),
##                     I0 = double(N),
##                     subsets = integer(sb),
##                     subsets.char = character(sb))#,                    PACKAGE="copula")

##     subsets <- bootstrap$subsets.char[2:sb]
##     subsets.binary <- bootstrap$subsets[2:sb]
##     dist.statistics.independence <- matrix(bootstrap$MA0,N,sb-1)

##     ## perform test
##     R <- .C("empirical_copula_test_rv_serial",
##              as.integer(x),
##              as.integer(n),
##              as.integer(p),
##              as.integer(q),
##              as.integer(m),
##              as.double(bootstrap$MA0),
##              as.double(bootstrap$I0),
##              as.integer(N),
##              as.integer(subsets.binary),
##              MA = double(sb - 1),
##              I = double(1),
##              pval = double(sb - 1),
##              fisher = double(1),
##              tippett = double(1),
##              Ipval = double(1)) #,             PACKAGE="copula")

##     ## compute critical values at the alpha level
##     beta <- (1 - alpha)^(1 / (sb - 1))
##     critical <- numeric(0)
##     from <- 1
##     for (j in 1:(m-1))
##     {
##         to <- from + choose(p-1,j) - 1
##         crit <- sort(as.double(dist.statistics.independence[,from:to]))[round(beta * N * (to - from + 1))]
##         critical <- c(critical,rep(crit,choose(p-1,j)))
##         from <- to + 1
##     }

##     test <- list(subsets=subsets,statistics=R $MA, critical.values=critical,
##                  pvalues = R $pval, fisher.pvalue=R $fisher, tippett.pvalue=R $tippett, alpha=alpha,
##                  beta=beta, global.statistic=R $I, global.statistic.pvalue=R $Ipval,
##                  MAn.sim=dist.statistics.independence)

##     class(test) <- "empcop.test"

##     return(test)
## }

## ## Computes the sum of binomial coefficients
## binom.sum <- function(n,k)
## {
##     bs <- 1
##     for (i in 1:k)
##         bs <- bs + choose(n,i)
##     bs
## }
