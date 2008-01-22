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
##     res<- .C("empirical_copula_test_rv_serial",
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

##     test <- list(subsets=subsets,statistics=res$MA, critical.values=critical,
##                  pvalues = res$pval, fisher.pvalue=res$fisher, tippett.pvalue=res$tippett, alpha=alpha,
##                  beta=beta, global.statistic=res$I, global.statistic.pvalue=res$Ipval,
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
