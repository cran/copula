#################################################################
##   Copula R package by Jun Yan Copyright (C) 2007
##
##   Farlie-Gumbel-Morgenstern multivariate copula class 
##   Copyright (C) 2007 Ivan Kojadinovic <ivan.kojadinovic@univ-nantes.fr>
##
##   This program is free software; you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation; either version 2 of the License, or
##   (at your option) any later version.
##
##   This program is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
##   You should have received a copy of the GNU General Public License along
##   with this program; if not, write to the Free Software Foundation, Inc.,
##   51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
##
#################################################################
##
##  Multivariate independence test based on the empirical 
##  copula process as proposed by Christian Genest and Bruno 
##  Rémillard (2004), Test 13:2, pages 335-369.

##  Ivan Kojadinovic, May 2007

##############################################################################

## Computes the sum of binomial coefficients
binom.sum <- function(n,k)
{
    bs <- 1
    for (i in 1:k)
        bs <- bs + choose(n,i)
    bs
}

##############################################################################

# simulate the distribution of TA for A of cardinality 2 to m under independence

empcop.simulate <- function(n,p,m,N=2000)
{
    if (!is.numeric(n) || as.integer(n) < 2)
        stop("n should be an integer greater than 2")
    if (!is.numeric(p) || as.integer(p) < 2)
        stop("p should be an integer greater than 2")
    if (!is.numeric(m) || as.integer(m) < 2 || as.integer(m) > as.integer(p))
        stop("m should be an integer greater than 2 and smaller than p")
    if (!is.numeric(N) || as.integer(N) < 100)
        stop("N should be an integer greater than 100")

    sb <- binom.sum(p,m)

    cat(paste("The simulation can take a long time as",sb-p-1,"(number of subsets) times",N,"(number of repetitions) statistics are to be computed.\n"))
  
    res<- .C("simulate_empirical_copula",
             n = as.integer(n),
             N = as.integer(N),
             p = as.integer(p),
             m = as.integer(m),
             TA0 = double(N * (sb - p - 1)),
             subsets = integer(sb),
             subsets.char = character(sb),
             PACKAGE="copula")

    d <- list(sample.size = n,
              data.dimension = p,
              max.card.subsets = m,
              number.repetitions = N,
              subsets = res$subsets.char[(p+2):sb],
              subsets.binary = res$subsets[(p+2):sb],
              dist.statistics.independence = matrix(res$TA0,N,sb-p-1))
    
    class(d) <- "empcop.distribution"
           
    return(d)
}

##############################################################################

empcop.test <- function(x, d, m, alpha=0.05)
{
    if (!is.numeric(x <- as.matrix(x)))
        stop("data should be numerical")
    if (sum(is.na(x)) > 0)
        stop("data cannot contain missing values")
    if (nrow(x) < 2)
        stop("data should contain more than 2 rows")
    if (class(d) != "empcop.distribution")
        stop("d should be obtained by means of the function empcop.simulate")
    if (ncol(x) != d$data.dimension)
        stop("d was not obtained from simulations based on data whose dimension is equal to that of x")
    if (!is.numeric(m) || as.integer(m) < 2 || as.integer(m) > as.integer(d$max.card.subsets))
        stop(paste("m should be an integer greater than 2 and smaller than",d$max.card.subsets))
    if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1)
        stop("the significance level alpha is not properly set")

    p <- ncol(x)
    N <- d$number.repetitions
      
    ## transform data to ranks
    for (j in 1:p)
        x[,j] <- rank(x[,j])

    sb <- binom.sum(p,m)

    ## perform test
    res<- .C("empirical_copula_test",
             as.integer(x),
             nrow(x),
             as.integer(p),
             as.integer(m),
             as.double(d$dist.statistics.independence),
             as.integer(N),
             as.integer(d$subsets.binary),
             TA = double(sb - p - 1),
             pval = double(sb - p - 1),
             fisher = double(1),
             tippett = double(1),
             PACKAGE="copula")
    
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
    
    test <- list(subsets=d$subsets[1:(sb - p - 1)],statistics=res$TA, critical.values=critical,
                 pvalues = res$pval, fisher.pvalue=res$fisher, tippett.pvalue=res$tippett, alpha=alpha, beta=beta)

    class(test) <- "empcop.test"
    return(test)
}

##############################################################################

dependogram <- function(test, pvalues=FALSE)
{
    if (class(test) != "empcop.test")
        stop("'test' should be obtained by means of the function empcop.test")
    if (!is.logical(pvalues))
        stop("'pvalues' should be a boolean")

    l <- length(test$statistics)
    if (pvalues == TRUE)
    {
        plot(c(1,1:l),c(0,test$pvalues), type="h", ylab="p-value per subset", xlab="subsets",
             axes=FALSE, main="Dependogram")
        axis(1,at=1:l,label=test$subsets)
        axis(2)
        points(1:l,rep(1 - test$beta, l),pch=20)
    }
    else
    {
        plot(c(1,1:l),c(0,test$statistics), type="h", ylab="statistic per subset", xlab="subsets",
             axes=FALSE, main="Dependogram")
        axis(1,at=1:l,label=test$subsets)
        axis(2)
        points(1:l,test$critical.values,pch=20)
    }
}

##############################################################################
