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


##' Test of exchangeability for bivariate EV copulas based on the Pickands
##' or CFG estimators -- see SJS paper
##'
##' @title Test of exchangeability based on An
##' @param x the data
##' @param N number of multiplier replications
##' @param estimator "Pickands" or "CFG"
##' @param ties logical indicating whether ties are present
##' @param ties.method passed to pobs
##' @param derivatives based on "An" or "Cn"
##' @param m grid size
##' @return an object of class 'htest'
##' @author Ivan Kojadinovic
exchEVTest <- function(x, N = 1000, estimator = c("CFG", "Pickands"),
                       ties = NA, ties.method = eval(formals(rank)$ties.method), m = 100,
                       derivatives = c("Cn", "An")) {

    ## Checks
    stopifnot(N >= 1L && m >= 5L)
    if(!is.matrix(x)) {
        warning("coercing 'x' to a matrix.")
        stopifnot(is.matrix(x <- as.matrix(x)))
    }
    estimator <- match.arg(estimator)
    ties.method <- match.arg(ties.method)
    derivatives <- match.arg(derivatives)

    ## Make pseudo-observations
    n <- nrow(x)
    u <- pobs(x, ties.method = ties.method)

    ## Make grid
    g <- seq(1/m, 0.5, len = m)

    ## Compute the test statistic
    s <- .C(evsymtest_stat,
            as.double(-log(u[,1])),
            as.double(-log(u[,2])),
            as.integer(n),
            as.double(g),
            as.integer(m),
            as.integer(estimator == "CFG"),
            stat = double(1))$stat

    ## Ties: by default, if at least one column has at least one duplicated entry
    if (is.na(ties <- as.logical(ties))) {
	ties <- any(apply(x, 2, anyDuplicated))
        if (ties)
            warning("argument 'ties' set to TRUE")
    }


    ## Compute approximate realizations under the null
    ## If there are ties, use an adapted bootstrap
    if (ties) {

        ## Get ties structure from initial sample
        ir <- apply(u, 2, function(y) rank(sort(y)))

        ## One replication
        one <- function() {
            ## The bootstrapped sample
            u.b <- u
            for (i in 1:n) {
                s <- sample(1:2)
                u.b[i,1] <- u[i,s[1]]
                u.b[i,2] <- u[i,s[2]]
            }
            ## Apply tie structure
            for (i in 1:2) {
                u.b <- u.b[order(u.b[,i]),]
                u.b[,i] <- u.b[ir[,i], i]
            }
            ## Compute pseudo-observations
            u.b <- pobs(u.b, ties.method = ties.method)

            ## Compute the test statistic
            .C(evsymtest_stat,
               as.double(-log(u.b[,1])),
               as.double(-log(u.b[,2])),
               as.integer(n),
               as.double(g),
               as.integer(m),
               as.integer(estimator == "CFG"),
               stat = double(1))$stat
        }
        ## N replications
        s0 <- replicate(N, one())
    }
    else { # If there are no ties, use the multiplier bootstrap
        if (derivatives == "Cn")
            s0 <- .C(evsymtest,
                     as.double(u[,1]),
                     as.double(u[,2]),
                     as.integer(n),
                     as.double(g),
                     as.integer(m),
                     as.integer(estimator == "CFG"),
                     as.integer(N),
                     s0 = double(N))$s0
        else
            s0 <- .C(evsymtest_derA,
                     as.double(u[,1]),
                     as.double(u[,2]),
                     as.integer(n),
                     as.double(g),
                     as.integer(m),
                     as.integer(estimator == "CFG"),
                     as.integer(N),
                     s0 = double(N))$s0
    }

    structure(class = "htest",
              list(method = paste("Test of exchangeability for bivariate extreme-value copulas with argument 'estimator' set to '", estimator, "' and argument 'm' set to ", m, sep=""),
                   statistic = c(statistic = s),
                   p.value = (sum(s0 >= s) + 0.5) / (N + 1),
                   data.name = deparse(substitute(x))))
}

##' Test of exchangeability for bivariate copulas based on the
##' empirical copula -- see SJS paper
##'
##' @title Test of exchangeability based on Cn
##' @param x the data
##' @param N number of multiplier replications
##' @param m grid size; if 0, use pseudo-observations
##' @param ties logical indicating whether ties are present
##' @param ties.method passed to pobs
##' @return an object of class 'htest'
##' @author Ivan Kojadinovic
exchTest <- function(x, N = 1000, ties = NA,
                     ties.method = eval(formals(rank)$ties.method), m = 0) {

    ## Checks
    stopifnot(N >= 1L && m >= 0L)
    if(!is.matrix(x)) {
        warning("coercing 'x' to a matrix.")
        stopifnot(is.matrix(x <- as.matrix(x)))
    }

    ## Make pseudo-observations
    n <- nrow(x)
    u <- pobs(x, ties.method = ties.method)

    ## Make grid
    if (m > 0L) {
        xis <- yis <- seq(1/m, 1 - 1/m, len = m)
        g <- as.matrix(expand.grid(xis, yis, KEEP.OUT.ATTRS = FALSE))
        ng <- m^2
    } else {
        g <- u
        ng <- n
    }

    ## Compute the test statistic
    s <- .C(exchtestCn_stat,
            as.double(u[,1]),
            as.double(u[,2]),
            as.integer(n),
            as.double(g[,1]),
            as.double(g[,2]),
            as.integer(ng),
            stat = double(1))$stat

    ## Ties: by default, if at least one column has at least one duplicated entry
    if (is.na(ties <- as.logical(ties))) {
	ties <- any(apply(x, 2, anyDuplicated))
        if (ties)
            warning("argument 'ties' set to TRUE")
    }

    ## Compute approximate realizations under the null
    ## If there are ties, use an adapted bootstrap
    if (ties) {

        ## Get ties structure from initial sample
        ir <- apply(u, 2, function(y) rank(sort(y)))

        ## One replication
        one <- function() {
            ## The bootstrapped sample
            u.b <- u
            for (i in 1:n) {
                s <- sample(1:2)
                u.b[i,1] <- u[i,s[1]]
                u.b[i,2] <- u[i,s[2]]
            }
            ## Apply tie structure
            for (i in 1:2) {
                u.b <- u.b[order(u.b[,i]),]
                u.b[,i] <- u.b[ir[,i], i]
            }
            ## Compute pseudo-observations
            u.b <- pobs(u.b, ties.method = ties.method)

            ## Make grid if necessary
            if (m == 0L)
                g <- u.b

            ## Compute the test statistic from the bootstrapped sample
            .C(exchtestCn_stat,
               as.double(u.b[,1]),
               as.double(u.b[,2]),
               as.integer(n),
               as.double(g[,1]),
               as.double(g[,2]),
               as.integer(ng),
               stat = double(1))$stat
        }
        ## N replications
        s0 <- replicate(N, one())
    }
    else # If there are no ties, use the multiplier bootstrap
        s0 <- .C(exchtestCn,
                 as.double(u[,1]),
                 as.double(u[,2]),
                 as.integer(n),
                 as.double(g[,1]),
                 as.double(g[,2]),
                 as.integer(ng),
                 as.integer(N),
                 s0 = double(N))$s0

    structure(class = "htest",
              list(method = paste("Test of exchangeability for bivariate copulas with argument 'm' set to",m),
                   statistic = c(statistic = s),
                   p.value = (sum(s0 >= s) + 0.5) / (N + 1),
                   data.name = deparse(substitute(x))))
}


##' Test of radial symmetry based on the empirical copula
##'
##' @title Test of radial symmetry based on Cn
##' @param x the data
##' @param N number of multiplier replications
##' @param ties logical indicating whether ties are present
##' @return an object of class 'htest'
##' @author Ivan Kojadinovic
radSymTest <- function(x, N = 1000, ties = NA) {

    ## Checks
    stopifnot(N >= 1L)
    if(!is.matrix(x)) {
        warning("coercing 'x' to a matrix.")
        stopifnot(is.matrix(x <- as.matrix(x)))
    }

    n <- nrow(x)
    p <- ncol(x)

    ## Make pseudo-observations
    u <- pobs(x)

    ## Make grid
    #g <- u
    #ng <- n

    ## Compute the test statistic
    s <- .C(radsymtestCn_stat,
            as.double(u),
            as.integer(n),
            as.integer(p),
            as.double(u),
            as.integer(n),
            stat = double(1))$stat

    ## Ties: by default, if at least one column has at least one duplicated entry
    if (is.na(ties <- as.logical(ties))) {
	ties <- any(apply(x, 2, anyDuplicated))
        if (ties)
            warning("argument 'ties' set to TRUE")
    }

    ## If ties, get ties structure from initial sample
    if (ties)
        ir <- apply(u, 2, function(y) rank(sort(y)))


    ## One replication
    one <- function() {

        ## The bootstrapped sample
        u.b <- u
        sel <- (runif(n) < 0.5)
        u.b[sel,] <- 1 - u[sel,]

        if (ties)
            ## Apply tie structure
            for (i in 1:p) {
                u.b <- u.b[order(u.b[,i]),]
                u.b[,i] <- u.b[ir[,i], i]
            }

        ## Compute pseudo-observations
        u.b <- pobs(u.b)

        ## Make grid
        #g <- u.b

        ## Compute the test statistic from the bootstrapped sample
        .C(radsymtestCn_stat,
           as.double(u.b),
           as.integer(n),
           as.integer(p),
           as.double(u.b),
           as.integer(n),
           stat = double(1))$stat
    }

    ## N replications
    s0 <- replicate(N, one())

    structure(class = "htest",
              list(method = "Test of radial symmetry based on the empirical copula",
                   statistic = c(statistic = s),
                   p.value = (sum(s0 >= s) + 0.5) / (N + 1),
                   data.name = deparse(substitute(x))))
}
