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


##' @title Computing Conditional Copulas C_{j|1,..,j-1}(u_j | u_1,..,u_{j-1})
##' @param u A data matrix in [0,1]^(n, d) of U[0,1]^d samples
##' @param copula An object of class Copula
##' @param indices A vector of indices j in {1,..,d} for which
##'        C_{j|1,..,j-1}(u_j | u_1,..,u_{j-1}) is computed
##' @param log A logical indicating whether the log-transform is computed
##' @param drop logical indicating whether one-column matrices are returned
##'        as vectors
##' @param ... Additional arguments
##' @return An (n, |indices|) matrix U of supposedly multivariate uniformly
##'         distributed realizations (or the log of
##'         the result if log = TRUE)
##' @author Marius Hofert and Martin Maechler
##' @note Hidden; call via public cCopula() (for having arguments tested)
rosenblatt <- function(u, copula, indices = 1:dim(copula), log = FALSE, drop = TRUE, ...)
{
    if (!is.matrix(u)) # minimal checking here!
        u <- rbind(u, deparse.level = 0L)
    n <- nrow(u) # sample size

    if(is(copula, "indepCopula")) {

        U <- if(log) log(u) else u
        U[, indices, drop = drop]

    } else if(is(copula, "moCopula")) {

        if(dim(copula) != 2)
            stop("cCopula() is currently only available for bivariate moCopula objects")
        u1 <- u[,1, drop = FALSE]
        u2 <- u[,2, drop = FALSE]
        a <- attr(copula, "parameters") # alpha_1, alpha_2
        ## Fix u1 and consider C_2(u_2|u_1), which is given by (1-a[1])*u1^(-a[1])*u2
        ## if u2 < u1^(a[1]/a[2]) and by u2^(1-a[2]) if u2 >= u1^(a[1]/a[2]). Note that
        ## C_2(u_2|u_1) jumps in u1^(a[1]/a[2]) from (1-a[1])*u1. to u1. (see below)
        ## which is why we implement the generalized distributional transform to
        ## get uniformity in the Rosenblatt transformed sample. Also note that
        ## if u2 == u1^(a[1]/a[2]), we need to work with all.equal() due to exact
        ## equality not picking up all samples on the singular component (see, e.g.,
        ## moCopula(c(0.9, 0.75))).
        ii1 <- 0 <= u2 & u2 < u1^(a[1]/a[2])
        ii2 <- u1^(a[1]/a[2]) < u2 & u2 <= 1
        ## Determine singular component
        all.all.equal <- function(x)
            apply(x, 1, function(x.) isTRUE(all.equal(x.[1],x.[2], check.attributes = FALSE)))
        ii3 <- all.all.equal(cbind(u1^a[1], u2^a[2]))
        U <- matrix(, nrow = length(u1), ncol = 2)
        U[,1] <- u1
        U[ii1, 2] <- (1-a[1]) * u1[ii1]^(-a[1]) * u2[ii1]
        U[ii2, 2] <- u2[ii2]^(1-a[2])
        u1. <- u1[ii3]^(a[1]*(1/a[2] - 1))
        U[ii3, 2] <- runif(sum(ii3), min = (1-a[1])*u1., max = u1.) # singular component (do generalized distributional transform idea)
        if(log) U <- log(U)
        U[, indices, drop = drop]

    } else if(is(copula, "normalCopula")) { # Gauss copula (see, e.g., Cambou, Hofert, Lemieux)

        max.ind <- tail(indices, n = 1) # maximal index
        P <- getSigma(copula) # (d, d)-matrix
        x <- qnorm(u[, 1:max.ind, drop = FALSE]) # compute all 'x' we need
        C.j <- function(j) # C_{j|1,..,j-1}(u_j | u_1,..,u_{j-1})
        {
            if(j == 1) {
                if(log) log(u[,1]) else u[,1]
            } else {
                ij1 <- seq_len(j - 1L)
                P. <- P[j, ij1, drop = FALSE] %*% solve(P[ij1, ij1, drop = FALSE]) # (1, j-1) %*% (j-1, j-1) = (1, j-1)
                mu.cond <- as.numeric(P. %*% t(x[, ij1, drop = FALSE])) # (1, j-1) %*% (j-1, n) = (1, n) = n
                P.cond <- P[j,j] - P. %*% P[ij1, j, drop = FALSE] # (1, 1) - (1, j-1) %*% (j-1, 1) = (1, 1)
                pnorm(x[,j], mean = mu.cond, sd = sqrt(P.cond), log.p = log)
            }
        }

        ## Return
        if(drop || n > 1)
            vapply(indices, C.j, numeric(n))
        else # n == 1, !drop
            rbind(vapply(indices, C.j, 1.))#, deparse.level = 0L)

    } else if(is(copula, "tCopula")) { # t copula (see, e.g., Cambou, Hofert, Lemieux)

        max.ind <- tail(indices, n = 1) # maximal index
        P <- getSigma(copula) # (d, d)-matrix
        nu <- getdf(copula) # degrees of freedom
        x <- qt(u[, 1:max.ind, drop = FALSE], df = nu) # compute all 'x' we need
        C.j <- function(j) # C_{j|1,..,j-1}(u_j | u_1,..,u_{j-1})
        {
            if(j == 1) {
                if(log) log(u[,1]) else u[,1]
            } else {
                P1.inv <- solve(P[1:(j-1), 1:(j-1), drop = FALSE])
                x1 <- x[, 1:(j-1), drop=FALSE]
                g  <- vapply(1:n, function(i) x1[i, ,drop = FALSE] %*% P1.inv %*%
                                              t(x1[i, ,drop = FALSE]), numeric(1))
                P.inv <- solve(P[1:j, 1:j, drop = FALSE])
                s1 <- sqrt((nu + j - 1) / (nu + g))
                s2 <- (x1 %*% P.inv[1:(j-1), j, drop = FALSE]) / sqrt(P.inv[j,j])
                lres <- pt(s1 * (sqrt(P.inv[j, j]) * x[,j, drop = FALSE] + s2),
                           df = nu+j-1, log.p = TRUE)
                if(log) lres else exp(lres)
            }
        }

        ## Return
        if(drop || n > 1)
            vapply(indices, C.j, numeric(n))
        else # n == 1, !drop
            rbind(vapply(indices, C.j, 1.))#, deparse.level = 0L)

    } else if((NAC <- is(copula, "outer_nacopula")) ||
              is(copula, "archmCopula")) { # (nested) Archimedean copulas

        ## Dealing with NACs and the two classes of Archimedean copulas
	if(NAC) {
	    if(length(copula@childCops))
		stop("Currently, only Archimedean copulas are supported")
            ## outer_nacopula but with no children => an AC => continue
	    cop <- copula@copula # class(cop) = "acopula"
	    th <- cop@theta
	} else { # class(cop) = "archmCopula"
	    th <- copula@parameters
	    cop <- getAcop(copula) # => class(cop) = "acopula" but without parameter or dim
	}
	stopifnot(cop@paraConstr(th, dim(copula)))

        ## Compute conditional probabilities C_{j|1,..,j-1}(u_j | u_1,..,u_{j-1})
	psiI  <- cop@iPsi(u, theta = th) # (n, d) matrix of psi^{-1}(u)
	psiI. <- t(apply(psiI, 1, cumsum)) # corresponding (n, d) matrix of row sums
        ## Note: C_{j|1,..,j-1}(u_j | u_1,...,u_{j-1})
        ##       = \psi^{(j-1)}(\sum_{k=1}^j \psi^{-1}(u_k)) /
        ##         \psi^{(j-1)}(\sum_{k=1}^{j-1} \psi^{-1}(u_k))
        C.j <- function(j) # C_{j|1,..,j-1}(u_j | u_1,..,u_{j-1})
        {
            if(j == 1) {
                if(log) log(u[,1]) else u[,1]
            } else {
		tt <- as.vector(psiI.[, c(j, j-1)])
		if(th < 0) { ## allowed for dim==2 and  {AMH, Clayton, Frank}:
		    D <- cop@absdPsi(tt, theta = th, degree = j-1, log = FALSE)
		    rat <- D[1:n] / D[(n+1):(2*n)]
		    if(log) log(rat) else rat
		} else {
		    logD <- cop@absdPsi(tt, theta = th, degree = j-1, log = TRUE)
		    res <- logD[1:n] - logD[(n+1):(2*n)]
		    if(log) res else exp(res)
		}
            }
        }

        ## Return
        if(drop || n > 1)
            vapply(indices, C.j, numeric(n))
        else # n == 1, !drop
            rbind(vapply(indices, C.j, 1.))#, deparse.level = 0L)

    } else if(is(copula, "rotCopula")) {

        ## Compute the Rosenblatt transform of a rotCopula object
        ## Note: Similar to the pCopula method for rotCopula objects
        u.flip <- apply.flip(u, copula@flip) # u with columns suitably flipped so that cCopula(u.flip) = rotated conditional copula
        res <- cCopula(u.flip, copula = copula@copula, indices = indices, log = log, drop = drop)
        ## Note: The first component of 'u' should be the original one (not changed)
        ##       as the 'change' (going to u.flip) was only done to compute the
        ##       conditional copula functions of a rotated copula.
        if(1 %in% indices) res[,1] <- u[,1]
        res

    } else if(is(copula, "mixCopula")) {

        ## Check
        ## Note: For d > 2, note that Schmitz' formula for conditional copulas is a
        ##       fraction of weighted sums of copulas and thus not equal to a weighted
        ##       sum of fractions (unless d = 2 in which case the denominators are all 1).
        if(dim(copula) != 2)
            stop("cCopula() is currently only available for bivariate mixCopula objects")

        ## Compute the Rosenblatt transform of a bivariate mixture copula
        ## Note: Similar to the pCopula method for mixCopula objects.
        w <- copula@w # vector of mixture weights
        m <- length(w) # number of copulas in the mixture
        len.ind <- length(indices)
        w.cond.cop <- vapply(1:m, FUN = function(k) # (n, len.ind, m)-array
            w[k] * cCopula(u, copula = copula@cops[[k]], indices = indices,
                           inverse = FALSE, log = FALSE, drop = FALSE),
                           FUN.VALUE = matrix(NA_real_, nrow = n,
                                              ncol = len.ind), ...)
        wccop <- array(w.cond.cop, dim = c(n, len.ind, m)) # make sure to indeed have an array (e.g., if n = 1)
        res <- apply(wccop, 1:2, sum) # aggregate over all components of the mixCopula object
        if(log) res <- log(res)

        ## Return (note: correctly returns a vector if 'drop = TRUE' and n = 1)
        if(drop && n == 1) drop(res) else res

    } else {
	stop("Not yet implemented for copula class ", class(copula))
    }
}


##' @title Computing Conditional Copula Quantile Functions C^-_{j|1,..,j-1}(u_j | u_1,..,u_{j-1})
##' @param u A data matrix in [0,1]^(n, d) of (pseudo-/copula-)observations
##' @param copula An object of class Copula
##' @param indices A vector of indices j in {1,..,d} for which
##'        C^-_{j|1,..,j-1}(u_j | u_1,..,u_{j-1}) is computed
##' @param log A logical indicating whether the log-transform is computed
##' @param drop logical indicating whether one-column matrices are returned
##'        as vectors
##' @param ... Additional arguments
##' @return An (n, |indices|) matrix U of copula distributed samples
##'         (or the log of the result if log = TRUE)
##' @author Marius Hofert and Martin Maechler
##' @note Call this via cCopula() (for having arguments tested)
iRosenblatt <- function(u, copula, indices = 1:dim(copula), log = FALSE, drop = TRUE, ...)
{
    if (!is.matrix(u)) # minimal checking here!
        u <- rbind(u, deparse.level = 0L)
    n <- nrow(u) # sample size
    if(is(copula, "indepCopula")) {

        U <- if(log) log(u) else u

    } else if(is(copula, "moCopula")) {

        if(dim(copula) != 2)
            stop("cCopula(, inverse = TRUE) is currently only available for bivariate moCopula objects")
        u1 <- u[,1, drop = FALSE]
        u2 <- u[,2, drop = FALSE]
        a <- attr(copula, "parameters") # alpha_1, alpha_2
        ## See Cambou et al. (2017)
        u1. <- u[1]^(a[1] * (1/a[2]-1))
        ii1 <- 0 <= u2 & u2 <= (1-a[1]) * u1.
        ii2 <- (1-a[1]) * u1. < u2 & u2 < u1.
        ii3 <- u1. <= u2 & u2 <= 1
        U <- matrix(, nrow = length(u1), ncol = 2)
        U[,1] <- u1
        U[ii1,2] <- u1[ii1]^a[1] * u2[ii1] / (1-a[1])
        U[ii2,2] <- u1[ii2]^(a[1]/a[2])
        U[ii3,2] <- u2[ii3]^(1/(1-a[2]))
        if(log) U <- log(U)

    } else if(is(copula, "normalCopula")) { # Gauss copula (see, e.g., Cambou, Hofert, Lemieux)

        P <- getSigma(copula) # (d, d)-matrix, always symmetric
        U <- u # consider u as U[0,1]^d
        max.ind <- tail(indices, n = 1) # maximal index
        x <- qnorm(u[, 1:max.ind, drop = FALSE]) # will be updated
        for(j in seq_len(max.ind)) {
            if(j == 1) {
                if(log) U[,1] <- log(U[,1]) # adjust the first column for 'log'
            } else { # j >= 2
                ji <- seq_len(j-1L)
                P. <- P[j, ji, drop = FALSE] %*% solve(P[ji, ji, drop = FALSE]) # (1, j-1) %*% (j-1, j-1) = (1, j-1)
                mu.cond <- as.numeric(tcrossprod(P., x[, ji, drop = FALSE])) # (1, j-1) %*% (j-1, n) = (1, n) = n
                P.cond <- P[j, j] - P. %*% P[ji, j, drop = FALSE] # (1, 1) - (1, j-1) %*% (j-1, 1) = (1, 1)
                U[,j] <- pnorm(qnorm(u[, j], mean = mu.cond, sd = sqrt(P.cond)), log.p = log)
                x[,j] <- qnorm(if(log) exp(U[,j]) else U[,j]) # update x[,j]
            }
        }

    } else if(is(copula, "tCopula")) { # t copula (see, e.g., Cambou, Hofert, Lemieux)

        P <- getSigma(copula) # (d, d)-matrix, always symmetric
        nu <- getdf(copula) # degrees of freedom
        U <- u # consider u as U[0,1]^d
        max.ind <- tail(indices, n = 1) # maximal index
        x <- qt(u[, 1:max.ind, drop = FALSE], df = nu) # will be updated
        for(j in seq_len(max.ind)) {
            if(j == 1) {
                if(log) U[,1] <- log(U[,1]) # adjust the first column for 'log'
            } else { # j >= 2
                ji <- seq_len(j-1L)
                ## compute solve(), as we need it n times:
                P1.inv <- solve(P[ji, ji, drop = FALSE])
                x1 <- x[, ji, drop = FALSE]
                g  <- vapply(1:n, function(i) { x1i <- x1[i, , drop = FALSE] #  X1 P1i X1' :
                                                x1i %*% tcrossprod(P1.inv, x1i) },
                             numeric(1))
                P.inv <- solve(P[1:j, 1:j, drop = FALSE])
                sPj <- sqrt(P.inv[j, j])
                s1 <- sqrt((nu + j - 1) / (nu + g))
                s2 <- x1 %*% P.inv[ji, j, drop = FALSE]
                U[,j] <- pt((qt(u[, j], df = nu+j-1)/s1 - s2/sPj) / sPj,
                            df = nu, log.p = log)
                x[,j] <- qt(if(log) exp(U[,j]) else U[,j], df = nu) # update x[,j]
            }
        }

    } else if((NAC <- is(copula, "outer_nacopula")) ||
              is(copula, "archmCopula")) { # (nested) Archimedean copulas

        ## Dealing with NACs and the two classes of Archimedean copulas
        if(NAC) {
	    if(length(copula@childCops))
		stop("Currently, only Archimedean copulas are supported")
            ## outer_nacopula but with no children => an AC => continue
	    cop <- copula@copula # class(cop) = "acopula"
	    th <- cop@theta
	} else { # class(cop) = "archmCopula"
	    th <- copula@parameters
	    cop <- getAcop(copula) # => class(cop) = "acopula" but without parameter or dim
	}
	stopifnot(cop@paraConstr(th, dim(copula)))

        ## Compute conditional quantiles C^-_{j|1,..,j-1}(u_j | u_1,..,u_{j-1})
        U <- u # u's are supposedly U[0,1]^d
        max.ind <- tail(indices, n = 1) # maximal index
        if(cop@name == "Clayton") { # Clayton case (explicit)
            sum. <- U[,1]^(-th)
            if(max.ind >= 2) {
                for(j in 2:max.ind) {
                    U[,j] <- log1p((1-j+1+sum.)*(u[,j]^(-1/(j-1+1/th)) - 1))/(-th)
                    eUj <- exp(U[,j])
                    sum. <- sum. + eUj^(-th)
                    if(!log) U[,j] <- eUj
                }
            }
        } else { # general case (non-Clayton)
            ## TODO: After the acopula and archmCopula classes are better merged, the
            ##       tedious conversion to acopula (see above) and then again to
            ##       archmCopula (below) is not required anymore.
            arCop <- archmCopula(cop@name, param = th)
            ## f(x) := C_{j|1,..,j-1}(x | u_{i1},..,u_{i j-1}) - u_{ij}
            if(max.ind >= 2) {
                f <- function(x, U.i1_to_jm1, u.ij)
                    rosenblatt(c(U.i1_to_jm1, x), copula = arCop, indices = j, drop=TRUE) - u.ij
                for(j in 2:max.ind) {
                    ## Precompute quantities from the jth column so that uniroot() is faster
                    U..1_to_jm1 <- U[, seq_len(j-1L), drop = FALSE] # matrix U_{., 1:(j-1)}
                    u..j <- u[,j] # vector u_{., j}
                    ## Iterate over samples and find the root
                    for(i in 1:n)
                        U[i,j] <- uniroot(f, interval = 0:1, U.i1_to_jm1 = U..1_to_jm1[i,],
                                          u.ij = u..j[i], ...)$root
                }
            }
            if(log) U <- log(U)
        }

    } else {
	stop("Not yet implemented for copula class ", class(copula))
    }

    ## Return
    U[, indices, drop = drop]
}


##' @title Conditional Copulas and Their Inverses
##' @param u A data matrix in [0,1]^(n, d) of U[0,1]^d samples if inverse = FALSE
##'        and ((pseudo-/copula-)observations if inverse = TRUE
##' @param copula object of class Copula
##' @param indices vector of indices j in {1,..,dim(copula)} for which
##'        C_{j|1,..,j-1}(u_j | u_1,..,u_{j-1}) (or its inverse
##'        C^-_{j|1,..,j-1}(u_j | u_1,..,u_{j-1}) if inverse = TRUE) is computed.
##' @param inverse logical indicating whether the inverse
##'        C^-_{j|1,..,j-1}(u_j | u_1,..,u_{j-1}) is computed (known as
##'        'conditional distribution method' for sampling)
##' @param log logical indicating whether the log-transform is computed
##' @param drop logical indicating whether one-column matrices are returned
##'        as vectors
##' @param ... Additional arguments passed to the underlying
##'         rosenblatt() and iRosenblatt()
##' @return An (n, d) matrix U of supposedly U[0,1]^d realizations
##'         or copula samples [or the log of the result if log = TRUE]
##' @author Marius Hofert and Martin Maechler ('drop')
##' @note This is just a wrapper (including argument testing) of rosenblatt()
##'       and iRosenblatt();   ___HELP___ in ../man/cCopula.Rd
cCopula <-  function(u, copula, indices = 1:dim(copula), inverse = FALSE,
                     log = FALSE, drop = FALSE, ...)
{
    ## Argument checks
    if(!is.matrix(u)) u <- rbind(u, deparse.level = 0L)
    d <- ncol(u)
    stopifnot(0 <= u, u <= 1, d >= 2, is(copula, "Copula"),
              is.logical(inverse), is.logical(log))
    if(!missing(indices)) {
        if(!all(1 <= indices & indices <= dim(copula)))
            stop("'indices' have to be between 1 and the copula dimension.")
        if(is.unsorted(indices))
            stop("'indices' have to be unique and given in increasing order.")
        if(indices[length(indices)] > d)
            stop("The maximal index must be less than or equal to the number of columns of 'u'")
    }
    ## Call work horses
    if(inverse)
        iRosenblatt(u, copula=copula, indices=indices, log=log, drop=drop, ...)
    else rosenblatt(u, copula=copula, indices=indices, log=log, drop=drop, ...)
}
