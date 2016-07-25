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

##################################################################################
### Partial derivatives of the CDF wrt arguments
##################################################################################

setGeneric("dCdu", function(copula, u, ...) standardGeneric("dCdu"))

##' @title List of control arguments for grad() / jacobian() in package numDeriv
##' @param eps see ?grad
##' @param d see ?grad
##' @param zero.tol see ?grad
##' @param r see ?grad
##' @param v see ?grad
##' @param show.details see ?grad
##' @return a list of the above
##' @author Ivan Kojadinovic
gradControl <- function(eps = 1e-4, d = 1e-4,
                        zero.tol = sqrt(.Machine$double.eps / 7e-7),
                        r = 6, v = 2, show.details = FALSE) {
    list(eps = eps, d = d, zero.tol = zero.tol, r = r, v = v,
         show.details = show.details)
}


##' @title Returns 'side' vector for grad() / jacobian() when computing
##' partial derivatives at x in (lb, ub) or [lb, ub] or ...
##' @param x point where partial derivatives will be computed
##' @param dimx length of x
##' @param d parameter from grad(); see ?grad
##' @param lb lower bound for x values
##' @param ub upper bound for x values
##' @return the 'side' vector
##' @author Ivan Kojadinovic
sides <- function(x, dimx, d, lb, ub) {
    res <- rep(NA, dimx) # two-sided derivatives by default
    res[x - lb < d] <- 1 # right derivative
    res[ub - x < d] <- -1 # left derivative
    res
}

##' @title Basic implementation of dCdu based on numerical differentiation
##' @param copula an object of class "copula"
##' @param u evaluation point in [0,1]^d or matrix of evaluation points
##' @param method.args to be passed to grad()
##' @param ...
##' @return partial derivatives of the copula wrt arguments at u
##' @author Ivan Kojadinovic
dCduCopulaNum <- function(copula, u, method.args = gradControl(d = 1e-1), ...) {
    warning("Function dCdu() not implemented for copulas of class '", class(copula),
            "'; numerical differentiation used")
    logC <- function(x) log(pCopula(x, copula))
    dim <- copula@dimension
    res <- t(apply(u, 1, function(u.)
        grad(logC, u., side = sides(u., dim, method.args$d, 0, 1),
             method.args = method.args))) * pCopula(u, copula)
    res[res < 0] <- 0 # because 0 < dCdu < 1 uniformly
    res[res > 1] <- 1
    res
}

## For archmCopula objects
dCduArchmCopula <- function(copula, u, ...) {
    ## TODO: For ACs, the following is better than 'fixed formulas' => use it
    ## require(copula)
    ## n <- 10
    ## d <- 4
    ## family <- "Gumbel"
    ## th <- 2
    ## u <- matrix(runif(n*d), ncol=d)
    ## cop <- onacopulaL(family, nacList=list(th, seq_len(d)))
    ## iPsi.u <- cop@copula@iPsi(u, theta=th)
    iPsi.u <- iPsi(copula, u)
    d <- copula@dimension
    j <- ceiling(d/2)
    ## NEED function absdPsi to make below work for any Archimedean copula in the same way that iPsi exists
    ## dCdu <- sapply(seq_len(d), function(j) exp(cop@copula@absdPsi(rowSums(iPsi.u), theta=th, degree=1, log=TRUE) - cop@copula@absdPsi(iPsi.u[,j], theta=th, degree=1, log=TRUE)))
}

## For copulas with explicit cdf
## Warning: This function assumes symmetry in u
dCduExplicitCopula <- function(copula, u, ...) {
    dim <- copula@dimension
    mat <- matrix(NA_real_, nrow(u), dim)
    algNm <- paste(class(copula)[1], "cdfDerWrtArg.algr", sep=".")
    ## FIXME: fails if class(copula) *extends* (but is not identical to) say "normalCopula"
    if(exists(algNm)) {
	alpha <- copula@parameters # typically used in 'eval(der.cdf.u, *)' below
        der.cdf.u <- get(algNm)[dim]
        unames0 <- paste0("u",1:dim)
        for (j in 1:dim) {
            unames <- unames0; unames[1] <- unames0[j]; unames[j] <- unames0[1]
            colnames(u) <- unames
            mat[,j] <- eval(der.cdf.u, data.frame(u))
        }
    } else warning("Function dCdu() not implemented for copulas of class '",
                   class(copula), "'")
    mat
}


## This function is used for khoudrajiCopula objects
dCduIndepCopula <- function(copula, u, ...) {
  dim <- copula@dimension
  mat <- matrix(0, nrow(u), dim)
  for (j in 1:dim) {
    mat[,j] <- apply(u[, -j, drop=FALSE], 1, prod)
  }
  mat
}

setMethod("dCdu", signature("copula"), dCduCopulaNum)
setMethod("dCdu", signature("joeCopula"), dCduCopulaNum)
setMethod("dCdu", signature("archmCopula"), dCduExplicitCopula)
setMethod("dCdu", signature("plackettCopula"), dCduExplicitCopula)
setMethod("dCdu", signature("evCopula"), dCduExplicitCopula)
setMethod("dCdu", signature("gumbelCopula"), dCduExplicitCopula)
setMethod("dCdu", signature("indepCopula"), dCduIndepCopula)

## For ellipCopula objects
dCduEllipCopula <- function(copula, u, ...) {
    dim <- copula@dimension
    sigma <- getSigma(copula)

    ## quantile transformation
    v <- switch(class(copula),
		"normalCopula" = qnorm(u),
		"tCopula" = {
                    df <- getdf(copula) # (needed further)
                    qt(u, df=df)
                },
                stop("not implemented for class ", class(copula)))
    n <- nrow(u)
    mat <- matrix(0,n,dim)

    for (j in 1:dim) {
	s <- sigma[-j,-j] - sigma[-j,j,drop=FALSE] %*% sigma[j,-j,drop=FALSE]

	switch(class(copula),
	       "normalCopula" =
		   {
		       if (dim == 2) {
			   rho <- copula@parameters[1] # IK: single parameter only
			   mat[,j] <- pnorm(v[,-j], rho * v[,j], sqrt(1 - rho^2))
		       }
		       else
			   for (i in 1:n)
			       mat[i,j] <- pmvnorm(lower = rep(-Inf, dim - 1), upper = v[i,-j],
						   mean = v[i,j] * sigma[-j,j],
						   sigma = drop(s))
		   },
	       "tCopula" =
		   {
		       if (dim == 2) {
			   rho <- copula@parameters[1] # IK: single parameter only
			   mat[,j] <-  pt(sqrt((df+1)/(df+v[,j]^2)) / sqrt(1 - rho^2)
					  * (v[,-j] - rho * v[,j]), df=df+1)
		       }
		       else {
			   if(df != as.integer(df))
			       stop("'df' is not integer; therefore, dCdu() cannot be computed yet")
			   for (i in 1:n)
			       mat[i,j] <- pmvt(lower = rep(-Inf, dim - 1),
						upper = drop(sqrt((df+1)/(df+v[i,j]^2)) *
							     (v[i,-j] - v[i,j] * sigma[-j,j])),
						sigma = s, df = df + 1)
		       }

		   })
    }
    mat
}

setMethod("dCdu", signature("ellipCopula"), dCduEllipCopula)

##################################################################################
### Plackett formula for elliptical copulas
##################################################################################

setGeneric("plackettFormulaDim2", function(copula, x) standardGeneric("plackettFormulaDim2"))

##' @title Derivative of the bivariate standard normal cdf wrt correlation
##' @param copula a bivariate normalCopula object
##' @param x evaluation point
##' @return the derivative at x
##' @author Ivan Kojadinovic
plackettFormulaDim2NormalCopula <- function(copula, x) {
    rho <- copula@parameters[1] # JY: single parameter only
    ir2 <- 1 - rho^2
    as.matrix(exp(-(x[,1]^2 + x[,2]^2 - 2 * rho * x[,1] * x[,2]) /
                   (2 * ir2)) /
              (2 * pi * sqrt(ir2)))
}

setMethod("plackettFormulaDim2", signature("normalCopula"), plackettFormulaDim2NormalCopula)

##' @title Derivative of the bivariate standard t cdf  wrt correlation
##' @param copula a bivariate tCopula object
##' @param x evaluation point
##' @return the derivative at x
##' @author Ivan Kojadinovic
plackettFormulaDim2TCopula <- function(copula, x) {
    rho <- copula@parameters[1] #JY: single parameter only
    ir2 <- 1 - rho^2
    df <- getdf(copula)
    as.matrix((1 + (x[,1]^2 + x[,2]^2 - 2 * rho * x[,1] * x[,2]) /
               (df * ir2))^(-df / 2) / (2 * pi * sqrt(ir2)))
}

setMethod("plackettFormulaDim2", signature("tCopula"), plackettFormulaDim2TCopula)


## JY: is rho a single parameter?
setGeneric("plackettFormula",  function(copula, dim, rho, s, m, x, i, j) standardGeneric("plackettFormula"))

## For the multivariate standard normal cdf
plackettFormulaNormalCopula <- function(copula, dim, rho, s, m, x, i, j) {
    exp(-(x[i]^2 + x[j]^2 - 2 * rho * x[i] * x[j]) /
         (2 * (1 - rho^2))) / (2 * pi * sqrt(1 - rho^2)) *
        (if (dim == 3) pnorm(drop((x[-c(i,j)] - m %*% x[c(i,j)])/sqrt(s)))
         else pmvnorm(lower = rep(-Inf, dim - 2),
                      upper = drop(x[-c(i,j)] - m %*% x[c(i,j)]),
                      sigma = s))
}

setMethod("plackettFormula", signature("normalCopula"), plackettFormulaNormalCopula)

## For the multivariate standard t cdf
plackettFormulaTCopula <- function(copula, dim, rho, s, m, x, i, j) {
    stopifnot(dim >= 3)
    df <- getdf(copula)
    if(df != as.integer(df) && dim > 3)
	stop("'df' is not integer; therefore, plackettFormula() cannot be computed yet")
    term <- 1 + (x[i]^2 + x[j]^2 - 2 * rho * x[i] * x[j]) / (df * (1 - rho^2))
    ## return:
    term^(-df / 2) / (2 * pi * sqrt(1 - rho^2)) *
             if (dim == 3) pt(drop((x[-c(i,j)] - m %*% x[c(i,j)]) / sqrt(term * s)), df= df)
             else pmvt(df = df, lower = rep(-Inf, dim - 2),
                       upper = drop((x[-c(i,j)] - m %*% x[c(i,j)]) / sqrt(term)),
                       sigma = s)
}

setMethod("plackettFormula", signature("tCopula"), plackettFormulaTCopula)

##################################################################################
### Partial derivatives of the CDF wrt parameters
##################################################################################

setGeneric("dCdtheta", function(copula, u, ...) standardGeneric("dCdtheta"))

##' @title Basic implementation of dCdtheta based on numerical differentiation
##' @param copula an object of class "copula"
##' @param u evaluation point in [0,1]^d or matrix of evaluation points
##' @param method.args to be passed to grad()
##' @param ...
##' @return n by p, partial derivatives of the copula wrt parameters at u
##' @author Ivan Kojadinovic
dCdthetaCopulaNum <- function(copula, u, method.args = gradControl(d = 1e-1), ...) {
    warning("Function dCdtheta() not implemented for copulas of class '",
            class(copula), "'; numerical differentiation used")
    logC <- function(theta) {
        freeParam(copula) <- theta
        log(pCopula(u, copula))
    }
    theta <- getParam(copula)
    p <- length(theta)
    lb <- attr(theta, "param.lowbnd")
    ub <- attr(theta, "param.upbnd" )
    jacobian(logC, theta, side = sides(theta, p, method.args$d, lb, ub),
             method.args = method.args) * pCopula(u, copula)
}

## For copulas with explicit cdf
dCdthetaExplicitCopula <- function(copula, u, ...) {
    dim <- copula@dimension
    algNm <- paste(class(copula)[1], "cdfDerWrtPar.algr", sep=".")
    alpha <- getParam(copula) # typically used in 'eval(*)' below
    ## JY: alpha is a scalar parameter
    if(exists(algNm)) {
        der.cdf.alpha <- get(algNm)[dim]
        colnames(u) <- paste0("u", 1:dim)
        as.matrix(eval(der.cdf.alpha, data.frame(u)))
    } else {
        warning("Function dCdtheta() not implemented for copulas of class '",
                class(copula), "'")
        matrix(NA_real_, nrow(u), length(alpha)) # copula@parameters))
    }
}

## For evCopula objects
## JY: For single parameter only
dCdthetaEvCopula <- function(copula, u, ...) {
  loguv <- log(u[,1], u[,2])
  w <- log(u[,2]) / loguv
  ## return  der.cdf.alpha
  ## JY: single parameter only
  as.matrix(pCopula(u, copula) * loguv * dAdtheta(copula, w))
}

setMethod("dCdtheta", signature("copula"), dCdthetaCopulaNum)
setMethod("dCdtheta", signature("joeCopula"), dCdthetaCopulaNum)
setMethod("dCdtheta", signature("tevCopula"), dCdthetaCopulaNum)
setMethod("dCdtheta", signature("archmCopula"), dCdthetaExplicitCopula)
setMethod("dCdtheta", signature("plackettCopula"), dCdthetaExplicitCopula)
setMethod("dCdtheta", signature("evCopula"), dCdthetaExplicitCopula)
setMethod("dCdtheta", signature("gumbelCopula"), dCdthetaExplicitCopula)

## For ellipCopula objects
dCdthetaEllipCopula <- function(copula, u, ...) {
    dim <- copula@dimension

    ## quantile transformation
    v <- switch(class(copula),
		"normalCopula" = qnorm(u),
		"tCopula" = qt(u, df = getdf(copula)),
		stop("not implemented for class ", class(copula)))

    val <-
    if (dim == 2)
        plackettFormulaDim2(copula, v)
    else {
        n <- nrow(u)
        sigma <- getSigma(copula)

        if (copula@dispstr %in% c("ex","ar1")) { ## exchangeable or ar1
            rho <- copula@parameters[1] # IK: single parameter only
            r <- matrix(c(1,-rho,-rho,1),2,2)/(1 - rho^2)
            m <- sigma[-c(1,2),c(1,2)] %*% r
            s <- sigma[-c(1,2),-c(1,2)] -
                sigma[-c(1,2),c(1,2)] %*% r %*% sigma[c(1,2),-c(1,2)]

            mat <- matrix(0,n,1)

            if (copula@dispstr == "ex") ## exchangeable
              for (k in 1:n)
                for (j in 1:(dim-1))
                  for (i in (j+1):dim)
                      mat[k,1] <- mat[k,1] +
                          plackettFormula(copula, dim, rho, s, m, v[k,], i, j)
            else ## ar1
              for (k in 1:n)
                for (j in 1:(dim-1))
                  for (i in (j+1):dim)
                    mat[k,1] <- mat[k,1] + (i - j) * rho ^ (i - j - 1) *
				  plackettFormula(copula, dim, rho, s, m, v[k,], i, j)

            mat
        }
	else { # unstructured or toeplitz or ...
            mat <- matrix(0,n,dim*(dim-1)/2)
            l <- 1
            for (j in 1:(dim-1))
                for (i in (j+1):dim)
                {
                    rho <- sigma[i,j]
                    r <- matrix(c(1,-rho,-rho,1),2,2)/(1 - rho^2)
                    m <- sigma[-c(i,j),c(i,j)] %*% r
                    s <- sigma[-c(i,j),-c(i,j)] - m %*% sigma[c(i,j),-c(i,j)]

                    for (k in 1:n)
                        mat[k,l] <- plackettFormula(copula, dim, rho, s, m, v[k,], i, j)
                    l <- l + 1
                }
            if (copula@dispstr == "un") ## unstructured
                mat
            else if (copula@dispstr == "toep") {
                coef <- matrix(0, dim*(dim-1)/2, dim-1)
                for (k in 1:(dim-1)) {
                    m <- row(sigma) == col(sigma) + k
                    coef[,k] <- P2p(m)
                }
                mat %*% coef
            }
            else stop("Not implemented yet for the dispersion structure ",
                      copula@dispstr)
        }
    }
    free <- isFree(copula@parameters)
    if (.hasSlot(copula, "df.fixed")) free <- free[-length(free)]
    val[, free, drop = FALSE]
}

setMethod("dCdtheta", signature("ellipCopula"), dCdthetaEllipCopula)

##################################################################################
### Partial derivatives of the log PDF wrt arguments
##################################################################################

setGeneric("dlogcdu", function(copula, u, ...) standardGeneric("dlogcdu"))

##' @title Basic implementation of dlogcdu based on numerical differentiation
##' @param copula an object of class "copula"
##' @param u evaluation point in [0,1]^d or matrix of evaluation points
##' @param method.args to be passed to grad()
##' @param ...
##' @return partial derivatives of the log copula density wrt arguments at u
##' @author Ivan Kojadinovic
dlogcduCopulaNum <- function(copula, u, method.args = gradControl(d = 1e-5), ...) {
    warning("Function dlogcdu() not implemented for copulas of class '",
            class(copula), "'; numerical differentiation used")
    logc <- function(x) dCopula(x, copula, log = TRUE)
    dim <- copula@dimension
    t(apply(u, 1, function(u.)
        grad(logc, u., side = sides(u., dim, method.args$d, 0, 1),
             method.args = method.args)))
}

## For copulas with explicit densities
dlogcduExplicitCopula <- function(copula, u, ...) {
    dim <- copula@dimension
    algNm <- paste(class(copula)[1], "pdfDerWrtArg.algr", sep=".")
    mat <- matrix(NA_real_, nrow(u), dim)
    if(exists(algNm)) {
        der.pdf.u <- get(algNm)[dim]
	alpha <- copula@parameters # typically used in eval() # JY: scalar alpha
        unames0 <- paste0("u", 1:dim)
        for (j in 1:dim) {
            unames <- unames0; unames[1] <- unames0[j]; unames[j] <- unames0[1]
            colnames(u) <- unames
            mat[,j] <- eval(der.pdf.u, data.frame(u))
        }
    } else warning("Function dlogcdu() not implemented for copulas of class '",
                   class(copula), "'")
    mat / dCopula(u, copula)
}

setMethod("dlogcdu", signature("copula"), dlogcduCopulaNum)
setMethod("dlogcdu", signature("joeCopula"), dlogcduCopulaNum)
setMethod("dlogcdu", signature("tevCopula"), dlogcduCopulaNum)
setMethod("dlogcdu", signature("archmCopula"), dlogcduExplicitCopula)
setMethod("dlogcdu", signature("evCopula"), dlogcduExplicitCopula)
setMethod("dlogcdu", signature("plackettCopula"), dlogcduExplicitCopula)

## For normalCopula objects
dlogcduNormalCopula <- function(copula, u, ...) {
    v <- qnorm(u)
    (- v %*% solve(getSigma(copula)) + v) / dnorm(v)
}

setMethod("dlogcdu", signature("normalCopula"), dlogcduNormalCopula)

## For tCopula objects
dlogcduTCopula <- function(copula, u, ...) {
    df <- getdf(copula)
    v <- qt(u,df=df)
    w <- dt(v,df=df)
    m <- v %*% solve(getSigma(copula))
    - (df + copula@dimension) * m / ((df + rowSums(m * v)) * w) +
        (df + 1) * v / ((df +  v^2) * w)
}

setMethod("dlogcdu", signature("tCopula"), dlogcduTCopula)

##################################################################################
### Partial derivatives of the log PDF wrt parameters
##################################################################################

setGeneric("dlogcdtheta", function(copula, u, ...) standardGeneric("dlogcdtheta"))

##' @title Basic implementation of dlogcdu based on numerical differentiation
##' @param copula an object of class "copula"
##' @param u evaluation point in [0,1]^d or matrix of evaluation points
##' @param method.args to be passed to grad()
##' @param ...
##' @return partial derivatives of the log copula density wrt parameters at u
##' @author Ivan Kojadinovic
dlogcdthetaCopulaNum <- function(copula, u, method.args = gradControl(d = 1e-5), ...) {
    warning("Function dlogcdtheta() not implemented for copulas of class '",
            class(copula), "'; numerical differentiation used")
    logc <- function(theta) {
        freeParam(copula) <- theta
        dCopula(u, copula, log = TRUE)
    }
    theta <- getParam(copula)
    p <- length(theta)
    lb <- attr(theta, "param.lowbnd")
    ub <- attr(theta, "param.upbnd" )
    jacobian(logc, theta, side = sides(theta, p, method.args$d, lb, ub),
             method.args = method.args)
}

## For copulas with explicit densities
dlogcdthetaExplicitCopula <- function(copula, u, ...) {
    dim <- copula@dimension
    algNm <- paste(class(copula)[1], "pdfDerWrtPar.algr", sep=".")
    if(exists(algNm) && !is.null((der.pdf.alpha <- get(algNm)[dim])[[1]])) {
	alpha <- copula@parameters # typically used in val(.)
        colnames(u) <- paste0("u", 1:dim)
        as.matrix(eval(der.pdf.alpha, data.frame(u))) / dCopula(u, copula)
    } else {
        warning("Function dlogcdtheta() not implemented for copulas of class '",
                class(copula), "'")
        matrix(NA_real_, nrow(u), nFree(copula@parameters))
    }
}

setMethod("dlogcdtheta", signature("copula"), dlogcdthetaCopulaNum)
setMethod("dlogcdtheta", signature("joeCopula"), dlogcdthetaCopulaNum)
setMethod("dlogcdtheta", signature("tevCopula"), dlogcdthetaCopulaNum)
setMethod("dlogcdtheta", signature("archmCopula"), dlogcdthetaExplicitCopula)
setMethod("dlogcdtheta", signature("plackettCopula"), dlogcdthetaExplicitCopula)
setMethod("dlogcdtheta", signature("evCopula"), dlogcdthetaExplicitCopula)

## For ellipCopula objects
dlogcdthetaEllipCopula <- function(copula, u, ...) {
    dim <- copula@dimension

    ## quantile transformation
    v <- switch(clc <- class(copula), ## FIXME: fails if copula *extends* "tCopula" (e.g.)
		"normalCopula" = qnorm(u),
		"tCopula" = {
		    df <- getdf(copula)
		    qt(u, df=df)
		},
		## else:
		stop("class ", sQuote(clc), " not implemented"))
    val <-
    if (dim == 2) {
	rho <- copula@parameters[1] # JY: single rho parameter, free
	ir2 <- 1 - rho^2
        sv2 <- rowSums(v^2) # == v[,1]^2 + v[,2]^2
	as.matrix(switch(clc,
			 "normalCopula" =
                             (rho * ir2 - rho * sv2 + (rho^2 + 1) * v[,1] * v[,2])/ir2^2,
			 "tCopula" =
                             (1 + df) * rho / -ir2 + (2 + df) * (df * rho + v[,1] * v[,2])
			 / (df * ir2 + sv2 - 2 * rho * v[,1] * v[,2])))

    } else { ##  dim >= 3
        n <- nrow(u)
        sigma <- getSigma(copula)
        detsig <- det(sigma)
        invsig <- solve(sigma)

        if (copula@dispstr %in% c("ex","ar1")) { ## exchangeable or ar1
            rho <- copula@parameters[1] # JY: single rho parameter
            dersig <- matrix(1, dim, dim)
            if (copula@dispstr == "ex") ## ex
                diag(dersig) <- 0
            else ## ar1
                for (i in 1:dim)
                    for (j in 1:dim) {
                        ij <- abs(i - j)
                        dersig[i,j] <- ij * rho^(ij - 1)
                    }

            ## MM:  sum(diag(A %*% B)) == sum(A * t(B)) .. and B=dersig is symmetric here
            ## derdetsig <- detsig * sum(diag(invsig %*% dersig))
            derdetsig <- detsig *  sum(invsig * dersig)
            derinvsig <- - invsig %*% dersig %*% invsig
            firstterm <- derdetsig/detsig

	    mat <-
		switch(clc,
		       "normalCopula" =
                           - (firstterm + rowSums((v %*% derinvsig) * v))/2,
		       "tCopula" =
                           - (firstterm + (df + dim) * rowSums((v %*% derinvsig) * v)
                    / (df +  rowSums((v %*% invsig) * v)) ) / 2)
            as.matrix(mat)
        }
	else {# unstructured or toeplitz or ...
            mat <- matrix(0, n, dim*(dim - 1)/2)
            l <- 1
            for (j in 1:(dim-1))
                for (i in (j+1):dim) {
                    derdetsig <- 2 * det(sigma[-i,-j,drop=FALSE]) * (-1)^(i+j)
                    derinvsig <- - invsig[,i] %*% t(invsig[,j]) - invsig[,j] %*% t(invsig[,i])
                    firstterm <- derdetsig/detsig

		    mat[,l] <-
			switch(clc,
			       "normalCopula" =
                                   - (firstterm + rowSums((v %*% derinvsig) * v))/2,
			       "tCopula" =
                                   - (firstterm + (df + dim) * rowSums((v %*% derinvsig) * v)
                            / (df +  rowSums((v %*% invsig) * v)) ) / 2)
                    l <- l + 1
                }
            if (copula@dispstr == "un") ## unstructured
                mat
            else if (copula@dispstr == "toep") { ## toeplitz: dim-1 parameters
                coef <- matrix(0, dim*(dim-1)/2, dim-1)
                for (k in 1:(dim-1)) {
                    m <- row(sigma) == col(sigma) + k
                    coef[,k] <- P2p(m)
                }
                mat %*% coef
            }
            else stop("Not implemented yet for the dispersion structure ", copula@dispstr)
        }
    } ## dim >= 3
    free <- isFree(copula@parameters)
    if (.hasSlot(copula, "df.fixed")) free <- free[-length(free)]
    val[, free, drop = FALSE]
}

setMethod("dlogcdtheta", signature("ellipCopula"), dlogcdthetaEllipCopula)
