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


### partial derivatives of the CDF wrt arguments ###############################

setGeneric("derCdfWrtArgs", function(cop, u) standardGeneric("derCdfWrtArgs"))

## Warning: This function assumes symmetry in u
derCdfWrtArgsExplicitCopula <- function(cop, u)
  {
    p <- cop@dimension
    alpha <- cop@parameters
    der.cdf.u <- get(paste(class(cop)[1], "cdfDerWrtArg.algr", sep="."))[p]
    unames0 <- paste("u",1:p,sep="")
    mat <- matrix(0,nrow(u),p)
    for (j in 1:p)
      {
        unames <- unames0; unames[1] <- unames0[j]; unames[j] <- unames0[1]
        colnames(u) <- unames
        mat[,j] <- eval(der.cdf.u, data.frame(u))
      }
    return(mat)
  }

## this function is used for Khoudraji's device
derCdfWrtArgsIndepCopula <- function(cop, u) {
  p <- cop@dimension
  mat <- matrix(0, nrow(u), p)
  for (j in 1:p) {
    mat[,j] <- apply(u[, -j, drop=FALSE], 1, prod)
  }
  mat
}

setMethod("derCdfWrtArgs", signature("archmCopula"), derCdfWrtArgsExplicitCopula)
setMethod("derCdfWrtArgs", signature("plackettCopula"), derCdfWrtArgsExplicitCopula)
setMethod("derCdfWrtArgs", signature("evCopula"), derCdfWrtArgsExplicitCopula)
setMethod("derCdfWrtArgs", signature("gumbelCopula"), derCdfWrtArgsExplicitCopula)
setMethod("derCdfWrtArgs", signature("indepCopula"), derCdfWrtArgsIndepCopula)

derCdfWrtArgsEllipCopula <- function(cop, u)
{
    p <- cop@dimension
    sigma <- getSigma(cop)

    ## quantile transformation
    if (class(cop) == "normalCopula")
	v <- qnorm(u)
    else if (class(cop) == "tCopula")
    {
	df <- cop@df
	v <- qt(u,df=df)
    }
    else stop("not implemented")

    n <- nrow(u)
    mat <- matrix(0,n,p)

    for (j in 1:p)
    {
	s <- sigma[-j,-j] - sigma[-j,j,drop=FALSE] %*% sigma[j,-j,drop=FALSE]

	if (class(cop) == "normalCopula")
	    if (p == 2) {
		rho <- cop@parameters
		mat[,j] <- pnorm(v[,-j], rho * v[,j], sqrt(1 - rho^2))
	    }
	    else
		for (i in 1:n)
		    mat[i,j] <- pmvnorm(lower = rep(-Inf, p - 1), upper = v[i,-j],
					mean = v[i,j] * sigma[-j,j],
					sigma = drop(s))
	else if (class(cop) == "tCopula")
	    if (p == 2) {
		rho <- cop@parameters
		mat[,j] <-  pt(sqrt((df+1)/(df+v[,j]^2)) / sqrt(1 - rho^2)
			       * (v[,-j] - rho * v[,j]), df=df+1)
	    }
	    else {
		if(df != as.integer(df))
		    stop("'df' is not integer; therefore, derCdfWrtArgs() cannot be computed yet")
		for (i in 1:n)
		    mat[i,j] <- pmvt(lower = rep(-Inf, p - 1),
				     upper = drop(sqrt((df+1)/(df+v[i,j]^2)) *
				     (v[i,-j] - v[i,j] * sigma[-j,j])),
				     sigma = s, df = df + 1)
	    }

    }
    mat
}

setMethod("derCdfWrtArgs", signature("ellipCopula"), derCdfWrtArgsEllipCopula)


### Plackett formula for elliptical copulas ####################################

setGeneric("plackettFormulaDim2", function(cop, x) standardGeneric("plackettFormulaDim2"))

plackettFormulaDim2NormalCopula <- function(cop, x)
  {
    rho <- cop@parameters
    ir2 <- 1 - rho^2
    as.matrix(exp(-(x[,1]^2 + x[,2]^2 - 2 * rho * x[,1] * x[,2]) /
                  (2 * ir2)) /
              (2 * pi * sqrt(ir2)))
  }

setMethod("plackettFormulaDim2", signature("normalCopula"), plackettFormulaDim2NormalCopula)

plackettFormulaDim2TCopula <- function(cop, x)
  {
    rho <- cop@parameters
    ir2 <- 1 - rho^2
    df <- cop@df
    as.matrix((1 + (x[,1]^2 + x[,2]^2 - 2 * rho * x[,1] * x[,2]) /
               (df * ir2))^(-df / 2) / (2 * pi * sqrt(ir2)))
  }

setMethod("plackettFormulaDim2", signature("tCopula"), plackettFormulaDim2TCopula)

setGeneric("plackettFormula",  function(cop, p, rho, s, m, x, i, j) standardGeneric("plackettFormula"))

plackettFormulaNormalCopula <- function(cop, p, rho, s, m, x, i, j)
{
    exp(-(x[i]^2 + x[j]^2 - 2 * rho * x[i] * x[j]) /
               (2 * (1 - rho^2))) / (2 * pi * sqrt(1 - rho^2)) *
           (if (p == 3) pnorm(drop((x[-c(i,j)] - m %*% x[c(i,j)])/sqrt(s)))
           else pmvnorm(lower = rep(-Inf, p - 2),
                        upper = drop(x[-c(i,j)] - m %*% x[c(i,j)]),
                        sigma = s))
}

setMethod("plackettFormula", signature("normalCopula"), plackettFormulaNormalCopula)

plackettFormulaTCopula <- function(cop, p, rho, s, m, x, i, j)
{
    stopifnot(p >= 3)
    df <- cop@df
    if(df != as.integer(df) && p > 3)
	stop("'df' is not integer; therefore, plackettFormula() cannot be computed yet")
    term <- 1 + (x[i]^2 + x[j]^2 - 2 * rho * x[i] * x[j]) / (df * (1 - rho^2))
    ## return:
    term^(-df / 2) / (2 * pi * sqrt(1 - rho^2)) *
	if (p == 3) pt(drop((x[-c(i,j)] - m %*% x[c(i,j)]) / sqrt(term * s)), df= df)
	else pmvt(df = df, lower = rep(-Inf, p - 2),
		  upper = drop((x[-c(i,j)] - m %*% x[c(i,j)]) / sqrt(term)),
		  sigma = s)
}

setMethod("plackettFormula", signature("tCopula"), plackettFormulaTCopula)


### partial derivatives of the CDF wrt parameters ##############################

setGeneric("derCdfWrtParams", function(cop, u) standardGeneric("derCdfWrtParams"))


derCdfWrtParamsExplicitCopula <- function(cop, u)
  {
    p <- cop@dimension
    alpha <- cop@parameters
    colnames(u) <- paste("u",1:p,sep="")
    der.cdf.alpha <- get(paste(class(cop)[1], "cdfDerWrtPar.algr", sep="."))[p]
    return(as.matrix(eval(der.cdf.alpha, data.frame(u))))
  }

derCdfWrtParamsEvCopula <- function(cop, u) {
  alpha <- cop@parameters
  loguv <- log(u[,1], u[,2])
  w <- log(u[,2]) / loguv
  der.cdf.alpha <- pcopula(cop, u) * loguv * derAfunWrtParam(cop, w)
  return(as.matrix(der.cdf.alpha))
}

setMethod("derCdfWrtParams", signature("archmCopula"), derCdfWrtParamsExplicitCopula)
setMethod("derCdfWrtParams", signature("plackettCopula"), derCdfWrtParamsExplicitCopula)
setMethod("derCdfWrtParams", signature("evCopula"), derCdfWrtParamsExplicitCopula)
setMethod("derCdfWrtParams", signature("gumbelCopula"), derCdfWrtParamsExplicitCopula)

derCdfWrtParamsEllipCopula <- function(cop, u)
  {
    p <- cop@dimension

    ## quantile transformation
    if (class(cop) == "normalCopula")
      v <- qnorm(u)
    else if (class(cop) == "tCopula")
      v <- qt(u,df=cop@df)
    else stop("not implemented")

    if (p == 2)
      plackettFormulaDim2(cop, v)
    else
      {
        n <- nrow(u)
        sigma <- getSigma(cop)

        if (cop@dispstr %in% c("ex","ar1")) ## exchangeable or ar1
          {
            rho <- cop@parameters
            r <- matrix(c(1,-rho,-rho,1),2,2)/(1 - rho^2)
            m <- sigma[-c(1,2),c(1,2)] %*% r
            s <- sigma[-c(1,2),-c(1,2)] - sigma[-c(1,2),c(1,2)] %*% r %*% sigma[c(1,2),-c(1,2)]

            mat <- matrix(0,n,1)

            if (cop@dispstr == "ex") ## exchangeable
              for (k in 1:n)
                for (j in 1:(p-1))
                  for (i in (j+1):p)
                    mat[k,1] <- mat[k,1] + plackettFormula(cop, p, rho, s, m, v[k,], i, j)
            else ## ar1
              for (k in 1:n)
                for (j in 1:(p-1))
                  for (i in (j+1):p)
                    mat[k,1] <- mat[k,1] + (i - j) * rho ^ (i - j - 1) *
                      plackettFormula(cop, p, rho, s, m, v[k,], i, j)

            return(mat)
          }
        else # unstructured or toeplitz
          {
            mat <- matrix(0,n,p*(p-1)/2)
            l <- 1
            for (j in 1:(p-1))
              for (i in (j+1):p)
                {
                  rho <- sigma[i,j]
                  r <- matrix(c(1,-rho,-rho,1),2,2)/(1 - rho^2)
                  m <- sigma[-c(i,j),c(i,j)] %*% r
                  s <- sigma[-c(i,j),-c(i,j)] - sigma[-c(i,j),c(i,j)] %*% r %*% sigma[c(i,j),-c(i,j)]

                  for (k in 1:n)
                    mat[k,l] <- plackettFormula(cop, p, rho, s, m, v[k,], i, j)
                  l <- l + 1
                }
            if (cop@dispstr == "un") ## unstructured
              return(mat)
            else ## toeplitz: p-1 parameters
              {
                coef <- matrix(0,p*(p-1)/2,p-1)
                for (k in 1:(p-1))
                  {
                    m <- row(sigma) == col(sigma) + k
                    coef[,k] <- m[lower.tri(m)]
                  }
                return(mat %*% coef)
              }
          }
      }
  }

setMethod("derCdfWrtParams", signature("ellipCopula"), derCdfWrtParamsEllipCopula)


### Partial derivatives of the PDF wrt arguments for ellipCopula: DIVIDED BY PDF

setGeneric("derPdfWrtArgs", function(cop, u) standardGeneric("derPdfWrtArgs"))

derPdfWrtArgsExplicitCopula <- function(cop, u)
  {
    p <- cop@dimension
    alpha <- cop@parameters
    der.pdf.u <- get(paste(class(cop)[1], "pdfDerWrtArg.algr", sep="."))[p]
    unames0 <- paste("u",1:p,sep="")
    mat <- matrix(0,nrow(u),p)
    for (j in 1:p)
      {
        unames <- unames0; unames[1] <- unames0[j]; unames[j] <- unames0[1]
        colnames(u) <- unames
        mat[,j] <- eval(der.pdf.u, data.frame(u))
      }
    return(mat)
  }

setMethod("derPdfWrtArgs", signature("archmCopula"), derPdfWrtArgsExplicitCopula)
setMethod("derPdfWrtArgs", signature("plackettCopula"), derPdfWrtArgsExplicitCopula)
setMethod("derPdfWrtArgs", signature("evCopula"), derPdfWrtArgsExplicitCopula)
setMethod("derPdfWrtArgs", signature("gumbelCopula"), derPdfWrtArgsExplicitCopula)

derPdfWrtArgsNormalCopula <- function(cop, u)
  {
    v <- qnorm(u)
    return( ( - v %*% solve(getSigma(cop)) + v)/ dnorm(v) )
  }

setMethod("derPdfWrtArgs", signature("normalCopula"), derPdfWrtArgsNormalCopula)

derPdfWrtArgsTCopula <- function(cop, u)
  {

    df <- cop@df
    v <- qt(u,df=df)
    w <- dt(v,df=df)
    m <- v %*% solve(getSigma(cop))

    return( - (df + cop@dimension) * m / ((df + rowSums(m * v)) * w) +
           (df + 1) * v / ((df +  v^2) * w) )
  }

setMethod("derPdfWrtArgs", signature("tCopula"), derPdfWrtArgsTCopula)


### Partial derivatives of the PDF wrt parameters for ellipCopula: DIVIDED BY PDF

setGeneric("derPdfWrtParams", function(cop, u) standardGeneric("derPdfWrtParams"))

derPdfWrtParamsExplicitCopula <- function(cop, u)
  {
    p <- cop@dimension
    alpha <- cop@parameters
    colnames(u) <- paste("u",1:p,sep="")
    der.pdf.alpha <- get(paste(class(cop)[1], "pdfDerWrtPar.algr", sep="."))[p]
    return(as.matrix(eval(der.pdf.alpha, data.frame(u))))
  }

setMethod("derPdfWrtParams", signature("archmCopula"), derPdfWrtParamsExplicitCopula)
setMethod("derPdfWrtParams", signature("plackettCopula"), derPdfWrtParamsExplicitCopula)
setMethod("derPdfWrtParams", signature("evCopula"), derPdfWrtParamsExplicitCopula)
setMethod("derPdfWrtParams", signature("gumbelCopula"), derPdfWrtParamsExplicitCopula)

derPdfWrtParamsEllipCopula <- function(cop, u)
  {
    p <- cop@dimension

    ## quantile transformation
    if (class(cop) == "normalCopula")
      v <- qnorm(u)
    else if (class(cop) == "tCopula")
      {
        df <- cop@df
        v <- qt(u,df=df)
      }
    else stop("not implemented")

    if (p == 2)
      {
        rho <- cop@parameters
        if (class(cop) == "normalCopula")
          return(as.matrix((rho * (1 - rho^2) - rho * (v[,1]^2 + v[,2]^2) +
                            (rho^2 + 1) * v[,1] * v[,2])/(1 - rho^2)^2))
        else if (class(cop) == "tCopula")
          return(as.matrix((1 + df) * rho / (rho^2 - 1) + (2 + df) * (df * rho + v[,1] * v[,2])
                           / (df * (1 - rho^2) + v[,1]^2 + v[,2]^2 - 2 * rho * v[,1] * v[,2])))
      }
    else
      {
        n <- nrow(u)
        sigma <- getSigma(cop)
        detsig <- det(sigma)
        invsig <- solve(sigma)

        if (cop@dispstr %in% c("ex","ar1")) ## exchangeable or ar1
          {
            rho <- cop@parameters

            dersig <- matrix(1,p,p)
            if (cop@dispstr == "ex") ## ex
              diag(dersig) <- 0
            else ## ar1
              for (i in 1:p)
                for (j in 1:p)
                  dersig[i,j] <- abs(i - j) * rho^(abs(i - j) - 1)

            derdetsig <- detsig * sum(diag(invsig %*% dersig))
            derinvsig <- - invsig %*% dersig %*% invsig
            firstterm <- derdetsig/detsig

            if (class(cop) == "normalCopula")
              mat <- - (firstterm + rowSums((v %*% derinvsig) * v))/2
            else if (class(cop) == "tCopula")
              mat <- - (firstterm + (df + p) * rowSums((v %*% derinvsig) * v)
                        / (df +  rowSums((v %*% invsig) * v)) ) / 2
            return(as.matrix(mat))
          }
        else ## unstructured or toeplitz
          {
            mat <- matrix(0,n,p*(p-1)/2)
            l <- 1
            for (j in 1:(p-1))
              for (i in (j+1):p)
                {
                  derdetsig <- 2 * det(sigma[-i,-j,drop=FALSE]) * (-1)^(i+j)
                  derinvsig <- - invsig[,i] %*% t(invsig[,j]) - invsig[,j] %*% t(invsig[,i])
                  firstterm <- derdetsig/detsig
                  if (class(cop) == "normalCopula")
                    mat[,l] <- - (firstterm + rowSums((v %*% derinvsig) * v))/2
                  else if (class(cop) == "tCopula")
                    mat[,l] <- - (firstterm + (df + p) * rowSums((v %*% derinvsig) * v)
                                  / (df +  rowSums((v %*% invsig) * v)) ) / 2
                  l <- l + 1
                }
            if (cop@dispstr == "un") ## unstructured
              return(mat)
            else ## toeplitz: p-1 parameters
              {
                coef <- matrix(0,p*(p-1)/2,p-1)
                for (k in 1:(p-1))
                  {
                    m <- row(sigma) == col(sigma) + k
                    coef[,k] <- m[lower.tri(m)]
                  }
                return(mat %*% coef)
              }
          }
      }
  }

setMethod("derPdfWrtParams", signature("ellipCopula"), derPdfWrtParamsEllipCopula)


### dcopula wrapper for influence coefficients #################################

setGeneric("dcopwrap",  function(cop, u, ...) standardGeneric("dcopwrap"))

dcopwrapExplicitCopula <- function(cop, u)
  {
    return(dcopula(cop,u))
  }

setMethod("dcopwrap", signature("archmCopula"), dcopwrapExplicitCopula)
setMethod("dcopwrap", signature("plackettCopula"), dcopwrapExplicitCopula)
setMethod("dcopwrap", signature("evCopula"), dcopwrapExplicitCopula)
setMethod("dcopwrap", signature("gumbelCopula"), dcopwrapExplicitCopula)

dcopwrapEllipCopula <- function(cop, u)
  {
    return(1)
  }

setMethod("dcopwrap", signature("ellipCopula"), dcopwrapEllipCopula)
