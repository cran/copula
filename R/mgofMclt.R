#################################################################################
##
##   R package Copula by Jun Yan and Ivan Kojadinovic Copyright (C) 2008
##
##   This file is part of the R package copula.
##
##   The R package copula is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package copula is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
##   You should have received a copy of the GNU General Public License
##   along with the R package copula. If not, see <http://www.gnu.org/licenses/>.
##
#################################################################################


#########################################################
## partial derivatives of the CDF wrt arguments
#########################################################

derCdfWrtArgs <- function(cop, u)
  {
    UseMethod("derCdfWrtArgs")
  }

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

setMethod("derCdfWrtArgs", signature("archmCopula"), derCdfWrtArgsExplicitCopula)
setMethod("derCdfWrtArgs", signature("plackettCopula"), derCdfWrtArgsExplicitCopula)


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
          if (p == 2)
            {
              rho <- cop@parameters
              mat[,j] <- pnorm(v[,-j], rho * v[,j], sqrt(1 - rho^2))
            }
          else
            for (i in 1:n)
              mat[i,j] <- pmvnorm(lower = rep(-Inf, p - 1),
                                  upper = v[i,-j],
                                  mean = v[i,j] * sigma[-j,j],
                                  sigma = drop(s))
        else if (class(cop) == "tCopula") 
          if (p == 2)
            {
              rho <- cop@parameters
              mat[,j] <-  pt(sqrt((df+1)/(df+v[,j]^2)) / sqrt(1 - rho^2)
                             * (v[,-j] - rho * v[,j]), df=df+1)
            }
          else
            for (i in 1:n)
              mat[i,j] <- pmvt(lower = rep(-Inf, p - 1),
                               upper = drop(sqrt((df+1)/(df+v[i,j]^2)) *
                                 (v[i,-j] - v[i,j] * sigma[-j,j])),
                               sigma = s, df = df + 1)
        
      }
    return(mat)    
  }

setMethod("derCdfWrtArgs", signature("ellipCopula"), derCdfWrtArgsEllipCopula)

#########################################################
## Plackett formula for elliptical copulas
#########################################################

plackettFormulaDim2 <- function(cop, x)
  {
    UseMethod("plackettFormulaDim2")
  }

plackettFormulaDim2NormalCopula <- function(cop, x)
  {
    rho <- cop@parameters
    return(as.matrix(exp(-(x[,1]^2 + x[,2]^2 - 2 * rho * x[,1] * x[,2])
                             /(2 * (1 - rho^2)))/(2 * pi * sqrt(1 - rho^2))))
  }

setMethod("plackettFormulaDim2", signature("normalCopula"), plackettFormulaDim2NormalCopula)

plackettFormulaDim2TCopula <- function(cop, x)
  {
    rho <- cop@parameters
    df <- cop@df
    return(as.matrix((1 + (x[,1]^2 + x[,2]^2 - 2 * rho * x[,1] * x[,2]) /
                          (df * (1 - rho^2)))^(-df / 2) / (2 * pi * sqrt(1 - rho^2))))
  }

setMethod("plackettFormulaDim2", signature("tCopula"), plackettFormulaDim2TCopula)

plackettFormula <- function(cop, p, rho, s, m, x, i, j)
  {
    UseMethod("plackettFormula")
  }

plackettFormulaNormalCopula <- function(cop, p, rho, s, m, x, i, j)
  {
    return(exp(-(x[i]^2 + x[j]^2 - 2 * rho * x[i] * x[j]) /
               (2 * (1 - rho^2))) / (2 * pi * sqrt(1 - rho^2))
           * if (p == 3) pnorm(drop((x[-c(i,j)] - m %*% x[c(i,j)])/sqrt(s)))
           else pmvnorm(lower = rep(-Inf, p - 2),
                        upper = drop(x[-c(i,j)] - m %*% x[c(i,j)]),
                        sigma = s))                 
  }

setMethod("plackettFormula", signature("normalCopula"), plackettFormulaNormalCopula)

plackettFormulaTCopula <- function(cop, p, rho, s, m, x, i, j)
  {
    df <- cop@df
    term <- 1 + (x[i]^2 + x[j]^2 - 2 * rho * x[i] * x[j]) / (df * (1 - rho^2))
    return(term^(-df / 2) / (2 * pi * sqrt(1 - rho^2)) *
           if (p == 3) pt(drop((x[-c(i,j)] - m %*% x[c(i,j)]) / sqrt(term * s)), df= df)
           else pmvt(df = df, lower = rep(-Inf, p - 2), 
                     upper = drop((x[-c(i,j)] - m %*% x[c(i,j)]) / sqrt(term)),
                     sigma = s))
  }

setMethod("plackettFormula", signature("tCopula"), plackettFormulaTCopula)


#########################################################
## partial derivatives of the CDF wrt parameters
#########################################################

derCdfWrtParams <- function(cop, u)
  {
    UseMethod("derCdfWrtParams")
  }

derCdfWrtParamsExplicitCopula <- function(cop, u)
  {
    p <- cop@dimension
    alpha <- cop@parameters
    colnames(u) <- paste("u",1:p,sep="")
    der.cdf.alpha <- get(paste(class(cop)[1], "cdfDerWrtPar.algr", sep="."))[p]
    return(as.matrix(eval(der.cdf.alpha, data.frame(u))))
  }

setMethod("derCdfWrtParams", signature("archmCopula"), derCdfWrtParamsExplicitCopula)
setMethod("derCdfWrtParams", signature("plackettCopula"), derCdfWrtParamsExplicitCopula)


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
    
#########################################################
## partial derivatives of the PDF wrt arguments
## for ellipCopula: DIVIDED BY PDF 
#########################################################

derPdfWrtArgs <- function(cop, u)
  {
    UseMethod("derPdfWrtArgs")
  }

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


#########################################################
## partial derivatives of the PDF wrt parameters
## for ellipCopula: DIVIDED BY PDF 
#########################################################

derPdfWrtParams <- function(cop, u)
  {
    UseMethod("derPdfWrtParams")
  }

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

#########################################################
## dcopula wrapper for influence coefficients
#########################################################

dcopwrap <- function(cop, u)
  {
    UseMethod("dcopwrap")
  }

dcopwrapExplicitCopula <- function(cop, u)
  {
    return(dcopula(cop,u))
  }

setMethod("dcopwrap", signature("archmCopula"), dcopwrapExplicitCopula)
setMethod("dcopwrap", signature("plackettCopula"), dcopwrapExplicitCopula)

dcopwrapEllipCopula <- function(cop, u)
  {
    return(1)
  }

setMethod("dcopwrap", signature("ellipCopula"), dcopwrapEllipCopula)


#########################################################
## influence coefficients
#########################################################

influCoef <- function(cop,u,M)
  {
    p <- cop@dimension

    ## influence: second part
    ## integrals computed from M realizations by Monte Carlo
    v <- rcopula(cop,M)
    dcop <- dcopwrap(cop,v) ## wrapper
    influ0 <- derPdfWrtParams(cop,v)/dcop
    derArg <- derPdfWrtArgs(cop,v)/dcop

    influ <- vector("list",p)
    for (i in 1:p)
        influ[[i]] <- influ0 * derArg[,i]

    ## expectation
    q <- length(cop@parameters)
    e <- crossprod(influ0)
    e <- e/M
    
    return(solve(e) %*% t(derPdfWrtParams(cop,u)/dcopwrap(cop,u) - add.influ(u,v,influ,q)))
  }

#########################################################
## second part of influence coefficients
#########################################################

add.influ <- function(u, v, influ, q)
{
  M <- nrow(v)
  p <- ncol(v)
  n <- nrow(u)
  
  o <- matrix(0,M,p)
  ob <- matrix(0,n,p)
  for (i in 1:p)
    {
      o[,i] <- order(v[,i], decreasing=TRUE)
      ob[,i] <- ecdf(v[,i])(u[,i]) * M
    }
  
  out <- matrix(0,n,q)
  for (i in 1:p)
      out <- out + rbind(rep(0,q),apply(influ[[i]][o[,i],,drop=FALSE],2,cumsum))[M + 1 - ob[,i],,drop=FALSE] / M -
        #matrix(apply(influ[[i]] * v[,i],2,mean),n,q,byrow=TRUE)
        matrix(colMeans(influ[[i]] * v[,i]),n,q,byrow=TRUE)
  return(out)
}

#########################################################
## goodness-of-fit test
#########################################################

## cop is a copula of the desired family whose parameters, if necessary, will be used
## as starting values in fitCopula

gofMCLT.PL <- function(cop, x, N, m, grid, G, R, M)
  {     
    n <- nrow(x)
    p <- ncol(x)

    ## make pseudo-observations
    u <- apply(x,2,rank)/(n+1)

    ## fit the copula
    cop <- fitCopula(u,cop,method="mpl",estimate.variance=FALSE)@copula

    ## perform R replications
    stat <- pval <- numeric(R)
    for (i in 1:R)
      { 
        
        ## grid points where to evaluate the process
        if (grid == "h0")  ## from the H0 copula 
          g <- rcopula(cop,G)
        else if (grid == "po") ## pseudo-observations
          {
            g <- u
            G <- n
          }
        else stop("Invalid grid")
        
        pcop <- pcopula(cop,g)
        
        ## compute the test statistic
        stat[i] <- .C("cramer_vonMises_2",
                      as.integer(p),
                      as.double(u),
                      as.integer(n),
                      as.double(g),
                      as.integer(G),
                      as.double(pcop),
                      stat = double(1),
                      PACKAGE="copula")$stat
        
        ## generate realizations under H0
        x0 <- rcopula(cop,m)
        
        s0 <- .C("multiplier",
                 as.integer(p),
                 as.double(x0),
                 as.integer(m),
                 as.double(g), 
                 as.integer(G), 
                 as.double(pcop), 
                 as.double(derCdfWrtArgs(cop,g)),
                 as.double(derCdfWrtParams(cop,g) %*% influCoef(cop,x0,M)),
                 as.integer(N),
                 s0 = double(N),
                 PACKAGE="copula")$s0
        
        pval[i] <- (sum(s0 >= stat[i])+0.5)/(N+1)
        
      }
    
    return(list(statistic=median(stat), pvalue=median(pval),
                sd.pvalues=sd(pval), parameters=cop@parameters))
  }

                      
