setClass("gumbelCopula",
         representation = representation("archmCopula"),
         contains = list("copula", "archmCopula")
         )

#### genFun and related functions
genFunGumbel <- function(copula, u) {
  alpha <- copula@parameters[1]
  ( - log(u))^alpha
}

genInvGumbel <- function(copula, s) {
  alpha <- copula@parameters[1]
  exp( -s^(1 / alpha) )
}

genDer1Gumbel <- function(copula, u) {
  alpha <- copula@parameters[1]
  -((-log(u))^(alpha - 1) * (alpha * (1/u))) 
}


genDer2Gumbel <- function(copula, u) {
  alpha <- copula@parameters[1]
  (-log(u))^(alpha - 1) * (alpha * (1/u^2)) + (-log(u))^((alpha - 1) - 1) * ((alpha - 1) * (1/u)) * (alpha * (1/u))
}

gumbelCopula <- function(param, dim = 2) {
  ## get expressions of cdf and pdf
  cdfExpr <- function(n) {
    expr <- "( - log(u1))^alpha"
    for (i in 2:n) {
      cur <- paste( "(-log(u", i, "))^alpha", sep = "")
      expr <- paste(expr, cur, sep=" + ")
    }
    expr <- paste("exp(- (", expr, ")^ (1/alpha))")
    parse(text = expr)
  }
  
  pdfExpr <- function(cdf, n) {
    val <- cdf
    for (i in 1:n) {
      val <- D(val, paste("u", i, sep=""))
    }
    val
  }

  cdf <- cdfExpr(dim)
  pdf <- pdfExpr(cdf, dim)
  val <- new("gumbelCopula",
             dimension = dim,
             parameters = param[1],
             exprdist = c(cdf = cdf, pdf = pdf),
             param.names = "param",
             param.lowbnd = 1,
             param.upbnd = Inf,
             message = "Gumbel copula family; Archimedean copula; Extreme value copula")
  val
}


rgumbelCopula <- function(copula, n) {
  ## frailty is stable(1,0,0) with 1/alpha
  dim <- copula@dimension
  alpha <- copula@parameters[1]
  b <- 1/alpha
  ## stable (b, 1), 0 < b < 1, Chambers, Mallows, and Stuck 1976, JASA, p.341
  fr <- rPosStable(n, b)
  fr <- matrix(fr, nrow=n, ncol=dim)
  ## now gumbel copula
  val <- matrix(runif(dim * n), nrow = n)
  genInv(copula, - log(val) / fr)
}


pgumbelCopula <- function(copula, u) {
  if (is.vector(u)) u <- matrix(u, ncol = dim)
  cdf <- copula@exprdist$cdf
  dim <- copula@dimension
  for (i in 1:dim) assign(paste("u", i, sep=""), u[,i])
  alpha <- copula@parameters[1]
  eval(cdf)
}

dgumbelCopula <- function(copula, u) {
  if (is.vector(u)) u <- matrix(u, nrow = 1)
  pdf <- copula@exprdist$pdf
  dim <- copula@dimension
  for (i in 1:dim) assign(paste("u", i, sep=""), u[,i])
  alpha <- copula@parameters[1]
  eval(pdf)
}


kendallsTauGumbelCopula <- function(copula) {
  alpha <- copula@parameters[1]
  1 - 1/alpha
}

setMethod("rcopula", signature("gumbelCopula"), rgumbelCopula)
setMethod("pcopula", signature("gumbelCopula"), pgumbelCopula)
setMethod("dcopula", signature("gumbelCopula"), dgumbelCopula)

setMethod("genFun", signature("gumbelCopula"), genFunGumbel)
setMethod("genInv", signature("gumbelCopula"), genInvGumbel)

setMethod("genDer1", signature("gumbelCopula"), genDer1Gumbel)
setMethod("genDer2", signature("gumbelCopula"), genDer2Gumbel)

setMethod("kendallsTau", signature("gumbelCopula"), kendallsTauGumbelCopula)

