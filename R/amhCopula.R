setClass("amhCopula",
         representation = representation("archmCopula"),
         contains = list("copula", "archmCopula")
         )


genFunAmh <- function(copula, u) {
  alpha <- copula@parameters[1]
  log((1 - alpha * (1 - u)) / u)
}

genInvAmh <- function(copula, s) {
  alpha <- copula@parameters[1]
  (1 - alpha) / (exp(s) - alpha)
}

genFunDer1Amh <- function(copula, u) {
  eval(amhCopula.genfun.expr[1], list(u=u, alpha=copula@parameters[1]))
}
genFunDer2Amh <- function(copula, u) {
  eval(amhCopula.genfun.expr[2], list(u=u, alpha=copula@parameters[1]))
}


amhCopula <- function(param, dim = 2) {
  ## get expressions of cdf and pdf
  cdfExpr <- function(n) {
    expr <-   "log((1 - alpha * (1 - u1)) / u1)"
    for (i in 2:n) {
      ui <- paste("u", i, sep="")
      cur <- gsub("u1", ui, expr)
      expr <- paste(expr, cur, sep=" + ")
    }
    expr <- gsub("s", expr, "(1 - alpha) / (exp(s) - alpha)")
    parse(text = expr)
  }
  
  pdfExpr <- function(cdf, n) {
    val <- cdf
    for (i in 1:n) {
      val <- D(val, paste("u", i, sep=""))
    }
    val
  }

  if (dim > 2 && param[1] < 0)
    stop("param can be negative only for dim = 2")
  cdf <- cdfExpr(dim)
  if (dim <= 6)  pdf <- pdfExpr(cdf, dim)
  else pdf <- NULL

  val <- new("amhCopula",
             dimension = dim,
             parameters = param[1],
             exprdist = c(cdf = cdf, pdf = pdf),
             param.names = "param",
             param.lowbnd = -1,
             param.upbnd = 1,
             message = "Amh copula family; Archimedean copula")
  val
}


pamhCopula <- function(copula, u) {
  dim <- copula@dimension
  if (is.vector(u)) u <- matrix(u, ncol = dim)
  alpha <- copula@parameters[1]
  if (abs(alpha) <= .Machine$double.eps^.9) return (apply(u, 1, prod))
  for (i in 1:dim) assign(paste("u", i, sep=""), u[,i])
  u1 * u2 / (1 - alpha * (1 - u1) * (1 - u2))
}


damhCopula <- function(copula, u) {
  dim <- copula@dimension
  if (is.vector(u)) u <- matrix(u, ncol = dim)
  alpha <- copula@parameters[1]
  if (abs(alpha) <= .Machine$double.eps^.9) return (rep(1, nrow(u)))
  for (i in 1:dim) assign(paste("u", i, sep=""), u[,i])
  (-1 + alpha^2*(-1 + u1 + u2 - u1*u2) - alpha*(-2 + u1 + u2 + u1*u2)) / (-1 + alpha*(-1 + u1)*(-1 + u2))^3 
}


ramhCopula <- function(copula, n) {
  warning("This function for amhCopula needs to be fixed")
  return(NULL)
  dim <- copula@dimension
  alpha <- copula@parameters[1]
  val <- matrix(runif(n * dim), nrow = n)
  if (abs(alpha) <= 100 * .Machine$double.eps)
    return (val)  ## the limit is independence
  ## Johnson (1987, p.362). Typo V_2 and p?
  u1 <- runif(n)
  u2 <- runif(n)
  b <- 1 - u1
  A <- -alpha * (2 * b * u2 + 1) + 2 * alpha^2 * b^2 * u2 + 1
  B <- alpha^2 * (4 * b^2 * u2 - 4 * b * u2 + 1) + alpha * (4 * b * u2 - 4 * b + 2) + 1
  v <- 2 * u2 * (alpha * b - 1)^2 / (A + sqrt(B))
  cbind(u1, b)
}

kendallsTauAmhCopula <- function(copula, ...) {
  alpha <- copula@parameters[1]
  ## Nelsen (1999, p.139)
  (3 * alpha - 2) / 3 / alpha - 2 / 3 * (1 - 1/alpha)^2 * log(1 - alpha)
}


setMethod("rcopula", signature("amhCopula"), ramhCopula)
setMethod("pcopula", signature("amhCopula"), pamhCopula)
setMethod("dcopula", signature("amhCopula"), damhCopula)

setMethod("genFun", signature("amhCopula"), genFunAmh)
setMethod("genInv", signature("amhCopula"), genInvAmh)
setMethod("genFunDer1", signature("amhCopula"), genFunDer1Amh)
setMethod("genFunDer2", signature("amhCopula"), genFunDer2Amh)

setMethod("kendallsTau", signature("amhCopula"), kendallsTauAmhCopula)
## setMethod("spearmansRho", signature("amhCopula"), spearmansRhoAmhCopula)
## setMethod("tailIndex", signature("amhCopula"), tailIndexAmhCopula)

## setMethod("calibKendallsTau", signature("amhCopula"), calibKendallsTauAmhCopula)

