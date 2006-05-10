setClass("frankCopula",
         representation = representation("archmCopula"),
#         validity = validCopula,
         contains = list("copula", "archmCopula")
         )


genFunFrank <- function(copula, u) {
  alpha <- copula@parameters[1]
  - log( (exp(- alpha * u) - 1) / (exp(- alpha) - 1))
}

genInvFrank <- function(copula, s) {
  alpha <- copula@parameters[1]
  -1/alpha * log(1 + exp(-s) * (exp(-alpha) - 1))
}

frankCopula <- function(param, dim = 2) {
  ## get expressions of cdf and pdf
  cdfExpr <- function(n) {
    expr <-   "- log( (exp(- alpha * u1) - 1) / (exp(- alpha) - 1))"
    for (i in 2:n) {
      cur <- paste("- log( (exp(- alpha * u", i, ") - 1) / (exp(- alpha) - 1))", sep="")
      expr <- paste(expr, cur, sep=" + ")
    }
    expr <- paste("-1/alpha * log(1 + exp(-(", expr, ")) * (exp(-alpha) - 1))")
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
  val <- new("frankCopula",
             dimension = dim,
             parameters = param[1],
             exprdist = c(cdf = cdf, pdf = pdf),
             param.names = "param",
             param.lowbnd = -Inf,
             param.upbnd = Inf,
             message = "Frank copula family; Archimedean copula")
  val
}

rfrankBivCopula <- function(copula, n) {
  val <- cbind(runif(n), runif(n))
  alpha <- copula@parameters[1]
  val[,2] <- -1/alpha * log(1 + val[,2] * (1 - exp(-alpha)) / (exp(-alpha * val[,1]) * (val[,2] - 1) - val[,2]))
  val
}

rfrankCopula <- function(copula, n) {
  dim <- copula@dimension
  alpha <- copula@parameters[1]
  if (abs(alpha) <= 100 * .Machine$double.eps)
    return (matrix(runif(n * dim), nrow = n))
  if (dim == 2) return (rfrankBivCopula(copula, n))
  ## the frailty is a log series distribution with a = 1 - exp(-alpha)
  fr <- rlogseries(n, 1 - exp(-alpha))
  fr <- matrix(fr, nrow = n, ncol = dim)
  val <- matrix(runif(dim * n), nrow = n)
  genInv(copula, - log(val) / fr)
}


pfrankCopula <- function(copula, u) {
  if (is.vector(u)) u <- matrix(u, ncol = dim)
  cdf <- copula@exprdist$cdf
  dim <- copula@dimension
  for (i in 1:dim) assign(paste("u", i, sep=""), u[,i])
  alpha <- copula@parameters[1]
  eval(cdf)
}

dfrankCopula <- function(copula, u) {
  if (is.vector(u)) u <- matrix(u, nrow = 1)
  pdf <- copula@exprdist$pdf
  dim <- copula@dimension
  for (i in 1:dim) assign(paste("u", i, sep=""), u[,i])
  alpha <- copula@parameters[1]
  val <- eval(pdf)
  val[apply(u, 1, function(v) any(v <= 0))] <- 0
  val[apply(u, 1, function(v) any(v >= 1))] <- 0
  val
}


kendallsTauFrankCopula <- function(copula, ...) {
  alpha <- copula@parameters[1]
  if (alpha == 0) return (0)
  - 1 + 4 / alpha * (debye(-alpha, 1) - 1)
}

spearmansRhoFrankCopula <- function(copula, ...) {
  alpha <- copula@parameters[1]
  if (alpha == 0) return (0)
  -1 + 12/alpha * (debye(-alpha, 2) - debye(-alpha, 1))
}

setMethod("rcopula", signature("frankCopula"), rfrankCopula)
setMethod("pcopula", signature("frankCopula"), pfrankCopula)
setMethod("dcopula", signature("frankCopula"), dfrankCopula)

setMethod("genFun", signature("frankCopula"), genFunFrank)
setMethod("genInv", signature("frankCopula"), genInvFrank)

setMethod("kendallsTau", signature("frankCopula"), kendallsTauFrankCopula)
setMethod("spearmansRho", signature("frankCopula"), spearmansRhoFrankCopula)

