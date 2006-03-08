validNormalCopula <- function(object) {
  validEllipCopula(object)
  ## can do more if needed here
}

setClass("normalCopula",
         representation = representation("ellipCopula"),
         validity = validNormalCopula,
         contains = list("copula", "ellipCopula")
         )


normalCopula <- function(param, dim = 2, dispstr = "ex") {
  pdim <- length(param)
  val <- new("normalCopula",
             dispstr = dispstr,
             dimension = dim,
             parameters = param,
             param.names = paste("rho", 1:pdim, sep="."),
             param.lowbnd = rep(-1, pdim),
             param.upbnd = rep(1, pdim),
             message = "Normal copula family",
             getRho = function(obj) {obj@parameters}
             )
  val
}



rnormalCopula <- function(copula, n) {
  dim <- copula@dimension
  sigma <- getSigma(copula)
  pnorm(rmvnorm(n, sigma = sigma))
}


pnormalCopula <- function(copula, u) {
  mycdf.vector <- function(x) {
    pmvnorm(lower = rep(-Inf, dim), upper = qnorm(x), sigma = sigma)
  }

  dim <- copula@dimension
  sigma <- getSigma(copula)
  if (is.vector(u)) u <- matrix(u, ncol = dim)
  u[u <= 0] <- 0
  u[u >= 1] <- 1
  val <- apply(u, 1, mycdf.vector)
  val
}

dnormalCopula <- function(copula, u) {
  dim <- copula@dimension
  sigma <- getSigma(copula)
  if (is.vector(u)) u <- matrix(u, ncol = dim)
  x <- qnorm(u)
  val <- dmvnorm(x, sigma = sigma) / apply(x, 1, function(v) prod(dnorm(v)))
  val[apply(u, 1, function(v) any(v <= 0))] <- 0
  val[apply(u, 1, function(v) any(v >= 1))] <- 0
  val
}

showNormalCopula <- function(object) {
  showCopula(object)
  if (object@dimension > 2) cat("dispstr: ", object@dispstr, "\n")
}

setMethod("rcopula", signature("normalCopula"), rnormalCopula)
setMethod("pcopula", signature("normalCopula"), pnormalCopula)
setMethod("dcopula", signature("normalCopula"), dnormalCopula)

setMethod("show", signature("normalCopula"), showNormalCopula)

setMethod("kendallsTau", signature("normalCopula"), kendallsTauEllipCopula)
setMethod("spearmansRho", signature("normalCopula"), spearmansRhoEllipCopula)
