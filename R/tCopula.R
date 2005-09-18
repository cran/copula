validTCopula <- function(object) {
  df <- object@df
  if (df <= 0) return ("df should be > 0")
  validEllipCopula(object)
}

setClass("tCopula",
         representation = representation("ellipCopula",
           df = "numeric"),
         validity = validTCopula,
         contains = list("copula", "ellipCopula")
         )


tCopula <- function(param, dim = 2, dispstr = "ex", df = 5) {
  pdim <- length(param)
  val <- new("tCopula",
             df = df,
             dispstr = dispstr,
             dimension = dim,
             parameters = param,
             param.names = paste("rho", 1:pdim, sep="."),
             param.lowbnd = rep(-1, pdim),
             param.upbnd = rep(1, pdim),
             message = "t copula family")
  val
}



rtCopula <- function(copula, n) {
  dim <- copula@dimension
  df <- copula@df
  sigma <- getSigma(copula)
  pt(rmvt(n, sigma = sigma, df = df), df = df)
}


ptCopula <- function(copula, u) {
  dim <- copula@dimension
  sigma <- getSigma(copula)
  df <- copula@df
  mycdf.vector <- function(x) {
    pmvt(lower = rep(-Inf, dim), upper = qt(x, df = df), sigma = sigma, df = df)
  }
  if (is.vector(u)) u <- matrix(u, ncol = dim)
  u[u <= 0] <- 0
  u[u >= 1] <- 1
  val <- apply(u, 1, mycdf.vector)
  val
}

dtCopula <- function(copula, u) {
  dim <- copula@dimension
  sigma <- getSigma(copula)
  df <- copula@df
  if (is.vector(u)) u <- matrix(u, ncol = dim)
  x <- qt(u, df)
  val <- dmst(x, Omega = sigma, alpha = rep(0, dim),  df = df) /
    apply(x, 1, function(v) prod(dt(v, df = df)))
  val[apply(u, 1, function(v) any(v <= 0))] <- 0
  val[apply(u, 1, function(v) any(v >= 1))] <- 0
  val
}

showTCopula <- function(object) {
  showCopula(object)
  if (object@dimension > 2) cat("dispstr: ", object@dispstr, "\n")
  cat("df: ", object@df, "\n")
}


setMethod("rcopula", signature("tCopula"), rtCopula)
setMethod("pcopula", signature("tCopula"), ptCopula)
setMethod("dcopula", signature("tCopula"), dtCopula)

setMethod("show", signature("tCopula"), showTCopula)
