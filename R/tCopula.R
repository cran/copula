validTCopula <- function(object) {
  param <- object@parameters
  df <- param[length(param)]
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
             dispstr = dispstr,
             dimension = dim,
             parameters = c(param, df),
             param.names = c(paste("rho", 1:pdim, sep="."), "df"),
             param.lowbnd = c(rep(-1, pdim), 0),
             param.upbnd = c(rep(1, pdim), Inf),
             message = "t copula family",
             getRho = function(obj) {
               param <- obj@parameters
               param[-length(param)]
             }
             )
  val
}

getdf <- function(object) {
  param <- object@parameters
  df <- param[length(param)]
  df
}

rtCopula <- function(copula, n) {
  dim <- copula@dimension
  df <- getdf(copula)
  sigma <- getSigma(copula)
  pt(rmvt(n, sigma = sigma, df = df), df = df)
}


ptCopula <- function(copula, u) {
  dim <- copula@dimension
  sigma <- getSigma(copula)
  df <- getdf(copula)
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
  df <- getdf(copula)
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
}


setMethod("rcopula", signature("tCopula"), rtCopula)
setMethod("pcopula", signature("tCopula"), ptCopula)
setMethod("dcopula", signature("tCopula"), dtCopula)

setMethod("show", signature("tCopula"), showTCopula)

setMethod("kendallsTau", signature("tCopula"), kendallsTauEllipCopula)
