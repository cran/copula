setClass("galambosCopula",
         representation = representation("evCopula"),
         contains = list("copula", "evCopula")
         )

AfunGalambos <- function(copula, w) {
  alpha <- copula@parameters[1]
  1 - (w^(-alpha) + (1 - w)^(-alpha))^(-1/alpha)
}

galambosCopula <- function(param) {
  ## dim = 2
  dim <- 2
  val <- new("galambosCopula",
             dimension = dim,
             parameters = param[1],
             param.names = "param",
             param.lowbnd = 0,
             param.upbnd = Inf,
             message = "Galambos copula family; Extreme value copula")
  val
}
  
pgalambosCopula <- function(copula, u) {
  dim <- copula@dimension
  if (is.vector(u)) u <- matrix(u, nrow = 1)
  for (i in 1:dim) assign(paste("u", i, sep=""), u[,i])
  alpha <- copula@parameters[1]
  eval(galambosCopula.algr$cdf)
}

dgalambosCopula <- function(copula, u) {
  dim <- copula@dimension
  if (is.vector(u)) u <- matrix(u, nrow = 1)
  for (i in 1:dim) assign(paste("u", i, sep=""), u[,i])
  alpha <- copula@parameters[1]
  eval(galambosCopula.algr$pdf)
}


rgalambosCopula <- function(copula, n) {
  u1 <- runif(n)
  v <- runif(n)
  alpha <- copula@parameters[1]
  deriv1 <- function(u1, u2) {
    eval(galambosCopula.algr$deriv1cdf)
  }
  eps <- .Machine$double.eps ^ 0.8  ## don't know a better way
  myfun <- function(u2, u1, v) {
    deriv1(u1, u2)/deriv1(u1, 1 - eps) - v
  }
  ## I don't understand the rejection method used by finmetrics yet, so
  u2 <- sapply(1:n, function(x) uniroot(myfun, c(eps, 1 - eps), v=v[x], u1=u1[x])$root)
  cbind(u1, u2)
}


setMethod("pcopula", signature("galambosCopula"), pgalambosCopula)
setMethod("dcopula", signature("galambosCopula"), dgalambosCopula)
setMethod("rcopula", signature("galambosCopula"), rgalambosCopula)

setMethod("Afun", signature("galambosCopula"), AfunGalambos)
