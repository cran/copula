setClass("plackettCopula",
         representation = representation("copula"),
         contains = list("copula")
         )

plackettCopula <- function(param) {
  ## dim = 2
  dim <- 2
  val <- new("plackettCopula",
             dimension = dim,
             parameters = param[1],
             param.names = "param",
             param.lowbnd = 0,
             param.upbnd = Inf,
             message = "Plackett copula family")
  val
}
  
pplackettCopula <- function(copula, u) {
  dim <- copula@dimension
  if (is.vector(u)) u <- matrix(u, nrow = 1)
  for (i in 1:dim) assign(paste("u", i, sep=""), u[,i])
  delta <- copula@parameters[1]
  eta <- delta - 1
  ## Joe (1997, p.141)
  0.5 / eta * (1 + eta * (u1 + u2) - ((1 + eta * (u1 + u2))^2 - 4 * delta * eta * u1 * u2)^0.5)
}

dplackettCopula <- function(copula, u) {
  dim <- copula@dimension
  if (is.vector(u)) u <- matrix(u, nrow = 1)
  for (i in 1:dim) assign(paste("u", i, sep=""), u[,i])
  delta <- copula@parameters[1]
  eta <- delta - 1
  ## Joe (1997, p.141)
  ((1 + eta * (u1 + u2))^2 - 4 * delta * eta * u1 * u2)^(- 3/2) * delta * (1 + eta * (u1 + u2 - 2 * u1 * u2))
}


rplackettCopula <- function(copula, n) {
  u1 <- runif(n)
  u2 <- runif(n)
  psi <- copula@parameters[1]
  ## Johnson (1987, p.193)
  a <- u2 * (1 - u2)
  A <- psi + a * (psi - 1)^2
  B <- 2 * a * (u1 * psi^2 + 1 - u1) + psi * (1 - 2 * a)
  D <- sqrt(psi * (psi + 4 * a * u1 * (1 - u1) * (1 - psi)^2))
  v <- (B - (1 - 2 * u2) * D) / 2 / A
  cbind(u1, v)
}


setMethod("pcopula", signature("plackettCopula"), pplackettCopula)
setMethod("dcopula", signature("plackettCopula"), dplackettCopula)
setMethod("rcopula", signature("plackettCopula"), rplackettCopula)
