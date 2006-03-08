##### kendall's tau for each copula

kendallsTau <- function(copula, ...) {
  ## bivariate association measurement
  UseMethod("kendallsTau")
}

kendallsTauEllipCopula <- function(copula) {
  rho <- copula@parameters[1]
  2 * asin(rho) /pi
}




## debye function is used for compute the assoc measure of frankCopula
debye <- function(x, k, ...) {
  Dk.integrand <- function(t) t^k / (exp(t) - 1)
  Dk.int <- function(x, k, ...) {
    integrate(Dk.integrand, 0, x, ...)$value
  }
  y <- abs(x)
  Dk <- k / y^k * sapply(y, Dk.int, k = k, ...)
  ifelse(x < 0, Dk <- Dk + k * y / (k + 1), Dk)
}


tau2paramClaytonCopula <- function(tau) {
  2 * tau / (1 - tau)
}


#### spearman's rho
spearmansRho <- function(copula, ...) {
  ## bivariate association measurement
  UseMethod("spearmansRho")
}

spearmansRhoEllipCopula <- function(copula) {
  rho <- copula@parameters[1]
  asin(rho / 2) * 6 / pi
}



#### tail dependence index
