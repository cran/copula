##### kendall's tau


tau2paramClaytonCopula <- function(tau) {
  2 * tau / (1 - tau)
}

kendallsTau <- function(copula, ...) {
  ## bivariate association measurement

}

kendallsTauEllipCopula <- function(copula) {
  rho <- copula@parameters[1]
  2 * asin(rho) /pi
}
  

kendallsTauClaytonCopula <- function(copula) {
  alpha <- copula@parameters[1]
  alpha / (alpha + 2)
}

kendallsTauGumbelCopula <- function(copula) {
  alpha <- copula@parameters[1]
  1 - 1/alpha
}


debye <- function(x, k, ...) {
  Dk.integrand <- function(t) t^k / (exp(t) - 1)
  y <- abs(x)
  Dk <- k / y^k * integrate(Dk.integrand, 0, y, ...)$value
  if (x < 0) Dk <- Dk + k * y / (k + 1)
  Dk
}

kendallsTauFrankCopula <- function(copula, ...) {
#   D1.integrand <- function(x) x / (exp(x) - 1)
#   D1fun <- function(x) {
#     y <- abs(x)
#     D1 <- integrate(D1.integrand, 0, y, ...)$value / y
#     if (x < 0) D1 <- D1 + y / 2
#     D1
#   }
  alpha <- copula@parameters[1]
  if (alpha == 0) return (0)
  - 1 + 4 / alpha * (debye(-alpha, 1) - 1)
}

#### spearman's rho

spearmansRhoFrankCopula <- function(copula, ...) {
  alpha <- copula@parameters[1]
  if (alpha == 0) return (0)
  -1 + 12/alpha * (debye(-alpha, 2) - debye(-alpha, 1))
}


#### tail dependence index
