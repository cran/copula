#################################################################################
##
##   R package Copula by Jun Yan and Ivan Kojadinovic Copyright (C) 2009
##
##   This file is part of the R package copula.
##
##   The R package copula is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package copula is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
##   You should have received a copy of the GNU General Public License
##   along with the R package copula. If not, see <http://www.gnu.org/licenses/>.
##
#################################################################################


AfunTawn <- function(copula, w) {
  alpha <- copula@parameters[1]
  A <- alpha * w^2 - alpha * w + 1
  ifelse(w == 0 | w == 1, 1, A)
}

AfunDerTawn <- function(copula, w) {
  alpha <- copula@parameters[1]
  ## deriv(expression(alpha * w^2 - alpha * w + 1), "w", hessian=TRUE)
  value <- eval(expression({
    .value <- alpha * w^2 - alpha * w + 1
    .grad <- array(0, c(length(.value), 1L), list(NULL, c("w")))
    .hessian <- array(0, c(length(.value), 1L, 1L), list(NULL, 
        c("w"), c("w")))
    .grad[, "w"] <- alpha * (2 * w) - alpha
    .hessian[, "w", "w"] <- alpha * 2
    attr(.value, "gradient") <- .grad
    attr(.value, "hessian") <- .hessian
    .value
  }), list(alpha = alpha, w = w))
  der1 <- c(attr(value, "gradient"))
  der2 <- c(attr(value, "hessian"))
  data.frame(der1 = der1, der2 = der2)
}

tawnCopula <- function(param) {
  ## dim = 2
  dim <- 2
  ## See Table 1 from Ghoudi, Khoudraji, and Rivest (1998, CJS, in french)
  cdf <- expression( u1 * u2 * exp( - alpha * log(u1) * log(u2) / log(u1 * u2)) )
  derCdfWrtU1 <- D(cdf, "u1")
  pdf <- D(derCdfWrtU1, "u2")

  val <- new("tawnCopula",
             dimension = dim,
             exprdist = c(cdf = cdf, pdf = pdf),
             parameters = param[1],
             param.names = "param",
             param.lowbnd = 0,
             param.upbnd = 1,
             message = "Tawn copula family; Extreme value copula")
  val
}
  
ptawnCopula <- function(copula, u) {
  dim <- copula@dimension
  if (is.vector(u)) u <- matrix(u, nrow = 1)
  for (i in 1:dim) assign(paste("u", i, sep=""), u[,i])
  alpha <- copula@parameters[1]
  c(eval(tawnCopula.cdf.algr[dim]))
}

dtawnCopula <- function(copula, u) {
  dim <- copula@dimension
  if (is.vector(u)) u <- matrix(u, nrow = 1)
  for (i in 1:dim) assign(paste("u", i, sep=""), u[,i])
  alpha <- copula@parameters[1]
  c(eval(tawnCopula.pdf.algr[dim]))
}

kendallsTauTawnCopula <- function(copula) {
  alpha <- copula@parameters[1]
  ## the range of tau is [0,  0.4183992]
  8 * atan(sqrt(alpha / (4 - alpha))) / sqrt(alpha * (4 - alpha)) - 2
}

calibKendallsTauTawnCopula <- function(copula, tau) {
  alpha <- 1
  taumax <- 8 * atan(sqrt(alpha / (4 - alpha))) / sqrt(alpha * (4 - alpha)) - 2
  bad <- (tau < 0 | tau >= taumax)
  if (any(bad)) warning("tau is out of the range [0, 0.4183992]")
  ifelse(tau <= 0, 0,
         ifelse(tau >= taumax, 1, calibKendallsTauCopula(copula, tau)))
}

tauDerTawnCopula <- function(copula) {
  alpha <- copula@parameters[1]
  ##  deriv(expression( 8 * atan(sqrt(alpha / (4 - alpha))) / sqrt(alpha * (4 - alpha)) - 2), "alpha")
  value <- eval(expression({
    .expr1 <- 4 - alpha
    .expr2 <- alpha/.expr1
    .expr3 <- sqrt(.expr2)
    .expr5 <- 8 * atan(.expr3)
    .expr6 <- alpha * .expr1
    .expr7 <- sqrt(.expr6)
    .value <- .expr5/.expr7 - 2
    .grad <- array(0, c(length(.value), 1L), list(NULL, c("alpha")))
    .grad[, "alpha"] <- 8 * (0.5 * ((1/.expr1 + alpha/.expr1^2) * 
        .expr2^-0.5)/(1 + .expr3^2))/.expr7 - .expr5 * (0.5 * 
        ((.expr1 - alpha) * .expr6^-0.5))/.expr7^2
    attr(.value, "gradient") <- .grad
    .value
  }), list(alpha = alpha))
  attr(value, "gradient")
}

spearmansRhoTawnCopula <- function(copula) {
  alpha <- copula@parameters[1]
  ## from Mathematica
  ## the range of rho is [0, 0.58743682]
  integ <- ( (8 - alpha) * alpha + 8 * sqrt( (8 - alpha) * alpha ) * atan(sqrt(alpha) / sqrt(8 - alpha)) ) / ( (8 - alpha)^2 * alpha )
  if(alpha == 0) 0 else 12 * integ - 3
}

calibSpearmansRhoTawnCopula <- function(copula, rho) {
  alpha <- 1
  rhomax <- 12 * ( (8 - alpha) * alpha + 8 * sqrt( (8 - alpha) * alpha ) * atan(sqrt(alpha) / sqrt(8 - alpha)) ) / ( (8 - alpha)^2 * alpha ) - 3
  bad <- (rho < 0 | rho >= rhomax)
  if (any(bad)) warning("rho is out of the range [0, 0.58743682]")
  ifelse(rho <= 0, 0,
         ifelse(rho >= rhomax, 1, calibSpearmansRhoCopula(copula, rho)))
}

rhoDerTawnCopula <- function(copula) {
  alpha <- copula@parameters[1]
  ## deriv(expression(12 * ( (8 - alpha) * alpha + 8 * sqrt( (8 - alpha) * alpha ) * atan(sqrt(alpha) / sqrt(8 - alpha)) ) / ( (8 - alpha)^2 * alpha ) - 3), "alpha")
  value <- eval(expression({
    .expr1 <- 8 - alpha
    .expr2 <- .expr1 * alpha
    .expr4 <- 8 * sqrt(.expr2)
    .expr5 <- sqrt(alpha)
    .expr6 <- sqrt(.expr1)
    .expr7 <- .expr5/.expr6
    .expr8 <- atan(.expr7)
    .expr11 <- 12 * (.expr2 + .expr4 * .expr8)
    .expr12 <- .expr1^2
    .expr13 <- .expr12 * alpha
    .expr16 <- .expr1 - alpha
    .value <- .expr11/.expr13 - 3
    .grad <- array(0, c(length(.value), 1L), list(NULL, c("alpha")))
    .grad[, "alpha"] <- 12 * (.expr16 + (8 * (0.5 * (.expr16 * 
        .expr2^-0.5)) * .expr8 + .expr4 * ((0.5 * alpha^-0.5/.expr6 + 
        .expr5 * (0.5 * .expr1^-0.5)/.expr6^2)/(1 + .expr7^2))))/.expr13 - 
        .expr11 * (.expr12 - 2 * .expr1 * alpha)/.expr13^2
    attr(.value, "gradient") <- .grad
    .value
  }), list(alpha = alpha))
  attr(value, "gradient")  
}
###########################################################################

setMethod("pcopula", signature("tawnCopula"), ptawnCopula)
setMethod("dcopula", signature("tawnCopula"), dtawnCopula)

setMethod("Afun", signature("tawnCopula"), AfunTawn)
setMethod("AfunDer", signature("tawnCopula"), AfunDerTawn)

setMethod("kendallsTau", signature("tawnCopula"), kendallsTauTawnCopula)
setMethod("spearmansRho", signature("tawnCopula"), spearmansRhoTawnCopula)

setMethod("calibKendallsTau", signature("tawnCopula"), calibKendallsTauTawnCopula)
setMethod("calibSpearmansRho", signature("tawnCopula"), calibSpearmansRhoTawnCopula)

setMethod("tauDer", signature("tawnCopula"), tauDerTawnCopula)
setMethod("rhoDer", signature("tawnCopula"), rhoDerTawnCopula)
