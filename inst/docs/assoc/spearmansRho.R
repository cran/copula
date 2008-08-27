########################################################################
## This script generates a spline approximation to map back and forth
## between a copula parameter and Kendall's tau for those copulas
## for which a closed form Kendall's tau formula does not exist.
########################################################################

library(copula)
source("genRcode.R")

## Taken from ../../../R/Classes.R
## numerical integration for Spearmans's rho
spearmansRhoCopula <- function(copula, eps = NULL, ...) {
  integrand <- function(u) pcopula(copula, u)
  if (is.null(eps)) .eps <- .Machine$double.eps^0.9
  else .eps <- eps
  lower <- c(.eps, .eps)
  upper <- c(1 - .eps, 1 - .eps)
  integ <- adapt(ndim=2, lower=lower, upper=upper, functn=integrand, ...)$value
  12 * integ - 3
}

## The derivative of Rho wrt alpha
spearmansRhoDerCopula <- function(copula, eps = NULL, ...) {
  pcopula.der <- function(u) {
    eval(D(copula@exprdist$cdf, "alpha"),
         list(u1 = u[1], u2 = u[2], alpha = copula@parameters))
  }
    
  integrand <- function(u) {
    pcopula.der(u)
  }
  if (is.null(eps)) .eps <- .Machine$double.eps^0.9
  else .eps <- eps
  lower <- c(.eps, .eps)
  upper <- c(1 - .eps, 1 - .eps)
  integ <- adapt(ndim=2, lower=lower, upper=upper, functn=integrand, ...)$value
  12 * integ 
}

## BIG QUESTION: What is the best way to construct the grid?

## clayton copula
theta.neg <- c(seq(-.20, -.001, by=.001))
theta.pos <- c(seq(.001, 1, by=.001), seq(1.01, 15, by=.01))

rho.neg <- sapply(theta.neg, function(x) spearmansRho(claytonCopula(x)))
rho.pos <- sapply(theta.pos, function(x) spearmansRho(claytonCopula(x)))

rhoDer.neg <- sapply(theta.neg, function(x) spearmansRhoDerCopula(claytonCopula(x)))
rhoDer.pos <- sapply(theta.pos, function(x) spearmansRhoDerCopula(claytonCopula(x)))

x <- c(theta.neg, 0, theta.pos)
y <- c(rho.neg, 0, rho.pos)
xp <- c(theta.neg, theta.pos)
z <- c(rhoDer.neg, rhoDer.pos)

write.table(cbind(x=x, y=y), row.names=FALSE, file="clayton.rho.xy")
write.table(cbind(x=xp, y=z), row.names=FALSE, file="clayton.rhoDer.xy")

genApproxFun("clayton.rho.xy", "claytonCopula.rho.R", "spearmansRhoClaytonCopula.tr", "calibSpearmansRhoClaytonCopula.tr")
genApproxFun("clayton.rhoDer.xy", "claytonCopula.rhoDer.R", "spearmansRhoDerClaytonCopula.tr", NULL)

## gumbel copula
theta <- c(seq(1.001, 5, by=.001), seq(5.01, 15, by=.01))
rho <- sapply(theta, function(x) spearmansRho(gumbelCopula(x)))
rhoDer <- sapply(theta, function(x) spearmansRhoDerCopula(gumbelCopula(x)))

write.table(cbind(x=c(1, theta), y=c(0, rho)), row.names=FALSE, file="gumbel.rho.xy")
write.table(cbind(x=theta, y=rhoDer), row.names=FALSE, file="gumbel.rhoDer.xy")

genApproxFun("gumbel.rho.xy", "gumbelCopula.rho.R", "spearmansRhoGumbelCopula.tr", "calibSpearmansRhoGumbelCopula.tr")
genApproxFun("gumbel.rhoDer.xy", "gumbelCopula.rhoDer.R", "spearmansRhoDerGumbelCopula.tr", NULL)

