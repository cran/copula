########################################################################
## This script generates a spline approximation to map back and forth
## between a copula parameter and Kendall's tau for those copulas
## for which a closed form Kendall's tau formula does not exist.
########################################################################

library(copula)
source("genRcode.R")

## Taken from ../../../R/Classes.R
## numerical integration for Kendall's tau
kendallsTauCopula <- function(copula, eps = NULL, ...) {
  integrand <- function(u) pcopula(copula, u) * dcopula(copula, u)
  if (is.null(eps)) .eps <- .Machine$double.eps^0.9
  else .eps <- eps
  lower <- c(.eps, .eps)
  upper <- c(1 - .eps, 1 - .eps)
  integ <- adapt(ndim=2, lower=lower, upper=upper, functn=integrand, ...)$value
  4 * integ - 1
}

## The derivative of Tau wrt alpha
kendallsTauDerCopula <- function(copula, eps = NULL, ...) {
  pcopula.der <- function(u) {
    eval(D(copula@exprdist$cdf, "alpha"),
         list(u1 = u[1], u2 = u[2], alpha = copula@parameters))
  }
  dcopula.der <- function(u) {
    eval(D(copula@exprdist$pdf, "alpha"),
         list(u1 = u[1], u2 = u[2], alpha = copula@parameters))
  }
    
  integrand <- function(u) {
    dcopula(copula, u)^2 * pcopula.der(u) + pcopula(copula, u) * dcopula.der(u)
  }
  if (is.null(eps)) .eps <- .Machine$double.eps^0.9
  else .eps <- eps
  lower <- c(.eps, .eps)
  upper <- c(1 - .eps, 1 - .eps)
  integ <- adapt(ndim=2, lower=lower, upper=upper, functn=integrand, ...)$value
  4 * integ   
}

## plackett copula
## mycop <- plackettCopula(1) ## independence copula
## transform the parameters
## x = tanh(log(theta)), such that x is in [-1,1]
## theta = exp(atanh(x)), such that theta is in [0, Inf]
## BIG QUESTION: What is the best way to construct the grid?

## for large values of theta, use simulation
theta <- 2^(seq(3, 10, by=0.5))
kend <- theta
nsample <- 40000
for (i in 1:length(kend)) kend[i] <- cor(rcopula(plackettCopula(theta[i]), nsample), method="kendall")[1,2]

## for smaller values of theta, use numerical integration
th <- exp(atanh(seq(0.001, 0.925, len=100)))
kd <- th
for (i in 1:length(kd)) kd[i] <- kendallsTauCopula(plackettCopula(th[i]))

x <- sort(c(th, theta))
y <- sort(c(kd, kend))
tx <- tanh(log(x))

tx.plackett <- sort(c(-1, -tx, 0, tx, 1))
y.plackett <- sort(c(-1, -y, 0, y, 1))

##write.table(cbind(tx.plackett, y.plackett), row.names=FALSE, file="plackett.tau.xy")
x.plackett <- exp(atanh(tx.plackett))
x.plackett[length(x.plackett)] <- .Machine$double.xmax^.0125
write.table(cbind(x.plackett, y.plackett), row.names=FALSE, file="plackett.tau.xy")


## generate approxfun into an R file
genApproxFun("plackett.tau.xy", "plackettCopula.asso.R", "kendallsTauPlackettCopula.tr", "calibKendallsTauPlackettCopula.tr")

## for dTau.dAlpha
th <- exp(atanh(seq(0.001, 0.999, len=200)))
dTau.dAlpha <- sapply(th, function(x) kendallsTauDerCopula(plackettCopula(x)))

x <- c(1, th, 1000)
y <- c(2/9, dTau.dAlpha, 0)

write.table(cbind(x, y), row.names=FALSE, file="plackett.tau.der.xy")

genApproxFun("plackett.tau.der.xy", "plackettCopula.dTau.PosLogAlp.R", "kendallsTauDerPosLogAlpPlackettCopula.tr", NULL)
