library(copula)
## library(fCopulae)

#### preparation for a grid
n <- 11
u <- seq(0, 1, length=n)
v <- seq(0, 1, length=n)

#### dim = 2
umat <- as.matrix(expand.grid(u1=u, u2=u))

#### all copulas give the same tau except galambos, amh
tau <- 0.5

## frankCopula 
theta.fr <- calibKendallsTau(frankCopula(0), tau)
dcop <- dcopula(frankCopula(param=theta.fr, dim = 2), umat)
round(matrix(dcop, ncol = n), 3)

## claytonCopula
theta.cl <- calibKendallsTau(claytonCopula(1), tau)
dcop <- dcopula(claytonCopula(param=theta.cl, dim = 2), umat)
round(matrix(dcop, ncol = n), 3)

## gumbelCopula
theta.gu <- calibKendallsTau(gumbelCopula(1), tau)
dcop <- dcopula(gumbelCopula(param=theta.gu, dim = 2), umat)
round(matrix(dcop, ncol = n), 3)

## normalCopula
theta.n <- calibKendallsTau(normalCopula(0), tau)
dcop <- dcopula(normalCopula(param=theta.n, dim = 2), umat)
round(matrix(dcop, ncol = n), 3)

## tCopula
theta.t <- calibKendallsTau(tCopula(0, df=10), tau)
dcop <- dcopula(tCopula(param=theta.n, dim = 2, df=10), umat)
round(matrix(dcop, ncol = n), 3)

## galambosCopula
#theta.ga <- calibKendallsTau(galambosCopula(1), tau) ## takes too long to find
theta.ga <- 10
kendallsTau(galambosCopula(theta.ga))
dcop <- dcopula(galambosCopula(param=theta.ga), umat)
round(matrix(dcop, ncol = n), 3)

## amhCopula
#theta.amh <- calibKendallsTau(galambosCopula(1), tau) ## tau not in the range
theta.amh <- 0.9
kendallsTau(amhCopula(theta.amh))
dcop <- dcopula(amhCopula(param=theta.amh), umat)
round(matrix(dcop, ncol = n), 3)

##############################################################
################ for package fCopulae
## dcop.f <- darchmCopula(umat, alpha=2, type=1, alternative=T)
## round(matrix(dcop.f, ncol = n), 3)
## dcop.fCopulae <- darchmCopula(umat, alpha=theta.fr, type = 5, alternative = TRUE)
## round (matrix(dcop.fCopulae, ncol = n), 3)

## dim = 3, frankCopula
umat <- expand.grid(u1=u, u2=u, u3=u)
dcop <- dcopula(frankCopula(param=theta.fr, dim = 3), umat)
round(array(dcop, c(n,n,n)), 3)

## dim = 4; doesn't work well; some NaN
umat <- expand.grid(u1=u, u2=u, u3=u, u4=u)
dcop <- dcopula(frankCopula(param=theta.fr, dim = 4), umat)
round(array(dcop, c(n,n,n,n)), 3)
