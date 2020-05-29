## ---- message=FALSE-----------------------------------------------------------
require(copula)
source(system.file("Rsource", "AC-Liouville.R", package="copula"))
set.seed(271)

## -----------------------------------------------------------------------------
n <- 1000
theta <- 0.59
d <- 3
U <- rACsimplex(n, d=d, theta=theta, Rdist="Gamma")
cor(U, method="kendall")

## ---- fig.align="center", fig.width=6, fig.height=6---------------------------
par(pty="s")
pairs(U, gap=0, cex=0.5)

## -----------------------------------------------------------------------------
n <- 2000
theta <- 0.6
alpha <- c(1, 5, 20)
U <- rLiouville(n, alpha=alpha, theta=theta, Rdist="Gamma")
cor(U, method="kendall")

## ---- fig.align="center", fig.width=6, fig.height=6---------------------------
par(pty="s")
pairs(U, gap=0, cex=0.5)

## -----------------------------------------------------------------------------
n <- 1000
theta <- 0.59
alpha <- c(1, 3, 4)
U <- rACLiouville(n, alpha=alpha, theta=theta, family="Clayton")
cor(U, method="kendall")

## ---- fig.align="center", fig.width=6, fig.height=6---------------------------
par(pty="s")
pairs(U, gap=0, cex=0.5)

