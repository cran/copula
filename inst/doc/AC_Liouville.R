## ----prelim, echo=FALSE-------------------------------------------------------
## lower resolution - less size  (default dpi = 72):
knitr::opts_chunk$set(dpi = 48)

## ----pkg+sourc, message=FALSE-------------------------------------------------
require(copula)
source(system.file("Rsource", "AC-Liouville.R", package="copula"))
set.seed(271)

## ----rACsimp------------------------------------------------------------------
n <- 1000
theta <- 0.59
d <- 3
U <- rACsimplex(n, d=d, theta=theta, Rdist="Gamma")
cor(U, method="kendall")

## ----pairs-rACsimp, fig.align="center", fig.width=6, fig.height=6-------------
par(pty="s")
pairs(U, gap=0, pch=".") # or cex=0.5

## ----Liouville----------------------------------------------------------------
n <- 2000
theta <- 0.6
alpha <- c(1, 5, 20)
U <- rLiouville(n, alpha=alpha, theta=theta, Rdist="Gamma")
cor(U, method="kendall")

## ----pairs-Liouville, fig.align="center", fig.width=6, fig.height=6-----------
par(pty="s")
pairs(U, gap=0, pch=".") # or cex=0.5

## ----ACLiou-------------------------------------------------------------------
n <- 1000
theta <- 0.59
alpha <- c(1, 3, 4)
U <- rACLiouville(n, alpha=alpha, theta=theta, family="Clayton")
cor(U, method="kendall")

## ----pairs-ACLiou, fig.align="center", fig.width=6, fig.height=6--------------
par(pty="s")
pairs(U, gap=0, pch=".") # or cex=0.5

