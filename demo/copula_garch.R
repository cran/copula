## Copyright (C) 2012 Marius Hofert, Ivan Kojadinovic, Martin Maechler, and Jun Yan
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## The copula--GARCH model

require(copula)
require(rugarch)


### 1) Simulate data from two ARMA(1,1)-GARCH(1,1) processes with dependent innovations

## simulate innovation distribution
n <- 200 # sample size
d <- 2 # dimension
nu <- 3 # degrees of freedom for t
tau <- 0.5 # Kendall's tau
th <- iTau(ellipCopula("t", df=nu), tau) # corresponding parameter
cop <- ellipCopula("t", param=th, dim=d, df=nu) # define copula object
set.seed(271) # set seed
U <- rCopula(n, cop) # sample the copula
Z <- qnorm(U) # adjust margins

## Simulate joint ARMA(1,1)-GARCH(1,1) process with these innovations
##
## Recall: ARMA(p_1,q_1)-GARCH(p_2,q_2) model:
##         X_t = mu_t + sigma_t * Z_t
##        mu_t = mu + \sum_{k=1}^{p_1} \phi_k  * (X_{t-k}-\mu) +
##                    \sum_{k=1}^{q_1} \theta_k* (X_{t-k}-\mu_{t-k})
##   sigma_t^2 = \alpha_0 + \sum_{k=1}^{p_2} \alpha_k* (X_{t-k}-\mu_{t-k})^2 +
##                          \sum_{k=1}^{q_2} \beta_k * sigma_{t-k}^2
## NB: alternatively use   X_t - mu_t = sigma_t * Z_t   in last two eq.
##
## set parameters
fixed.p <- list(mu  = 1,
                ar1 = 0.5,
                ma1 = 0.3,
                omega = 2, # alpha_0 (conditional variance intercept)
                alpha1= 0.4,
                beta1 = 0.2)
varModel <- list(model = "sGARCH", garchOrder=c(1,1)) # standard GARCH
uspec <- ugarchspec(varModel, mean.model = list(armaOrder=c(1,1)),
                    fixed.pars = fixed.p,
                    distribution.model = "norm") # conditional innovation density
## note: ugarchpath(): simulate from a spec; ugarchsim(): simulate from a fitted object
X <- ugarchpath(uspec,
                n.sim= n, # simulated path length
                m.sim= d, # number of paths to simulate
                custom.dist=list(name="sample", distfit=Z)) # passing sample (n x d)-matrix
str(X, max.level=2) # => @path$sigmaSim, $seriesSim, $residSim
matplot(X@path$sigmaSim,  type="l") # plot of sigma's (conditional standard deviations)
matplot(X@path$seriesSim, type="l") # plot of X's
matplot(X@path$residSim,  type="l") # plot of Z's
plot(pobs(X@path$residSim)) # plot of Z's pseudo-observations => seem fine


### 2) Fit procedure based on the simulated data ###############################

## fit ARMA(1,1)-GARCH(1,1) process to X
## remove 'fixed.pars' from specification to be able to fit
uspec <- ugarchspec(varModel, mean.model = list(armaOrder=c(1,1)),
                    distribution.model = "norm")
fit <- apply(X@path$seriesSim, 2, function(x) ugarchfit(uspec, x))
str(fit, max.level=3)
str(fit[[1]], max.level=2) # for first time series
stopifnot(identical(fit[[1]]@fit$residuals, residuals(fit[[1]]@fit))) # => the same

## check residuals
Z. <- sapply(fit, function(fit.) residuals(fit.@fit))
U. <- pobs(Z.)
plot(U., # plot of Z's pseudo-observations => seem fine
     xlab=expression(italic(hat(U)[1])),
     ylab=expression(italic(hat(U)[2])))

## fit a t copula to the residuals Z
fitcop <- fitCopula(ellipCopula("t", dim=2), data=U., method="mpl")
rbind(est = fitcop@estimate, true = c(th, nu)) # hat{rho}, hat{nu}; close to th, nu


### 3) Simulate from the fitted model ##########################################

## simulate from the fitted copula model
U.. <- rCopula(n, fitcop@copula)
Z.. <- qnorm(U..)

## simulate from the fitted time series model
X..sim <- lapply(1:d, function(j)
                 ugarchsim(fit[[j]], n.sim=n, m.sim=1,
                           custom.dist=list(name="sample",
                           distfit=Z..[,j, drop=FALSE]))@simulation)
str(X..sim, max.level=3)
X..Z <- sapply(X..sim, `[[`, "residSim")
X.. <- sapply(X..sim, `[[`, "seriesSim")

plot(X..Z, main="residSim"); abline(h=0,v=0, lty=2, col=adjustcolor(1, .5))
plot(X.., main="seriesSim"); abline(h=0,v=0, lty=2, col=adjustcolor(1, .5))
matplot(pobs(X..), type="l")
plot(pobs(X..), main="pobs(series..)"); rect(0,0,1,1, border=adjustcolor(1, 1/2))
