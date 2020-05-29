## ---- message = FALSE---------------------------------------------------------
require(copula)
require(rugarch)

## -----------------------------------------------------------------------------
## Simulate innovations
n <- 200 # sample size
d <- 2 # dimension
nu <- 3 # degrees of freedom for t
tau <- 0.5 # Kendall's tau
th <- iTau(ellipCopula("t", df = nu), tau) # corresponding parameter
cop <- ellipCopula("t", param = th, dim = d, df = nu) # define copula object
set.seed(271) # reproducibility
U <- rCopula(n, cop) # sample the copula
nu. <- 3.5 # degrees of freedom for the t margins
Z <- sqrt((nu.-2)/nu.) * qt(U, df = nu.) # margins must have mean 0 and variance 1 for ugarchpath()!

## ---- fig.align = "center", fig.width = 7.5, fig.height = 6-------------------
## Fix parameters for the marginal models
fixed.p <- list(mu  = 1,
                ar1 = 0.5,
                ma1 = 0.3,
                omega  = 2, # alpha_0 (conditional variance intercept)
                alpha1 = 0.4,
                beta1  = 0.2)
meanModel <- list(armaOrder = c(1,1))
varModel <- list(model = "sGARCH", garchOrder = c(1,1)) # standard GARCH
uspec <- ugarchspec(varModel, mean.model = meanModel,
                    fixed.pars = fixed.p) # conditional innovation density (or use, e.g., "std")

## Simulate ARMA-GARCH models using the dependent innovations
## Note: ugarchpath(): simulate from a spec; ugarchsim(): simulate from a fitted object
X <- ugarchpath(uspec,
                n.sim = n, # simulated path length
                m.sim = d, # number of paths to simulate
                custom.dist = list(name = "sample", distfit = Z)) # passing (n, d)-matrix of *standardized* innovations

## Extract the resulting series
X. <- fitted(X) # X_t = mu_t + eps_t (simulated process)
sig.X <- sigma(X) # sigma_t (conditional standard deviations)
eps.X <- X@path$residSim # epsilon_t = sigma_t * Z_t (residuals)

## Basic sanity checks :
stopifnot(all.equal(X.,    X@path$seriesSim, check.attributes = FALSE),
          all.equal(sig.X, X@path$sigmaSim,  check.attributes = FALSE),
          all.equal(eps.X, sig.X * Z,        check.attributes = FALSE))

## Plot (X_t) for each margin
matplot(X., type = "l", xlab = "t", ylab = expression(X[t]~"for each margin"))

## -----------------------------------------------------------------------------
uspec <- ugarchspec(varModel, mean.model = meanModel, distribution.model = "std")
fit <- apply(X., 2, function(x) ugarchfit(uspec, data = x))

## ---- fig.align = "center", fig.width = 6, fig.height = 6---------------------
Z. <- sapply(fit, residuals, standardize = TRUE)
U. <- pobs(Z.)
par(pty = "s")
plot(U., xlab = expression(hat(U)[1]), ylab = expression(hat(U)[2]))

## -----------------------------------------------------------------------------
fitcop <- fitCopula(ellipCopula("t", dim = 2), data = U., method = "mpl")
nu. <- rep(nu., d) # marginal degrees of freedom; for simplicity using the known ones here
est <- cbind(fitted = c(fitcop@estimate, nu.), true = c(th, nu, nu.)) # fitted vs true
rownames(est) <- c("theta", "nu (copula)", paste0("nu (margin ",1:2,")"))
est

## ---- sim-fit-----------------------------------------------------------------
set.seed(271) # reproducibility
U.. <- rCopula(n, fitcop@copula)
Z.. <- sapply(1:d, function(j) sqrt((nu.[j]-2)/nu.[j]) * qt(U..[,j], df = nu.[j]))
## => Innovations have to be standardized for ugarchsim()
sim <- lapply(1:d, function(j)
    ugarchsim(fit[[j]], n.sim = n, m.sim = 1,
              custom.dist = list(name = "sample",
                                 distfit = Z..[,j, drop = FALSE])))

## ---- fig.align = "center", fig.width = 6, fig.height = 6---------------------
X.. <- sapply(sim, function(x) fitted(x)) # simulated series X_t (= x@simulation$seriesSim)
matplot(X.., type = "l", xlab = "t", ylab = expression(X[t]))

