## ----prelim, echo=FALSE-------------------------------------------------------
## lower resolution - less size  (default dpi = 72):
knitr::opts_chunk$set(dpi = 48)

## ----message = FALSE----------------------------------------------------------
library(copula) # for hierarchical frailties, generators etc.
library(mev) # for rmev()

## -----------------------------------------------------------------------------
## Parameters
d. <- c(2, 3) # copula sector dimensions
d <- sum(d.) # copula dimension
stopifnot(d >= 3) # for 3d plots
n <- 777 # sample size  (was 1000; save space for *.html)

## -----------------------------------------------------------------------------
##' @title Generated two hierarchical frailties (V0, V01, V02)
##' @param n sample size
##' @param family copula family
##' @param tau Kendall's tau
##' @return 3-column matrix containing the hierarchical frailties
##' @author Marius Hofert and Avinash Prasad
rV012 <- function(n, family, tau)
{
    stopifnot(n >= 1, is.character(family), length(tau) == 3, tau > 0)
    cop <- getAcop(family) # corresponding AC
    th <- iTau(cop, tau = tau) # copula parameters
    V0 <- cop@V0(n, theta = th[1]) # sample frailties V_0
    V01 <- cop@V01(V0, theta0 = th[1], theta1 = th[2]) # sample frailties V_01
    V02 <- cop@V01(V0, theta0 = th[1], theta1 = th[3]) # sample frailties V_02
    cbind(V0  = V0, V01 = V01, V02 = V02)
}

## -----------------------------------------------------------------------------
## Plot with possible export to PDF
mypairs <- function(x, pch = ".", file = NULL, width = 6, height = 6, ...)
{
    opar <- par(pty = "s")
    on.exit(par(opar))
    doPDF <- !is.null(file)
    if(doPDF) {
        pdf(file = file, width = width, height = height)
        stopifnot(require(crop))
    }
    pairs2(x, pch = pch, ...)
    if(doPDF) dev.off.crop(file)
}

## ----fig.align = "center", fig.width = 6, fig.height = 6----------------------
## Sample frailties (to recycle the frailties, we apply the Marshall--Olkin (MO)
## algorithm manually here)
tau <- 0.4 # Kendall's tau
family <- "Clayton" # frailty family
cop <- getAcop(family) # corresponding AC
th <- iTau(cop, tau = tau) # copula parameter
set.seed(271) # for reproducibility
V <- cop@V0(n, theta = th) # sample frailty
E <- matrix(rexp(n * d), ncol = d) # sample Exp(1)
U.AC <- cop@psi(E/V, theta = th) # MO

## Plot
mypairs(U.AC)

## ----fig.align = "center", fig.width = 6, fig.height = 6----------------------
## Generate Gumbel EVC sample with Exp(1) margins
tau.EVC <- 0.5 # Kendall's tau
family.EVC <- "Gumbel" # EVC family
th.EVC <- iTau(getAcop(family.EVC), tau = tau.EVC) # EVC parameter
cop.EVC <- onacopulaL(family.EVC, list(th.EVC, 1:d)) # EVC
set.seed(271) # for reproducibility
U.EVC <- rCopula(n, copula = cop.EVC) # sample EVC
E.EVC <- -log(U.EVC) # map the EVC to Exp(1) margins

## Combine with the (same) Clayton frailties (as before)
U.AXC <- cop@psi(E.EVC/V, theta = th)

## Plot
mypairs(U.AXC)

## ----fig.align = "center", fig.width = 6, fig.height = 6----------------------
## Generate hierarchical frailties
tau.N <- c(0.2, 0.4, 0.6) # Kendall's taus
family.N <- "Clayton" # frailty family
cop.N <- getAcop(family.N) # corresponding AC
th.N <- iTau(cop.N, tau = tau.N) # copula parameters
set.seed(271) # for reproducibility
V.N <- rV012(n, family = family.N, tau = tau.N) # generate corresponding frailties
V0  <- V.N[,"V0"]
V01 <- V.N[,"V01"]
V02 <- V.N[,"V02"]

## Combine with independent Exp(1)
U.NAC <- cbind(cop.N@psi(E[,1:d.[1]]/V01,     theta = th.N[2]),
               cop.N@psi(E[,(d.[1]+1):d]/V02, theta = th.N[3]))

## Plot
mypairs(U.NAC)

## ----fig.align = "center", fig.width = 6, fig.height = 6----------------------
## Generate samples
U.HAXC <- cbind(cop.N@psi(E.EVC[,1:d.[1]]/V01,     theta = th.N[2]),
                cop.N@psi(E.EVC[,(d.[1]+1):d]/V02, theta = th.N[3]))

## Plot
mypairs(U.HAXC)

## ----fig.align = "center", fig.width = 6, fig.height = 6----------------------
## Sampling from a hierarchical Gumbel copula
tau.HEVC <- c(0.2, 0.5, 0.7) # Kendall's taus
family.HEVC <- "Gumbel" # EVC family
cop.HEVC <- getAcop(family.HEVC) # corresponding base EVC
th.HEVC <- iTau(cop.HEVC, tau = tau.HEVC) # copula parameters
cop.HEVC <- onacopulaL(family.HEVC, list(th.HEVC[1], NULL,
                                     list(list(th.HEVC[2], 1:d.[1]),
                                          list(th.HEVC[3], (d.[1]+1):d)))) # hierarchical EVC
set.seed(271) # for reproducibility
U.HEVC <- rCopula(n, copula = cop.HEVC) # sample the HEVC
E.HEVC <- -log(U.HEVC) # map the HEVC samples to Exp(1) margins

## Combine the hierarchical Gumbel EVC with Exp(1) margins with the hierarchical frailties
U.HAXC.HEVC.same <- cbind(cop.N@psi(E.HEVC[,1:d.[1]]/V01,     theta = th.N[2]),
                          cop.N@psi(E.HEVC[,(d.[1]+1):d]/V02, theta = th.N[3]))

## Plot
mypairs(U.HAXC.HEVC.same)

## ----fig.align = "center", fig.width = 6, fig.height = 6----------------------
## Combine the hierarchical Gumbel EVC with Exp(1) margins with the hierarchical frailties
## in an 'asymmetric' way (hierarchical structures of the HEVC and the hierarchical frailties differ)
d.. <- rev(d.) # 3, 2
U.HAXC.HEVC.dffr <- cbind(cop.N@psi(E.HEVC[,1:d..[1]]/V01,     theta = th.N[2]), # first three components get V01
                          cop.N@psi(E.HEVC[,(d..[1]+1):d]/V02, theta = th.N[3])) # last two components get V02

## Plot
mypairs(U.HAXC.HEVC.dffr)

## ----fig.align = "center", fig.width = 6, fig.height = 6----------------------
## Correlation matrix (parameter matrix of the extremal t)
P <- matrix(0.7, ncol = d, nrow = d)
diag(P) <- 1

## Extremal t copula (EVC)
nu <- 3.5
set.seed(271)
U.EVC <- exp(-1/rmev(n, d = d, param = nu, sigma = P, model = "xstud")) # apply unit Frechet margins to get U

## Plot
mypairs(U.EVC)

## ----fig.align = "center", fig.width = 6, fig.height = 6----------------------
## Construct a hierarchical correlation matrix for the extremal t
P.h <- matrix(0.2, ncol = d, nrow = d)
P.h[1:d.[1], 1:d.[1]] <- 0.5
P.h[(d-d.[2]+1):d, (d-d.[2]+1):d] <- 0.7
diag(P.h) <- 1

## Hierarchical extremal t copula (HEVC)
X.et <- rmev(n, d = d, param = nu, sigma = P.h, model = "xstud")
U.HEVC <- exp(-1/X.et) # apply unit Frechet margins to get U

## Plot
mypairs(U.HEVC)

## ----fig.align = "center", fig.width = 6, fig.height = 6----------------------
## Construct samples
E.HEVC <- -log(U.HEVC) # map to Exp(1) margins
U.AXC.HEVC <- cop@psi(E.HEVC/V, theta = th)

## Plot
mypairs(U.AXC.HEVC)

## ----fig.align = "center", fig.width = 6, fig.height = 6----------------------
## Construct samples
E.EVC <- -log(U.EVC)
U.HAXC.EVC <- cbind(cop.N@psi(E.EVC[,1:d.[1]]/V01,     theta = th.N[2]),
                    cop.N@psi(E.EVC[,(d.[1]+1):d]/V02, theta = th.N[3]))

## Plot
mypairs(U.HAXC.EVC)

## ----fig.align = "center", fig.width = 6, fig.height = 6----------------------
## Construct samples
U.NAXC.HEVC.same <- cbind(cop.N@psi(E.HEVC[,1:d.[1]]/V01,     theta = th.N[2]),
                          cop.N@psi(E.HEVC[,(d.[1]+1):d]/V02, theta = th.N[3]))

## Plot
mypairs(U.HAXC.HEVC.same)

## ----fig.align = "center", fig.width = 6, fig.height = 6----------------------
## Construct samples
U.NAXC.HEVC.dffr <- cbind(cop.N@psi(E.HEVC[,1:d..[1]]/V01,     theta = th.N[2]), # first three components get V01
                          cop.N@psi(E.HEVC[,(d..[1]+1):d]/V02, theta = th.N[3])) # last two components get V02

## Plot
mypairs(U.HAXC.HEVC.dffr)

