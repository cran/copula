## ----message=FALSE------------------------------------------------------------
require(gsl) # for exponential integral
require(copula)
doPDF <- FALSE

## -----------------------------------------------------------------------------
## Tail integral of a variance gamma Levy process
## \bar{\nu}(x) = \int_x^\infty f(z) dz for the Levy density
## f(z) = (c/x)*exp(-lambda*x) for x>0, c=1/kappa and
## lambda=(sqrt(theta^2+2*sigma^2/kappa)-theta)/sigma^2
nu_bar_vargamma <- function(x, th, kap, sig) {
    lambda <- (sqrt(th^2+2*sig^2/kap)-th)/sig^2
    -expint_Ei(-lambda*x, give=FALSE)/kap
}

## -----------------------------------------------------------------------------
## Inverse of the tail integral of a variance gamma Levy process
## \bar{\nu}(x) = \int_x^\infty f(z) dz for the Levy density
## f(z) = (c/x)*exp(-lambda*x) for x>0, c=1/kappa and
## lambda=(sqrt(theta^2+2*sigma^2/kappa)-theta)/sigma^2
nu_bar_inv_vargamma <- function(Gamma, th, kap, sig, ...)
{
    max.val <- nu_bar_vargamma(.Machine$double.xmin, th=th, kap=kap, sig=sig)
    res <- numeric(length(Gamma))
    large <- Gamma >= max.val
    res[large] <- 0 # de facto indistinguishable from 0 anyways
    if(any(!large)) {
        lambda <- (sqrt(th^2+2*sig^2/kap)-th)/sig^2
        nu_bar_vargamma_minus <- function(x, z)
            -expint_Ei(-lambda*x, give=FALSE)/kap - z
        res[!large] <- vapply(Gamma[!large], function(Gamma.)
            uniroot(nu_bar_vargamma_minus, z=Gamma.,
                    interval=c(.Machine$double.xmin, 29), ...)$root, NA_real_)
    }
    res
}

## -----------------------------------------------------------------------------
## Transforming Gamma with variance-gamma Levy margins
hom_vargamma_Levy <- function(Gamma, th, kap, sig)
{
    U <- runif(nrow(Gamma)) # jump times
    ord <- order(U) # determine the order of the U's
    jump_time <- U[ord] # (sorted) jump times
    jump_size <- apply(Gamma, 2, function(y)
        nu_bar_inv_vargamma(y, th=th, kap=kap, sig=sig)) # (unsorted) jump sizes (apply inverses of marginal tail integrals)
    value <- apply(jump_size, 2, function(x) cumsum(x[ord])) # sort jump sizes according to U's and add them up => (L_t) at jump times
    list(jump_time=jump_time, value=value)
}

## -----------------------------------------------------------------------------
## \bar{\psi} for Clayton Levy copulas
psi_bar_Clayton <- function(t, theta) t^(-1/theta)

## -----------------------------------------------------------------------------
## V_{01} for nested Clayton Levy copulas
## Note: V_{01,k} | V_{0,k} ~ LS^{-1}[\bar{\psi}_{01}(.; V_{0,k})] with
##       \bar{\psi}_{01}(t; V_{0,k}) = \exp(-V_{0,k} t^{\theta_0/\theta_1})
##       = copGumbel@V01() (not copClayton@V01()!)
V01_nested_Clayton_Levy <- function(V0, theta0, theta1)
    copGumbel@V01(V0, theta0=theta0, theta1=theta1)

## -----------------------------------------------------------------------------
## Generate Gamma for a d-dimensional Clayton Levy copula with parameter theta
## Note: - Don't confuse the Clayton parameter theta with the parameter th
##         for the marginal tail integral (variance gamma)
##       - The advantage of a fixed truncation point Gamma^* is that one can
##         correct the bias introduced when cutting off small jumps by adding
##         a drift; see Asmussen, Rosinski (2001) for more details
##       - The best stopping criterion would be if we are sure that in each
##         dimension all generated Gammas which are <= Gamma^* (= Gamma.star)
##         form a sample of jump times of a homogeneous Poi(1) process on
##         [0, Gamma^*]; this could be tested.
##       - We go with a simpler stopping criterion here: Given a burn.in value
##         (integer), we stop (only) if in the last burn.in-many generated
##         Gammas each had at least one component larger than Gamma^*. So it's
##         unlikely that we still get such (uniformly) small Gammas
##         (<= Gamma^* in each component); large Gamma => small jump
##         => we correctly only truncate (small) jumps.
Gamma_Clayton_Levy <- function(d, theta, Gamma.star, burn.in)
{
    stopifnot(d >= 1, length(theta) == 1, theta > 0 , Gamma.star > 0,
              burn.in >= 1)
    Gamma <- matrix(, nrow=0, ncol=d)
    count <- 0
    W <- 0
    repeat {
        E <- rexp(d+1)
        W <- W + E[d+1]
        V <- (W/theta * gamma(1/theta))^theta # generate V = F^{-1}(W)
        G <- psi_bar_Clayton(E[1:d]/V, theta=theta) # Gamma
        Gamma <- rbind(Gamma, G) # update Gamma
        if(count >= burn.in) break # stopping criterion
        count <- if(any(G <= Gamma.star)) 1 else count + 1 # if there are still Gammas <= Gamma^*, keep generating Gammas
    }
    Gamma[Gamma > Gamma.star] <- Inf # => produce \bar{\mu}(.) = 0 (0-height jumps)
    Gamma
}

## -----------------------------------------------------------------------------
## Generate Gamma for a 4-dimensional nested Clayton Levy copula
Gamma_nested_Clayton_Levy <- function(theta, Gamma.star, burn.in)
{
    stopifnot(d >= 1, length(theta) == 3, theta > 0, min(theta[2:3]) >= theta[1],
              Gamma.star > 0, burn.in >= 1)
    d <- 4 # d must be 4 here; obviously, this could be generalized
    Gamma <- matrix(, nrow=0, ncol=d)
    count<- 0
    W <- 0
    repeat {
        E <- rexp(d+1)
        W <- W + E[d+1]
        V0 <- (W/theta[1] * gamma(1/theta[1]))^theta[1] # generate V_0 = F_0^{-1}(W)
        V01 <- V01_nested_Clayton_Levy(V0, theta0=theta[1], theta1=theta[2]) # generate V_{01}
        V02 <- V01_nested_Clayton_Levy(V0, theta0=theta[1], theta1=theta[3]) # generate V_{02}
        G <- c(psi_bar_Clayton(E[1:2]/V01, theta=theta[2]),
               psi_bar_Clayton(E[3:4]/V02, theta=theta[3])) # Gamma
        Gamma <- rbind(Gamma, G) # update Gamma
        if(count >= burn.in) break # stopping criterion
        count <- if(any(G <= Gamma.star)) 1 else count + 1 # if there are still Gammas <= Gamma^*, keep generating Gammas
    }
    Gamma[Gamma > Gamma.star] <- Inf # => produce \bar{\mu}(.) = 0 (0-height jumps)
    Gamma
}

## -----------------------------------------------------------------------------
## Plot Gammas
plot_Gamma <- function(Gamma, Gamma.star, file=NULL, ...)
{
    stopifnot(is.matrix(Gamma), (d <- ncol(Gamma)) >= 2,
              is.null(file) || is.character(file))
    palette <- colorRampPalette(c("black", "royalblue3", "darkorange2",
                                  "maroon3"), space="Lab")
    cols <- palette(d)
    ## cols <- adjustcolor(cols, alpha.f=0.1) # no improvement here
    doPDF <- !is.null(file)
    if(doPDF) pdf(file=file, width=7, height=7)
    plot(Gamma[,1], type="l", ylim=range(Gamma, finite=TRUE), # omit Inf
         log="y", xlab="k", ylab="", col=cols[1], ...)
    for(j in 2:d) lines(Gamma[,j], col=cols[j])
    abline(h=Gamma.star, lty=2, lwd=1.6)
    legend("bottomright", bty="n", lty=c(rep(1, d), 2), lwd=c(rep(1,d), 1.6),
           col=c(cols, "black"), as.expression( c(lapply(1:d, function(j)
               bquote(Gamma[list(k,.(j))])), list(bquote(Gamma*"*")))))
    if(doPDF) dev.off()
}

## -----------------------------------------------------------------------------
## Plot a multivariate Levy process
plot_Levy <- function(L, file=NULL, ...)
{
    stopifnot(is.matrix(L$value), (d <- ncol(L$value)) >= 2,
              length(L$jump_time)==nrow(L$value),
              is.null(file) || is.character(file))
    palette <- colorRampPalette(c("black", "royalblue3", "darkorange2",
                                  "maroon3"), space="Lab")
    cols <- palette(d) # d colors
    x_jump_time <- c(0, L$jump_time, 1) # extended jump times (for nicer plotting)
    x_L <- rbind(rep(0, d), L$value, L$value[nrow(L$value),]) # extended Levy process (for nicer plotting)
    doPDF <- !is.null(file)
    if(doPDF) pdf(file=file, width=7, height=7)
    plot(x_jump_time, x_L[,1], type="s", ylim=range(L),
         xlab="t", ylab=expression(bold(L)[t]), col=cols[1], ...)
    for(j in 2:d)
        lines(x_jump_time, x_L[,j], type="s", col=cols[j])
    legend("bottomright", bty="n", lty=rep(1, d), col=cols,
           legend=as.expression( lapply(1:d, function(j) bquote(L[list(t,.(j))]))))
    if(doPDF) dev.off()
}

## -----------------------------------------------------------------------------
## Marginal (variance gamma) parameters
th <- -0.2
kap <- 0.05
sig <- 0.3

## Truncation specification
Gamma.star <- 2000
burn.in <- 500

## ----fig.align="center", fig.width=7.5, fig.height=6--------------------------
## Gamma
theta <- 4 # theta
d <- 4 # dimension
set.seed(271)
system.time(Gamma <- Gamma_Clayton_Levy(d, theta=theta, Gamma.star=Gamma.star,
                                        burn.in=burn.in))
plot_Gamma(Gamma, Gamma.star=Gamma.star,
           file=if(doPDF) "fig_Gamma_positive_Clayton_Levy_copula.pdf" else NULL,
           main=expression(bold(group("(",Gamma[k],")")~
                           "for a positive Clayton Levy copula")))

## (L_t)
L <- hom_vargamma_Levy(Gamma, th=th, kap=kap, sig=sig)
plot_Levy(L, file=if(doPDF) "fig_L_with_positive_Clayton_Levy_copula.pdf" else NULL,
          main="Levy process with positive Clayton Levy copula")

## ----fig.align="center", fig.width=7.5, fig.height=6--------------------------
## Gamma
theta <- c(0.7, 4, 2) # theta_0, theta_1, theta_2
set.seed(271)
system.time(Gamma <- Gamma_nested_Clayton_Levy(theta, Gamma.star=Gamma.star,
                                               burn.in=burn.in))
## 15 seconds on 2015-fast platform
plot_Gamma(Gamma, Gamma.star=Gamma.star,
           file=if(doPDF) "fig_Gamma_positive_nested_Clayton_Levy_copula.pdf" else NULL,
           main=expression(bold(group("(",Gamma[k],")")~
                           "for a positive nested Clayton Levy copula")))

## (L_t)
L <- hom_vargamma_Levy(Gamma, th=th, kap=kap, sig=sig)
plot_Levy(L, file=if(doPDF) "fig_L_with_positive_nested_Clayton_Levy_copula.pdf" else NULL,
          main="Levy process with positive nested Clayton Levy copula")

