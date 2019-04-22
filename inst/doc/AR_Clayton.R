## ---- message=FALSE------------------------------------------------------
require(copula)
set.seed(271)

## ------------------------------------------------------------------------
mu <- 0
sigma <- 1
df <- 3
alpha <- 10

## ------------------------------------------------------------------------
rtls <- function(n, df, mu, sigma) sigma * rt(n,df) + mu
ptls <- function(x, df, mu, sigma) pt((x - mu)/sigma,df)
qtls <- function(u, df, mu, sigma) sigma * qt(u,df) + mu
dtls <- function(u, df, mu, sigma) dt((x - mu)/sigma,df)/sigma

## ------------------------------------------------------------------------
rclayton <- function(n, alpha) {
	u <- runif(n+1) # innovations
	v <- u
	for(i in 2:(n+1))
            v[i] <- ((u[i]^(-alpha/(1+alpha)) -1)*v[i-1]^(-alpha) +1)^(-1/alpha)
	v[2:(n+1)]
}
n <- 200
u <- rclayton(n, alpha = alpha)
u <- qtls(u, df=df, mu=mu, sigma=sigma)
y <- u[-n]
x <- u[-1]

## ------------------------------------------------------------------------
fitCopula(claytonCopula(dim=2),
          cbind(ptls(x,df,mu,sigma), ptls(y,df,mu,sigma)))

## ------------------------------------------------------------------------
## Identical margins
M2tlsI <- mvdc(claytonCopula(dim=2), c("tls","tls"),
               rep(list(list(df=NA, mu=NA, sigma=NA)), 2), marginsIdentical= TRUE)
(fit.id.mar <- fitMvdc(cbind(x,y), M2tlsI, start=c(3,1,1, 10)))

## Not necessarily identical margins
M2tls <- mvdc(claytonCopula(dim=2), c("tls","tls"),
              rep(list(list(df=NA, mu=NA, sigma=NA)), 2))
fitMvdc(cbind(x,y), M2tls, start=c(3,1,1, 3,1,1, 10))

## ---- fig.align="center", fig.width=6, fig.height=6----------------------
u.cond <- function(z, tau, df, mu, sigma, alpha)
    ((tau^(-alpha/(1+alpha)) -1) * ptls(z,df,mu,sigma)^(-alpha) + 1) ^ (-1/alpha)
y.cond <- function(z, tau, df, mu, sigma, alpha) {
    u <- u.cond(z, tau, df, mu, sigma, alpha)
    qtls(u, df=df, mu=mu, sigma=sigma)
}
plot(x, y)
title("True and estimated conditional quantile functions")
mtext(quote("for" ~~  tau == (1:5)/6))
z <- seq(min(y),max(y),len = 60)
for(i in 1:5) {
    tau <- i/6
    lines(z, y.cond(z, tau, df,mu,sigma, alpha))
    ## and compare with estimate:
    b <- fit.id.mar@estimate
    lines(z, y.cond(z, tau, df=b[1], mu=b[2], sigma=b[3], alpha=b[4]),
          col="red")
}

