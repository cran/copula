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


##------- FIXME{MM}: Shouldn't all these inherit from "htest" (and use its print() method)?
##        and ditto for all other ./*Tests.R  !!!


### EV test based on Cn ########################################################

evTestC <- function(x, N = 1000)
{
  ## make pseudo-observations
  p <- ncol(x)
  n <- nrow(x)
  u <- apply(x,2,rank)/(n+1)

  ## set r according to recommendations
  r <- 3:5
  nr <- length(r)

  ## grid = pseudo-observations
  m <- 0

  ## parameters
  offsetstat <- 0.75
  der2n <- TRUE

  ## make grid
  if (m > 0)
    {
      y <- seq(1/m, 1 - 1/m, len = m)
      v <- vector("list", p)
      for (i in 1:p)
        v[[i]] <- y
      g <- as.matrix(expand.grid(v))
      m <- nrow(g)
    }
  else
    {
      g <- u
      m <- n
    }

  out <- .C(evtest,
            as.double(u),
            as.integer(n),
            as.integer(p),
            as.double(g),
            as.integer(m),
            as.integer(N),
            as.double(1/r),
            as.integer(nr),
            s0 = double(N * nr),
            as.integer(der2n),
            as.double(offsetstat),
            stat = double(nr))[c("s0", "stat")]

  ## s <- out$stat
  ## s0 <- matrix(out$s0, ncol = nr, byrow = TRUE)
  ## pval <- apply(s0 >= matrix(s, nrow=N, ncol=nr, byrow=TRUE),
  ##               2, function(x) (sum(x) + 0.5) / (N + 1) )
  ## comb.s <- sum(s)
  ## comb.pval <- ( sum( apply(s0, 1, sum) >= comb.s ) + 0.5 ) / (N + 1)

  ## return(list(statistic=c(s, comb.s),
  ##             pvalue=c(pval, comb.pval), s0=s0))

  ## p-values
  s0 <- matrix(out$s0, ncol = nr, byrow = TRUE)
  comb.s <- sum(out$stat)
  comb.pval <- ( sum( apply(s0, 1, sum) >= comb.s ) + 0.5 ) / (N + 1)
  evt <- list(statistic=comb.s,pvalue=comb.pval)
  class(evt) <- "evTest"
  evt
}

print.evTest <- function(x, ...)
{
  cat("Statistic:", x$statistic,
      "with p-value", x$pvalue, "\n\n")
  invisible(x)
}


### EV test based on An ########################################################

evTestA <- function(x, N = 1000, derivatives = "An")
{
  ## make pseudo-observations
  n <- nrow(x)
  u <- apply(x,2,rank)/(n+1)

  ## make grid
  ## m = 0
  g <- u
  m <- n

  estimator <- "CFG"
  offset <- 0.5

  ## compute the test statistic
  s <- .C(evtestA_stat,
          as.double(u[,1]),
          as.double(u[,2]),
          as.integer(n),
          as.double(g[,1]),
          as.double(g[,2]),
          as.integer(m),
          as.integer(estimator == "CFG"),
          stat = double(1),
          as.double(offset))$stat

  if (derivatives == "Cn")
    s0 <- .C(evtestA,
             as.double(u[,1]),
             as.double(u[,2]),
             as.integer(n),
             as.double(g[,1]),
             as.double(g[,2]),
             as.integer(m),
             as.integer(estimator == "CFG"),
             as.integer(N),
             s0 = double(N))$s0
  else
    s0 <- .C(evtestA_derA,
             as.double(u[,1]),
             as.double(u[,2]),
             as.integer(n),
             as.double(g[,1]),
             as.double(g[,2]),
             as.integer(m),
             as.integer(estimator == "CFG"),
             as.integer(N),
             s0 = double(N))$s0

  evt <- list(statistic=s, pvalue=(sum(s0 >= s)+0.5)/(N+1))
  class(evt) <- "evTest"
  evt
}


### EV test based on An CFG - An Pickands ######################################

evTestAA <- function(x, N = 1000,  derivatives = "Cn", m = 100)
{
  ## make pseudo-observations
  n <- nrow(x)
  u <- apply(x,2,rank)/(n+1)

  ## make grid
  g <- seq(1/m, 1 - 1/m, len = m)

  ## compute the test statistic
  s <- .C(evTestAA_stat,
          as.double(-log(u[,1])),
          as.double(-log(u[,2])),
          as.integer(n),
          as.double(g),
          as.integer(m),
          stat = double(1))$stat

  if (derivatives == "Cn")
    s0 <- .C(evTestAA,
             as.double(u[,1]),
             as.double(u[,2]),
             as.integer(n),
             as.double(g),
             as.integer(m),
             as.integer(N),
             s0 = double(N))$s0
  else
    s0 <- .C(evTestAA_derA,
             as.double(u[,1]),
             as.double(u[,2]),
             as.integer(n),
             as.double(g),
             as.integer(m),
             as.integer(N),
             s0 = double(N))$s0

  structure(class = "evTest",
            list(statistic=s, pvalue=(sum(s0 >= s)+0.5)/(N+1),s0=s0))
}


### EV test based on K - Ben Ghorbal, Neslehova and Genest (2009)
### Canadian Journal of Statistics, volume 37
### Code provided by Johanna Neslehova

## internal functions
I <- function(X,i,j)
{
  if (X[i,1] <= X[j,1] && X[i,2] <= X[j,2]) {ind <- 1} else {ind <- 0}
  ind
}

zero.mat <- function(n)
{
  x <- rep(0,n)
  y <- rep(1,n)
  A <- NULL
  for(k in 1:n){
    A <- rbind(A,matrix(c(rep(y,(k-1)),x,rep(y,(n-k))),byrow=TRUE,nrow=n))
  }
  A
}

mult.row <- function(M,delete=TRUE)
{
  n <- ncol(M)
  out <- matrix(rep(as.vector(t(M)),n),ncol=n,byrow=TRUE) *
         matrix(as.vector(t(matrix(rep(as.vector(M),n),nrow=n))),byrow=TRUE,ncol=n)
  if(delete) out*zero.mat(n) else out
}


ind.matrix <- function(X)
{
  n <- nrow(X)
  x <- as.numeric(X[,1])
  y <- as.numeric(X[,2])
  fun <- function(x,y) as.numeric(!(x>y))
  M1 <- outer(x,x, FUN=fun)
  M2 <- outer(y,y, FUN=fun)
  M1*M2
}

## calculating the thetas for GKRstatistic
thetas <- function(X)
{
  mu <- numeric(7)
  theta <- numeric(8)
  n <- nrow(X)
  D <- matrix(1,nrow=n,ncol=n)-diag(n)
  M <- ind.matrix(X)*D #deleting diagonal entries
  N <- t(M)  #transpose of M
  O <- M%*%N
  P <- O*D  #deleting diagonal entries
  A <- matrix(apply(M,1,sum),ncol=1)
  B <- matrix(apply(M,2,sum),ncol=1)
  B.sq <- B*B
  MB <- M%*%B
  ab <- sum(A*B)
  asq <- sum(A*A)
  absq <- sum(A*B.sq)
  BMB <- sum(B*MB)
  AMB <- sum(A*MB)
  BsqMB <- sum(B.sq*MB)
  MN <- sum(O)
  MBMB <- sum(MB*MB)
  b1 <- sum(B)
  b2 <- sum(B.sq)
  b3 <- sum(B.sq*B)
  b4 <- sum(B.sq*B.sq)
  c <- numeric(8)
  mu[1] <- b1
  mu[2] <- b2 - mu[1]
  mu[3] <- b3 - 3*mu[2] - mu[1]
  mu[4] <- b4 - 6*mu[3] - 7*mu[2] - mu[1]
  theta[1] <- ab
  theta[2] <- asq-b1
  theta[3] <- absq-ab
  theta[4] <- BMB-2*theta[1]
  theta[5] <- AMB-b2-theta[2]-theta[1]
  theta[8] <- sum(P*P)-mu[2]
  theta[6] <- BsqMB - 2*absq-theta[3]-theta[4]
  theta[7] <- MBMB-b3-absq-AMB+b2-theta[8]-theta[3]-theta[5]
  mu[5] <- sum(B%*%t(B)) - 2*theta[1] -MN-theta[2]
  mu[6] <- b2*b1 - 2*BMB - 2*AMB  - (b1)^2 + 3*b2 + 4*ab + 2*asq -2*b1 - theta[3] - b3
  ##c[1] <- #BsqMB - 3*absq
  c[1] <- mu[6]#2*(sum((A*B)*MB) - sum((A*A)*B) - 2*sum(t(B)%*%M%*%A) + ab)
  c[2] <- 2*theta[7]#b1*b2 - (b1)^2 - 2*AMB + asq + 4*ab + 3*b2  - b1 - 2*BMB -b3 - absq
  c[3] <- BsqMB - 3*absq + 3*ab - BMB #BsqMB
  c[4] <- 2*(MBMB - b3 - AMB + b2 - absq) #2*(MBMB - b3 - absq)
  c[5] <- BsqMB - BMB #b1*MN - b3 -2*BMB
  c[6] <- BsqMB - absq #BsqMB
  c[7] <- b2*b1 - b3 - absq #b4
  mu[7] <- (b2)^2 - b4 - BsqMB - sum(c)
  return(list(mu=mu,theta=theta))
}

## calculation of the GKR test statistic Sn
Sn <- function(X)
{
  n <- nrow(X)
  D <- matrix(1,nrow=n,ncol=n)-diag(n)
  M <- ind.matrix(X)*D  #delete entries on the diagonal
  w <- sum(M)
  w.sq <- sum((M%*%t(M))*D)
  -1 + (8*w)/(n*(n-1)) - (9*w.sq)/(n*(n-1)*(n-2))
}

## Jackknife variance estimator THIS CODE CAN BE IMPROVED!!!!!!!
GKRJack <- function(X)
{
    stopifnot(is.numeric(n <- nrow(X)))
    Sn <- Sn(X)
    VSn <- vapply(1:n, function(i) Sn(X[-i, ,drop=FALSE]), 1.)
    list(Sn=Sn, var = (n-1)/n * sum(VSn - Sn)^2)
}

## Calculation of Sn and its finite sample and asymptotic variance;
## theta41 = theta[4], theta42 = theta[5], theta51 = theta[6] and theta52=theta[7]
## If variance="all", the function gives the asymptotic (var[1]), finite sample (var[2])
## and jackknife variance (var[3]) of sqrt(n)*Sn.
GKRstatistic <- function(X, variance=c("fsample","asymptotic","all"))
{
  variance <- match.arg(variance)
  n <- nrow(X)
  psi <- numeric(5)
  mu <- numeric(7)
  var <- NULL
  tmp <- thetas(X)
  mu <- tmp$mu
  w <- mu[1]
  w.sq <- mu[2]
  tau <- -1+(4*w)/(n*(n-1))
  Sn <- -1+(8*w)/(n*(n-1))-(9*w.sq)/(n*(n-1)*(n-2))
  theta <- tmp$theta
  mu[1] <- (mu[1])/(n*(n-1))
  mu[2] <- (mu[2])/(n*(n-1)*(n-2))
  mu[3] <- (mu[3])/(n*(n-1)*(n-2)*(n-3))
  mu[4] <- (mu[4])/(n*(n-1)*(n-2)*(n-3)*(n-4))
  mu[5] <- (mu[5])/(n*(n-1)*(n-2)*(n-3))
  mu[6] <- (mu[6])/(n*(n-1)*(n-2)*(n-3)*(n-4))
  mu[7] <- (mu[7])/(n*(n-1)*(n-2)*(n-3)*(n-4)*(n-5))
  psi[1] <- 64*(mu[2] - 4*mu[5]) - 144*(mu[3]-6*mu[6]) + 81*(mu[4] - 9*mu[7])
  psi[2] <- (2*theta[1]+theta[2])/(n*(n-1)*(n-2))
  psi[3] <- (theta[3] + 2*(theta[4]+theta[5]))/(n*(n-1)*(n-2)*(n-3))
  psi[4] <- (4*(theta[6]+theta[7]))/(n*(n-1)*(n-2)*(n-3)*(n-4))
  if((variance=="asymptotic") || (variance=="all")){
    var = c(var,psi[1]+64*psi[2]-144*psi[3]+81*psi[4])
  }
  if((variance == "fsample") || (variance=="all")){
    psi[5] <- (theta[8]+4*theta[3])/(n*(n-1)*(n-2)*(n-3))
    theta[1] <- (theta[1])/(n*(n-1)*(n-2))
    An <- mu[1]+(n-2)*(mu[2]+psi[2])+(n-2)*(n-3)*mu[5]
    Bn <- 2*mu[2]+2*theta[1]+(n-3)*(mu[3]+psi[3])+(n-3)*(n-4)*mu[6]
    Cn <- 2*mu[2] + (n-3)*(4*mu[3]+2*psi[5])+(n-3)*(n-4)*(mu[4]+psi[4])+(n-3)*(n-4)*(n-5)*mu[7]
    var.tmp <- (64*An)/(n*(n-1))-(144*Bn)/(n*(n-1))+(81*Cn)/(n*(n-1)*(n-2))-64*mu[5] + 144*mu[6] - 81*mu[7]
    var <- c(var,n*var.tmp)
  }
  if((variance=="all")){
    VSn <- numeric(n)
    for(i in 1:n){
      cond <- !(1:n == i)
      VSn[i] <- Sn(X[cond,])
    }
    var <- c(var,n*((n-1)/n)*(sum((VSn-rep(Sn,n))^2)))
  }
  return(list(Sn=Sn,tau=tau,mu=mu,psi=psi,var=var))
}


### Test function: GKR test ####################################################

evTestK <- function(x, method = c("fsample","asymptotic","jackknife"))
{
  method <- match.arg(method)
  n <- nrow(x)
  negvar <- FALSE
  tmp <- switch(method,
                "fsample"= GKRstatistic(x,variance="fsample"),
                "asymptotic"=GKRstatistic(x,variance="asymptotic"),
                "jackknife"=GKRJack(x)
                )
  var <- switch(method,
                "fsample"=tmp$var,
                "asymptotic"=tmp$var,
                "jackknife"=n*tmp$var
                )
  if(var < 0){
    message("Variance estimator less then zero, using jackknife instead")## <- MM: reactivated 2012-02-22
    negvar <- TRUE
    var <- n*GKRJack(x)$var
    method <- "jackknife"
  }

  Tn <- sqrt(n)*tmp$Sn/sqrt(var)
  ##calpha <- qnorm((1-alpha/2))
  ##reject <- (abs(Tn)>calpha)
  p.value <- pnorm(-abs(Tn))+pnorm(abs(Tn),lower.tail=FALSE)
  evt <- list(statistic=Tn, pvalue=p.value, negvar=negvar, method=method)
  class(evt) <- "evTest"
  evt
}
