#### Testing fitCopula
######################################

require(copula)

## one replicate
do1 <- function(cop,n) {
  x <- rCopula(n, cop)
  u <- pobs(x)
  f.itau <- fitCopula(cop, u, method="itau")
  f.irho <- fitCopula(cop, u, method="irho")
  f.mpl <- fitCopula(cop, u, method="mpl")
  f.ml <- fitCopula(cop, x, method="ml")
  cbind(est = c(itau = f.itau@estimate, irho = f.irho@estimate,
          mpl = f.mpl@estimate, ml   = f.ml@estimate),
        se = c(itau = sqrt(f.itau@var.est), irho = sqrt(f.irho@var.est),
          mpl = sqrt(f.mpl@var.est), ml = sqrt(f.ml@var.est)))
}

## N replicates with mean of estimates and se, and sd of estimates
testCop <- function(cop, n, N, verbose=TRUE) {
     if(verbose) { k <- 0; pb <- txtProgressBar(max=N, style = 3);
                   on.exit(close(pb)) } # setup progress bar and close at end
     sim <- replicate(N, {
          if(verbose) setTxtProgressBar(pb, (k <<- k + 1)) # update progress bar
          do1(cop,n)
     })
     m <- rowMeans(sim, dims=2)
     cbind(bias = abs(m[,"est"] - cop@parameters), ## abs of bias
           dSE  = abs(m[,"se"] - apply(sim[,1,],1,sd))) ## abs of mean of se minus sd of estimates
}

##' Test of the methods in fitCopula for a bivariate one-parameter copula family
##'
##' @title Test of the methods in fitCopula for a bivariate one-parameter copula family
##' @param cop is the copula family
##' @param tau.set is the set of tau values for which the parameters of cop are set
##' @param n.set is the set of n values used in the simulation
##' @param N is the number of repetitions for computing the bias and dSE
##' @param verbose logical indicating if progress bar is desired
##' @return an array bias and dSE in different tau and n scenarios
##' @author Martin
run1test <- function(cop, tau.set, n.set, N, verbose=TRUE) {
  theta <- structure(iTau(cop, tau.set),
		     names = paste("tau", format(tau.set), sep="="))
  names(n.set) <- paste("n",format(n.set),sep="=")
  setPar <- function(cop, value) { cop@parameters <- value ; cop }
  cop.set <- lapply(theta, setPar, cop=cop)
  sapply(cop.set,
         function(cop) {
              if(verbose) cat("copula par : ", format(cop@parameters),"\n")
              sapply(n.set, function(n)
                {
                     if(verbose) cat("  n = ",n,"\n", sep="")
                     testCop(cop, n, N=N, verbose=verbose)
                     ##-----
                }, simplify="array")
         }, simplify = "array")
}


##' Reshapes the results for processing by xyplot
##' @title Reshapes the results for processing by xyplot
##' @param res object returned by run1test
##' @param taunum logical indicating whether the taus are numeric or character
##' @return a matrix containing the reshaped results
##' @author Martin
reshape.res <- function(res, taunum=FALSE) {
  names(dimnames(res)) <- c("method","stat","n","tau")
  d <- cbind(expand.grid(dimnames(res),KEEP.OUT.ATTRS=FALSE), x=as.vector(res))
  d[,"n"  ] <- as.numeric(vapply(strsplit(as.character(d[, "n" ]),split="="),`[`,"",2))
  if (taunum)
    d[,"tau"] <- as.numeric(vapply(strsplit(as.character(d[,"tau"]),split="="),`[`,"",2))
  d
}

## plot the results
plots.res <- function(d, log=TRUE, auto.key=TRUE, ...) {
  require(lattice)
  xyplot(x ~ n | stat*tau, groups=method, data=d, type="b",
         scales = list(x=list(log=log), y=list(log=log)),
         auto.key=TRUE, ...)
}

## test code
system.time(
rr <- run1test(normalCopula(), tau.set=seq(0.2,0.8,by=0.2), n.set=c(25,50,100,200), N=200)
)
## ~ 400 seconds

d <- reshape.res(rr)

plots.res(d)# MM: desirable but ugly tick labelling

plots.res(d, log=FALSE)
