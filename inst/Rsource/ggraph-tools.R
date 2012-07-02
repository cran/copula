## Copyright (C) 2012 Marius Hofert and Martin Maechler
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


### Conditional copula functions ###############################################

##' Conditional copula function C(u2|u1) of u2 given u1
##'
##' @title Bivariate ("2") Conditional Copula function
##' @param u2 data vector (numeric(n))
##' @param u1 data vector (numeric(n))
##' @param family copula family
##' @param theta parameter (for ACs; for elliptical copulas, its rho)
##' @param ... additional args (e.g., df for t copulas)
##' @return C(u2|u1)
##' @author Marius Hofert and Martin Maechler
ccop2 <- function(u2, u1, family, theta, ...) {
    stopifnot(length(u1)==length(u2))
    switch(copFamilyClass(family),
           "elliptical"={
               switch(family,
                      "normal"={
                          pnorm((qnorm(u2)-theta*qnorm(u1))/sqrt(1-theta^2))
                      },
                      "t"={
                          stopifnot(hasArg(df))
                          df <- list(...)$df
                          qt.u1 <- qt(u1, df=df)
                          mu <- theta * qt.u1
                          sigma <- sqrt((df+qt.u1^2)*(1-theta^2)/(df+1))
                          pt((qt(u2, df=df)-mu)/sigma, df=df+1)
                      },
                      stop("not yet supported family"))
           },
           "nArchimedean"={
               ## cacopula(u, cop=onacopulaL(family, list(theta, 1:2)))
               cop <- getAcop(family)
               u <- cbind(u1, u2)
               psiI <- cop@psiInv(u, theta=theta)
	       exp(cop@psiDabs(rowSums(psiI), theta=theta, log=TRUE) -
		   cop@psiDabs(psiI[,1], theta=theta, log=TRUE))
           },
           stop("family ", family, " not yet supported"))
}

##' Compute pairwise Rosenblatt-transformed variables C(u[,i]|u[,j])) for the given
##' matrix u
##'
##' @title Compute pairwise Rosenblatt-transformed variables
##' @param u (n, d)-matrix (typically pseudo-observations, but also perfectly
##'        simulated data)
##' @param cop copula object used for the Rosenblatt transform (H_0;
##'        either nArchimedean or elliptical)
##' @param ... additional arguments passed to ccop2
##' @return (n,d,d)-array cu with cu[,i,j] containing C(u[,i]|u[,j]) for i!=j
##'         and u[,i] for i=j
##' @author Marius Hofert
pairwiseCcop <- function(u, cop, ...)
{
    if(!is.matrix(u)) u <- rbind(u)
    stopifnot((d <- ncol(u)) >= 2, 0 <= u, u <= 1,
	      d == dim(cop))

    ## (1) determine copula class and compute auxiliary results
    cls <- copClass(cop)
    switch(cls,
           "elliptical"={
               ## build correlation matrix from vector of copula parameters
               P <- matrix(0, nrow=d, ncol=d)
	       rho <- cop@parameters
	       P[lower.tri(P)] <- if(cop@df.fixed) rho else rho[-length(rho)]
               P <- P + t(P)
               diag(P) <- rep.int(1, d)
           },
           "nArchimedean"={
               ## build "matrix" of dependence parameters
               P <- nacPairthetas(cop)
           },
           stop("not yet supported copula object"))

    ## (2) compute pairwise C(u_i|u_j)
    family <- copFamily(cop) # determine copula family
    n <- nrow(u)
    cu.u <- array(NA_real_, dim=c(n,d,d), dimnames=list(C.ui.uj=1:n, ui=1:d, uj=1:d))
    ## cu.u[,i,j] contains C(u[,i]|u[,j]) for i!=j and u[,i] for i=j
    for(i in 1:d) { # first index C(u[,i]|..)
        for(j in 1:d) { # conditioning index C(..|u[,j])
            cu.u[,i,j] <- if(i==j) u[,i] else
            ccop2(u[,i], u[,j], family, theta=P[i,j], ...)
        }
    }
    cu.u
}


##' Compute matrix of pairwise tests for independence
##'
##' @title Compute matrix of pairwise tests for independence
##' @param cu.u (n,d,d)-array as returned by \code{pairwiseCcop()}
##' @param N.sim number of simulation \code{N} for \code{\link{indepTestSim}()}
##' @param verbose logical indicating if and how much progress info should be printed.
##' @param ... additional arguments passed to indepTestSim()
##' @return (d,d)-matrix of lists with test results (as returned by indepTest())
##'         for the pairwise tests of independence
##' @author Marius Hofert
pairwiseIndepTest <- function(cu.u, N.sim=256, verbose=TRUE, ...)
{
    require(copula)

    ## (1) simulate test statistic under independence
    stopifnot(length(dim. <- dim(cu.u)) == 3)
    stopifnot(dim.[2]==dim.[3])
    n <- dim.[1]
    d <- dim.[2]
    if(verbose) cat(sprintf("pairwiseIndepTest( (n=%d, d=%d)): indepTestSim(*, N=%d) .. ",
                            n,d, N.sim))
    indepTestDistObj <- indepTestSim(n, p=2, m=2, N=N.sim, ...)
    if(verbose) cat(" *Sim  done\n")

    ## (2) compute matrix of pairwise p-values for test of independence
    p <- matrix(list(), nrow=d, ncol=d)
    for(j in 1:d) { # colum
        if(verbose) if(verbose >= 2) cat("j = ", j, ";  i =", sep="") else cat(j, "")
        uj <- cu.u[,,j] # (n x d)
        for(i in 1:d) { # row
            if(verbose >= 2) cat(i,"")
            p[i,j][[1]] <- if(j==i) list(fisher.pvalue = NA)
            else indepTest(cbind(uj[,j], uj[,i]), d=indepTestDistObj)
        }
        if(verbose >= 2) cat("\n")
    }
    p
}

##' Extract p-values from a matrix of indepTest objects
##'
##' @title Extract p-values from a matrix of indepTest objects
##' @param piTest matrix of indepTest objects
##' @return matrix of p-values
##' @author Marius Hofert, Martin Maechler
##' Note: Kojadinovic, Yan (2010) mention that "fisher.pvalue" performs best;
##'       In d=2 (as we use in pairwiseIndepTest) all three methods are equal.
pviTest <- function(piTest){
    matrix(unlist(lapply(piTest, function(x) x$fisher.pvalue)),
           ncol=ncol(piTest))
}

##' Computing a global p-value
##'
##' @title Computing a global p-value
##' @param pvalues (matrix of pairwise) p-values
##' @param method vector of methods for pvalues (how the p-values are adjusted)
##' @param globalFun function determining how to compute a global p-value from a
##'        matrix of pairwise adjusted p-values
##' @return global p-values for each of the specified methods (of how to adjust the
##'         pairwise p-values)
##' @author Marius Hofert
gpviTest <- function(pvalues, method=p.adjust.methods, globalFun=min){
    pvalues <- pvalues[!is.na(pvalues)] # vector
    if(all(c("fdr","BH") %in% method))## p.adjust():  "fdr" is alias for "BH"
	method <- method[method != "fdr"]
    sapply(method, function(meth) globalFun(p.adjust(pvalues, method=meth)))
}

gpviTest0 <- function(pvalues) {
    pvalues <- pvalues[!is.na(pvalues)] # vector
    c("minimum" = min(pvalues),
      "global (Bonferroni/Holm)" = min(p.adjust(pvalues, method="holm")))
}

gpviString <- function(pvalues, sep = "   ", digits = 2) {
    pv <- gpviTest0(pvalues)
    paste0("p-values:", sep,
	  paste(paste(names(pv), sapply(pv, format.pval, digits=digits),
		      sep=": "),
		collapse=paste0(";",sep)))
}
