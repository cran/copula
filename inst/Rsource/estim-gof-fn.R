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


### is source()d from  demo(estimation.gof)  and also in some of the tests
##		       ------------------- ==> ../../demo/estimation.gof.R
##					       ~~~~~~~~~~~~~~~~~~~~~~~~~~~

##' measures user run time in milliseconds
utms <- function(x) 1000 * system.time(x)[[1]]
##' formats such utms() times "uniformly":
f.tms <- function(x) paste(round(x),"ms") # with a space (sep=" ") !

##' Demonstrates the fitting and goodness-of-fit capabilities for Archimedean
##' copulas
##'
##' @title Fitting and Goodness-Of-Fit for Archimedean copulas
##' @param n sample size
##' @param d dimension
##' @param simFamily Archimedean family to be sampled
##' @param tau degree of dependence of the sampled family in terms of Kendall's tau
##' @param n.bootstrap
##' @param include.K
##' @param n.MC if > 0 it denotes the sample size for SMLE
##' @param esti.method estimation method (see enacopula)
##' @param gof.method goodness-of-fit transformation (see gnacopula)
##' @param checkFamilies vector of Archimedean families to be used for gof
##' @param verbose
##' @param n.bootstrap, see gnacopula()
##' @return a numeric matrix ...
##' @author Marius Hofert and Martin Maechler
estimation.gof <- function(n, d, simFamily, tau,
			   n.bootstrap = 1, # dummy number of bootstrap replications
			   include.K = TRUE,
			   n.MC = if(esti.method=="smle") 10000 else 0,
			   esti.method = eval(formals(enacopula)$method),
			   gof.trafo =	 eval(formals(gnacopula)$trafo),
			   gof.method =	 eval(formals(gnacopula)$method),
			   checkFamilies = copula:::c_longNames,
			   verbose = TRUE)
{
    ## generate data
    copFamily <- getAcop(simFamily)
    theta <- copFamily@tauInv(tau)
    cop <- onacopulaL(simFamily, list(theta,1:d))
    if(verbose){
        r <- function(x) round(x,4) # for output
	cat("\n\n### Output for esti.method = \"",esti.method,
            "\" gof.trafo = \"",gof.trafo,"\", and gof.method = \"",gof.method,"\"\n\n",sep="")
    }
    u <- rnacopula(n,cop)

    ## estimation and gof
    esti.method <- match.arg(esti.method)
    gof.trafo <- match.arg(gof.trafo)
    gof.method <- match.arg(gof.method)
    n.cf <- length(checkFamilies)
    sig <- logical(n.cf)
    tau <- gof <- ute <- utg <- est <- numeric(n.cf)
    for(k in seq_len(n.cf)) {
        ## estimate the parameter with the provided method
        cop.hat <- onacopulaL(checkFamilies[k],list(NA,1:d))
        if(verbose) cat("Estimation and GoF for ",checkFamilies[k],":\n\n",sep="")
	ute[k] <- utms(est[k] <- enacopula(u, cop=cop.hat,
					   method=esti.method, n.MC=n.MC))
        tau[k] <- cop.hat@copula@tau(est[k])
        if(verbose){
            cat("   theta hat      = ",r(est[k]),"\n",
                "   tau hat        = ",r(cop.hat@copula@tau(est[k])),"\n",
                ## The exact string 'Time ' must appear at the beginning of line
                ## for 'R CMD Rdiff' to *not* look at differences there:
                "Time estimation   = ",f.tms(ute[k]),"\n", sep="")
	}
	cop.hat@copula@theta <- est[k]
        ## apply a goodness-of-fit test to the estimated copula
        ## {{ use rtrafo() or htrafo() if you want the transformed u }}
	utg[k] <-
	    utms(gof[k] <-
		 gnacopula(u, cop=cop.hat, n.bootstrap=n.bootstrap,
			   estimation.method=esti.method,
			   include.K=include.K, n.MC=n.MC,
			   trafo=gof.trafo, method=gof.method,
			   verbose=FALSE)$p.value)
	sig[k] <- (gof[k] < 0.05) # TRUE/FALSE <--> 1/0 -- Careful: may be  NA !
	if(verbose)
	    cat("   p-value	   = ",r(gof[k]),"\n",
		"   < 0.05	   = ", format(sig[k]), "\n",
		"Time GoF comp = ",f.tms(utg[k]),"\n\n", sep="")
    }

    ## result
    names(est) <- checkFamilies
    cbind(theta_hat=est, tau_hat  =tau,
	  timeEstim=ute,
	  P_value  =gof, "< 0.05" =sig,
	  timeGoF  =utg)
}
