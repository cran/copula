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


require(copula)
sessionInfo() # will change too often.. but if we need the info

### A faster, more checking version of demo(estimation.gof)
### that is, of ../demo/estimation.gof.R
##              ~~~~~~~~~~~~~~~~~~~~~~~~
## Note: This is only for proof of concept, the numbers chosen are not reasonable
##       for proper (estimation and) goodness-of-fit testing

source(system.file("Rsource", "estim-gof-fn.R", package="copula"))
## --> estimation.gof() etc

## Use GoF methods:
(gofTraf <- eval(formals(gnacopula)$trafo))
(gofMeth <- eval(formals(gnacopula)$method))

set.seed(1) # set seed

n <- 64 # sample size [small here for CPU reasons]
d <- 5 # dimension
tau <- 0.25 # Kendall's tau

### apply all procedures (to data from AMH) ####################################

simFamily <- "AMH"
cop <- getAcop(simFamily)
theta <- cop@tauInv(tau) # true parameter

## start the loop
cat("\n### data from ",simFamily," (n = ",n,", d = ",d,", theta = ",
    format(theta),", tau = ", format(tau),") ###\n\n",sep="")

## note: this might (still) take a while...
RR <- sapply(gofTraf, simplify="array", function(gt)
         {
             sapply(gofMeth, simplify="array", function(gm)
                    estimation.gof(n, d=d, simFamily=simFamily, tau=tau,
                                   n.bootstrap=1, # << "nonsense" for speed reasons;..
### for a particular method under consideration, please choose a larger number here, for example 1000
                                   include.K=TRUE, esti.method = "mle",
                                   gof.trafo=gt, gof.method=gm))
         })
str(RR)
dimnames(RR)

## Now print RR
options(digits=5)

## No times here...
RR[,c("theta_hat", "tau_hat", "P_value", "< 0.05"),,]

## ... but rather here, separately:
apply(RR[,c("timeEstim","timeGoF"),,], c(3,1,2), mean)


cat('Time elapsed: ', proc.time(),'\n') # for ''statistical reasons''

### Make sure the log-Likelihood demos run: ####################################

demo("logL-vis", package="copula")

cat('Time elapsed: ', proc.time(),'\n') # for ''statistical reasons''
