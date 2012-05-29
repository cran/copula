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
sessionInfo() # will change often.. but if we need the info, get it here

### A faster, more checking version of demo(estimation.gof)
### that is, of ../demo/estimation.gof.R
##              ~~~~~~~~~~~~~~~~~~~~~~~~
## Note: This is only for proof of concept, the numbers chosen are not reasonable
##       for proper (estimation and) goodness-of-fit testing

source(system.file("Rsource", "estim-gof-fn.R", package="copula"))
## --> estimation.gof() etc

## From source(system.file("test-tools.R", package = "Matrix")) :
showProc.time <- local({
    pct <- proc.time()
    function() { ## CPU elapsed __since last called__
	ot <- pct ; pct <<- proc.time()
	cat('Time elapsed: ', (pct - ot)[1:3],'\n')
    }
})

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

showProc.time()

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

showProc.time()

## Now print RR
options(digits=5)

## No times here...
RR[,c("theta_hat", "tau_hat", "P_value", "< 0.05"),,]

## ... but rather here, separately:
apply(RR[,c("timeEstim","timeGoF"),,], c(3,1,2), mean)

showProc.time()

### Make sure the log-Likelihood demos run: ####################################

demo("logL-vis", package="copula")

showProc.time()


### *Some* minimal  gofCopula() examples {the help page has \dontrun{} !}
set.seed(101)

## A two-dimensional data example -------
x <- rcopula(claytonCopula(3), 200)

for(fitMeth in c("mpl", "ml", "itau", "irho")) {
    cat("\nfit*( method = '", fitMeth,"')\n----------------------\n\n", sep="")
    print(gofCopula(gumbelCopula (1), x, N = 10, print. =0, method = fitMeth))
    print(gofCopula(claytonCopula(1), x, N = 10, print. =0, method = fitMeth))
}
showProc.time()

## The same using the multiplier approach -- "ml" is not allowed:
for(fitMeth in c("mpl", "itau", "irho")) {
    cat("\nfit*( method = '", fitMeth,"')\n----------------------\n\n", sep="")
    print(gofCopula(gumbelCopula (1), x, N = 10, print. =0, method = fitMeth, simulation="mult"))
    print(gofCopula(claytonCopula(1), x, N = 10, print. =0, method = fitMeth, simulation="mult"))
}

## A three-dimensional example  ------------------------------------
x <- rcopula(tCopula(c(0.5, 0.6, 0.7), dim = 3, dispstr = "un"),
             200)

gumbC <- gumbelCopula(1, dim = 3)
t.cop <- tCopula(rep(0, 3), dim = 3, dispstr = "un", df.fixed=TRUE)

showProc.time()

for(fitMeth in c("mpl", "ml", "itau", "irho")) {
    cat("\nfit*( method = '", fitMeth,"')\n----------------------\n\n", sep="")
    print(gofCopula(gumbC, x, N = 10, print. =0, method = fitMeth))
    print(gofCopula(t.cop, x, N = 10, print. =0, method = fitMeth))
}
showProc.time()

## The same using the multiplier approach -- "ml" is not allowed in general;
##  "itau" and "irho"  only  for  d = 2  (for now !)
for(fitMeth in c("mpl")) {
    cat("\nfit*( method = '", fitMeth,"')\n----------------------\n\n", sep="")
    print(gofCopula(gumbC, x, N = 10, method = fitMeth, simulation="mult"))
    print(gofCopula(t.cop, x, N = 10, method = fitMeth, simulation="mult"))
}
showProc.time()
