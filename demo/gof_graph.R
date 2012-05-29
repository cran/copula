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


### setup ######################################################################

require(copula)

copFile <- function(f) system.file("Rsource", f, package="copula", mustWork=TRUE)
source(copFile("wrapper.R"))
source(copFile("ggraph-tools.R"))
source(copFile("ggraph-graphics.R"))

setSeeds <- TRUE

### Example 1: 5d Gumbel copula ################################################

## setup
n <- 1000 # sample size
d <- 5 # dimension
family <- "Gumbel" # copula family
tau <- 0.5
if(setSeeds) set.seed(1)

## define and sample the copula (= H0 copula), build pseudo-observations
cop <- getAcop(family)
th <- cop@tauInv(tau) # correct parameter value
copH0 <- onacopulaL(family, list(th, 1:d)) # define H0 copula
U. <- pobs(rcop(n, cop=copH0))

## create array of pairwise copH0-transformed data columns
cu.u <- pairwiseCcop(U., copH0)
stopifnot(is.array(cu.u), dim(cu.u) == c(n,d,d)) # check

## compute pairwise matrix of p-values and corresponding colors
pwIT <- pairwiseIndepTest(cu.u, verbose=TRUE) # (d,d)-matrix of test results
round(pmat <- pviTest(pwIT), 3) # pick out p-values
cc <- pairsColList(pmat) # compute corresponding colors

## which pairs violate H0?
which(pmat < 0.05, arr.ind=TRUE)

## check whether p-values are uniform
if(d > 10){
    n. <- d*(d-1)
    qqplot(qunif(ppoints(n.)), sort(pmat), main = paste("n = ",n.))
    abline(0,1)
}


### Example 1: Plots ###########################################################

## 1) plain
pairsRosenblatt(cu.u, pvalueMat=pmat)

## 2) with title, no subtitle
pwRoto <- "Pairwise Rosenblatt transformed observations"
pairsRosenblatt(cu.u, pvalueMat=pmat, pch=".", main=pwRoto, sub=NULL)

## 3) with title and manual subtitle
gp <- format(gpviTest(pmat), digits=1, nsmall=1)
sub <- paste(names(gp), gp, sep=": ")
sub. <- paste(paste(sub[1:3], collapse=", "), "\n",
              paste(sub[4:7], collapse=", "), sep="")
pairsRosenblatt(cu.u, pvalueMat=pmat, pch=".", main=pwRoto, sub=sub., sub.line=5.4)

## 4) two-line title including expressions, and centered  --- JCGS, Fig.1(left) ---
title <- list(paste(pwRoto, "to test"),
              substitute(italic(H[0]:C~~bold("is Gumbel with"~~tau==tau.)),
                         list(tau.=tau)))
main.line <- c(4, 1.4)
pairsRosenblatt(cu.u, pvalueMat=pmat, pch=".",
                main=title, main.line=main.line, main.centered=TRUE)

## 5) omit panel borders
pairsRosenblatt(cu.u, pvalueMat=pmat, pch=".", panel.border=FALSE)

## 6) without axes
pairsRosenblatt(cu.u, pvalueMat=pmat, pch=".", axes=FALSE)

## 7) without axes and borders
pairsRosenblatt(cu.u, pvalueMat=pmat, pch=".", axes=FALSE, panel.border=FALSE)

## 8) make black colors of the dots less dominant
pairsRosenblatt(cu.u, pvalueMat=pmat, pch=".", col=adjustcolor("black", 0.5))

## 9) plot just colors (axis labels are automagically removed)
pairsRosenblatt(cu.u, pvalueMat=pmat, method="none")

## 10) also remove labels on the diagonal
pairsRosenblatt(cu.u, pvalueMat=pmat, method="none", labels="n")

## 11) Q-Q plots -- can, in general, better detect outliers
pairsRosenblatt(cu.u, pvalueMat=pmat, method="QQchisq", cex=0.2)
## pairsRosenblatt(cu.u, pvalueMat=pmat, method="QQchisq", cex=0.2,
##                 panel=function(x, y, ...){
##                     points(x, y, ...)
##                     qqline(y) ## only makes sense for the normal distribution
##                 })
## => TODO: maybe in the future with a more general function for qqline (supporting
##          other distributions)

## 12) P-P plots -- actually, MM sees *more* (though outliers are invisible)
pairsRosenblatt(cu.u, pvalueMat=pmat, method="PPchisq")
pairsRosenblatt(cu.u, pvalueMat=pmat, method="PPchisq",
                panel=function(x, y, ...){
                    points(x, y, ...)
                    abline(0, 1, col="green") # add straight line
                })


### Example 1: Boundary cases ##################################################

## Note: this is only for checking "boundary input cases", it does not make
##       sense given the data.

## 1) one pdiv, within range(pmat)
pmat <- matrix(c(
               NA, 0.1, 0.2, 0.3, 0.4,
               0.2, NA, 0.3, 0.4, 0.5,
               0.3, 0.4, NA, 0.5, 0.6,
               0.4, 0.5, 0.6, NA, 0.7,
               0.5, 0.6, 0.7, 0.8, NA
               ), nrow=5, ncol=5)
bgColList <- pairsColList(pmat, pdiv=0.2, signif.P=0.2,
                          colsBelow="green", colsAbove="orange")
pairsRosenblatt(cu.u, pvalueMat=pmat, pch=".", bgColList=bgColList)

## 2) one pdiv, pdiv < min(pmat)
bgColList <- pairsColList(pmat, pdiv=0.05, signif.P=0.05,
                          colsBelow="green", colsAbove="orange") # note: we don't even need colsBelow here
pairsRosenblatt(cu.u, pvalueMat=pmat, pch=".", bgColList=bgColList)

## 3) one pdiv, pdiv > max(pmat)
bgColList <- pairsColList(pmat, pdiv=0.9, signif.P=0.9, colsBelow="green", colsAbove="orange") # note: we don't even need colsAbove here
pairsRosenblatt(cu.u, pvalueMat=pmat, pch=".", bgColList=bgColList)

## 4) one pdiv, equal to all values of pmat
pmat <- matrix(0.05, nrow=5, ncol=5)
diag(pmat) <- NA
bgColList <- pairsColList(pmat, pdiv=0.05, signif.P=0.05, colsBelow="green", colsAbove="orange") # note: we don't even need colsBelow here
pairsRosenblatt(cu.u, pvalueMat=pmat, pch=".", bgColList=bgColList)
## => plot color key up to 1 (see pairsColList())

## 5) one p-value equal to 0; in this case we need pmin0 > 0
pmat[1,2] <- 0
bgColList <- pairsColList(pmat, signif.P=0.05, pmin0=1e-5)
pairsRosenblatt(cu.u, pvalueMat=pmat, pch=".", bgColList=bgColList)

## 6) specifically call pairsColList but forget to set pmin0 to check the error message
if(FALSE){
    bgColList <- pairsColList(pmat, signif.P=0.05)
    pairsRosenblatt(cu.u, pvalueMat=pmat, pch=".", bgColList=bgColList)
}


### Example 2: (2,3)-nested Gumbel copula ######################################

## setup
n <- 1000 # sample size
d <- 5 # dimension
family <- "Gumbel" # copula family
if(setSeeds) set.seed(2)

## define and sample the copula, build pseudo-observations
cop <- getAcop(family)
th <- cop@tauInv(tau <- c(0.2, 0.4, 0.6))
nacList <- list(th[1], NULL, list(list(th[2], 1:2), list(th[3], 3:d)))
copG <- copCreate(family, nacList=nacList)
U <- rcop(n, cop=copG)
U. <- pobs(U)

## define the H0 copula
th0 <- cop@tauInv(tau0 <- c(0.2, 0.4, 0.4)) # wrong 2nd-sector-parameter
nacList <- list(th0[1], NULL, list(list(th0[2], 1:2), list(th0[3], 3:d)))
copH0 <- onacopulaL(family, nacList)

## create array of pairwise copH0-transformed data columns
cu.u <- pairwiseCcop(U., copH0)

## compute pairwise matrix of p-values and corresponding colors
pwIT <- pairwiseIndepTest(cu.u, verbose=TRUE) # (d,d)-matrix of test results
round(pmat <- pviTest(pwIT), 3) # pick out p-values
cc <- pairsColList(pmat) # compute corresponding colors

## which pairs violate H0?
which(pmat < 0.05, arr.ind=TRUE)

## pairwise Rosenblatt plot
title <- list(paste(pwRoto, "to test"),
              substitute(italic(H[0]:C~~bold("is nested Gumbel with"~~
                                                      tau[0]==tau0*","
                                                      ~~tau[1]==tau1*","
                                                      ~~tau[2]==tau2)),
                         list(tau0=tau0[1], tau1=tau0[2], tau2=tau0[3])))
main.line <- c(4, 1.4)
## --- JCGS, Fig.1(right) ---
pairsRosenblatt(cu.u, pvalueMat=pmat, pch=".", main=title, main.line=main.line,
                main.centered=TRUE)


### Example 3: 5d t_4 copula (fixed/known d.o.f., estimated P) #################

## setup
n <- 1000 # sample size
d <- 5 # dimension
family <- "t" # copula family
df <- 4 # degrees of freedom
if(setSeeds) set.seed(4)

## define and sample the copula, build pseudo-observations
tau <- c(0.2, 0.4, 0.6)
r <- tauInv(tau, family=family)
P <- c(r[2], r[1], r[1], r[1], # upper triangle (without diagonal) of correlation "matrix"
             r[1], r[1], r[1],
                   r[3], r[3],
                         r[3])
copt4 <- copCreate(family, theta=P, d=d, dispstr="un", df=df, df.fixed=TRUE)
U <- rcop(n, cop=copt4)
U. <- pobs(U)

## define the H0 copula
## Note: that's the same result when using pseudo-observations since estimation via
##       tau is invariant strictly increasing transformations
P. <- nearPD(tauInv(cor(U., method="kendall"), family=family))$mat # estimate P
P. <- P.[lower.tri(P.)] # note: upper.tri() would lead to the wrong result due to the required ordering
plot(P, P., asp=1); abline(0,1, col=adjustcolor("gray",0.5)) # P. should be close to P
copH0 <- copCreate(family, theta=P., d=d, dispstr="un", df=df, df.fixed=TRUE)

## create array of pairwise copH0-transformed data columns
cu.u <- pairwiseCcop(U., copH0, df=df)

## compute pairwise matrix of p-values and corresponding colors
pwIT <- pairwiseIndepTest(cu.u, verbose=TRUE) # (d,d)-matrix of test results
round(pmat <- pviTest(pwIT), 3) # pick out p-values
cc <- pairsColList(pmat) # compute corresponding colors

## which pairs violate H0?
which(pmat < 0.05, arr.ind=TRUE)

## pairwise Rosenblatt plot
title <- list("Pairwise Rosenblatt transformed pseudo-observations",
              expression(bold("to test")~~italic(H[0]:C~~bold("is t")[4])))
main.line <- c(4, 1.4)
## --- JCGS, Fig.2(left) ---
pairsRosenblatt(cu.u, pvalueMat=pmat, pch=".", main=title, main.line=main.line,
                main.centered=TRUE)
## --- JCGS, Fig.2(right) ---
pairsRosenblatt(cu.u, pvalueMat=pmat, pch=".", main=title, main.line=main.line,
                method = "QQchisq", main.centered=TRUE)



### Example 4: SMI constituents ################################################

begSMI <- "2011-09-09"
endSMI <- "2012-03-28"
##-- read *public* data ------------------------------
require(zoo)# -> to access all the zoo methods
if(file.exists(SMI.dfil <- "SMI_yh.rda")) load(SMI.dfil) else {
    require(tseries)
    symSMI <- c("ABBN.VX","ATLN.VX","ADEN.VX","CSGN.VX","GIVN.VX","HOLN.VX",
    	    "BAER.VX","NESN.VX","NOVN.VX","CFR.VX", "ROG.VX", "SGSN.VX",
    	    "UHR.VX", "SREN.VX","SCMN.VX","SYNN.VX","SYST.VX","RIGN.VX",
    	    "UBSN.VX","ZURN.VX")
    lSMI <- sapply(symSMI, function(sym)
    	       get.hist.quote(instrument = sym, start= begSMI, end= endSMI,
    			      quote = "Close", provider = "yahoo",
    			      ## fails: retclass= "its", #-> POSIXct
    			      drop=TRUE))
    ## check if stock data have the same length for each company.
    sapply(lSMI, length)

    ## "concatenate" all:
    SMIo <- do.call(cbind, lSMI)
    ## and fill in the NAs
    SMI <- na.fill(SMIo, "extend")
    colnames(SMI) <- sub("\\.VX", "", colnames(SMI))

    save(symSMI, lSMI, SMI, file=SMI.dfil)
}

n <- nrow(SMI)
d <- ncol(SMI)

x <- diff(log(SMI))[-1,] # build log-returns
u <- pobs(x) # build pseudo-observations

## --- JCGS, Fig.3 ---
pairs(u, gap=0, pch=".", xaxt="n", yaxt="n", main="Pseudo-observations of the log-returns of the SMI",
      labels=as.expression( sapply(1:d, function(j) bquote(italic(hat(U)[.(j)]))) ))

tau <- cor(u, method="kendall") # estimate pairwise tau
P <- tauInv(tau, family="normal") # compute corresponding matrix of pairwise correlations (equal to family="t")

### Estimate (a) t-copula(s) with the approach of Demarta, McNeil (2005)

## compute a positive-definite estimator of P
P. <- nearPD(P, corr=TRUE)$mat
## image(P.) # nice (because 'P.' is a Matrix-pkg Matrix)

## estimate nu via MLE for given P
## nus <- seq(.5, 128, by=.5)
## nLLt.nu <- sapply(nus, nLLt, P=as.matrix(P.), u=u)
## for optimization, use s.th. like:
## optimize(nLLt, interval=c(.5, 128), P=as.matrix(P.), u=u)$minimum
## plot(nus, nLLt.nu, type="l", xlab=bquote(nu),
##      ylab=expression(-logL(nu)))
## plot(nus, nLLt.nu, type="l", xlab=bquote(nu),
##      ylab=expression(-logL(nu)), log = "xy")
## okay, ... points to a Gaussian copula

## define the H0 copula
## Note: that's the same result when using pseudo-observations since estimation via
##       tau is invariant under strictly increasing transformations
P.. <- P.[lower.tri(P.)] # note: upper.tri() would lead to the wrong result due to the required ordering!
copH0 <- copCreate("normal", theta=P.., d=ncol(P.), dispstr="un")

## create array of pairwise copH0-transformed data columns
cu.u <- pairwiseCcop(u, copH0)

if(setSeeds) set.seed(8)
## compute pairwise matrix of p-values and corresponding colors
if(file.exists(pwfil <- "pwIT5.rda")) load(pwfil) else {
  pwIT <- pairwiseIndepTest(cu.u, verbose=TRUE) # (d,d)-matrix of test results
  save(pwIT, file=pwfil)
}
pmat <- pviTest(pwIT) # pick out p-values
cc <- pairsColList(pmat) # compute corresponding colors

## which pairs violate H0?
## (ind <- which(pmat < 0.05, arr.ind=TRUE)) # => none!

{## testing *multivariate normality*
    require(mvnormtest)

    print( mshapiro.test(t(x)) )
    ## => p-value < 2.2e-16
    ## => Multivariate normal dependence structure seems fine, but
    ##  *not* multivariate normal distribution
}

## pairwise Rosenblatt plot --- JCGS, Fig.4 ---
title <- list("Pairwise Rosenblatt transformed pseudo-observations",
              expression(bold("to test")~~italic(H[0]:C~~bold("is Gaussian"))))
pairsRosenblatt(cu.u, pvalueMat=pmat, method="none", cex.labels=0.75,
                key.space=1.5,
                main.centered=TRUE, main=title, main.line=c(3, 0.4))
