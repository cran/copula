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


### Demo of the two-parameter Archimedean GIG copula ###########################

### Setup ######################################################################

source(system.file("Rsource", "GIG.R", package="copula"))

require(copula)
require(bbmle)
require(lattice)
require(grid)

do.profile <- FALSE # set this to TRUE to compute profile-likelihood plots (time-consuming)
## without the profile-likelihood plots, the demo takes about 5min on a standard laptop


### plot Kendall's tau as a function in theta for different nu's

th1 <- c(0, 0.1, 0.5, 1, 5, 10)
cols <- colorRampPalette(c("red", "orange", "darkgreen", "turquoise", "blue"),
                         space="Lab")(length(th1))
for(i in seq_along(th1))
    curve(tau.GIG(cbind(th1[i],x)), 1e-12, 2,
          main="Kendall's tau for the GIG family", ylim=c(0,1),
          xlab=expression(theta), ylab=expression(tau(nu,theta)), add=(i>1),
          lwd=1.4, col=cols[i])
label <- as.expression(lapply(1:length(th1), function(i) substitute(nu==nu., list(nu.=th1[i]))))
legend("topright", label, bty="n", lwd=1.4, col=cols)


### sampling and estimation ####################################################

## specify parameters
n <- 100 # sample size
d <- 10 # dimension
nu <- 0.2 # fix nu
tau <- 0.5 # => psi(t) = (1+t)^(-nu/2)besselK(theta*sqrt(1+t), nu=nu)/besselK(theta, nu=nu) with (nu,theta) = (0.2, 0.0838) (see below)

## adjustment for initial value
h <- c(0.15, 0.15) # h_-, h_+


### functions ##################################################################

##' Initial interval for GIG
##'
##' @title Initial interval for GIG
##' @param U (n x d)-matrix of simulated data
##' @param h non-negative auxiliary parameter for computing initial intervals
##' @param method "etau" via sample version of Kendall's tau
##'               "dmle.G" via DMLE of Gumbel
##' @return (2 x 2)-matrix containing the initial interval [1st row: lower,
##'         2nd row: upper; 2 parameters => 2 cols]
##' @author Marius Hofert
ii.GIG <- function(U, h, method=c("etau","dmle.G")){
    stopifnot(h >= 0, length(h) >= 2)
    I <- matrix(, nrow=2, ncol=2, dimnames=list(c("lower", "upper"), c("nu", "theta")))
    ## estimate Kendall's tau
    method <- match.arg(method)
    tau.hat <- switch(method,
                      "etau" = { # uses sample version of tau, more accurate but slower
                          tau.hat.mat <- cor(U, method="kendall")
                          mean(tau.hat.mat[upper.tri(tau.hat.mat)])
                      },
                      "dmle.G" = { # uses DMLE for Gumbel to estimate tau
                          Z <- apply(U, 1, max)
                          theta.hat.G <- log(ncol(U))/(log(length(Z))-log(sum(-log(Z))))
                          copGumbel@tau(theta.hat.G)
                      },
                      stop("wrong method:", method))
    ## compute largest value of theta (for upper left endpoint of the inital interval)
    stopifnot(tau.hat > 0)
    nu.min <- 0
    I[1,1] <- nu.min # smallest value for nu
    th.max <- iTau.GIG(max(tau.hat-h[1],0.005), theta=c(nu.min, NA))
    I[2,2] <- th.max[2] # largest value for theta
    ## compute smallest theta (for lower left endpoint of the inital interval)
    th.min <- iTau.GIG(min(tau.hat+h[2],0.995), theta=c(nu.min, NA)) # largest attainable tau with 1e-30 is one.m.eps=0.9602
    I[1,2] <- th.min[2]
    ## compute largest nu (for lower right endpoint of the inital interval)
    nu.max <- iTau.GIG(max(tau.hat-h[1],0.005), theta=c(NA, th.min[2]))
    I[2,1] <- nu.max[1]
    ## result
    I
}

##' -log-likelihood
##'
##' @title -log-likelihood
##' @param nu parameter of the generator/copula
##' @param theta parameter of the generator/copula
##' @param u (n x d)-matrix of simulated data
##' @return -sum(log(density))
##' @author Marius Hofert
nlogl.GIG <- function(nu, theta, u){
    if(!is.matrix(u)) u <- rbind(u)
    if((d <- ncol(u)) < 2) stop("u should be at least bivariate") # check that d >= 2
    -sum(dacopula.GIG(u, theta=c(nu, theta), n.MC=0, log=TRUE))
}

## vectorized version
nlogl.GIG. <- function(theta, u) nlogl.GIG(theta[1], theta=theta[2], u=u)


### estimation #################################################################

## determine theta such that tau is matched (for given nu)
theta <- iTau.GIG(tau, c(nu, NA))

## sample
set.seed(1000)
U <- rnacopula.GIG(n, d, theta)

## plot
splom2(U, cex=0.4, pscales=0, main=paste("Sample of size",n,
                              "from a GIG copula"))

## initial interval and value
I <- ii.GIG(U, h)
start <- colMeans(I)

## without profiling: optim with method="L-BFGS-B"
system.time(optim(par=start, method="L-BFGS-B",
                  fn=function(x) nlogl.GIG(x[1], theta=x[2], u=U),
                  lower=c(I[1,1], I[1,2]), upper=c(I[2,1], I[2,2])))

## with profiling: via mle (uses optimizer="optim" with method="L-BFGS-B")
nLL <- function(nu, theta) nlogl.GIG(nu, theta, u=U)
system.time(ml <- mle(nLL, method="L-BFGS-B",
                      start=list(nu=mean(I[,1]), theta=mean(I[,2])),
                      lower=c(nu=I[1,1], theta=I[1,2]),
                      upper=c(nu=I[2,1], theta=I[2,2])))
summary(ml)
str(ml@details)

## with profiling: via mle2 (uses optimizer="optim" with method="L-BFGS-B")
system.time(ml2 <- mle2(nlogl.GIG, data=list(u=U), method="L-BFGS-B",
                        start=list(nu=mean(I[,1]), theta=mean(I[,2])),
                        lower=c(nu=I[1,1], theta=I[1,2]),
                        upper=c(nu=I[2,1], theta=I[2,2])))
summary(ml2)
str(ml2@details)

## Note: run time is mainly determined by the evaluations of iPsi.GIG


### plots ######################################################################

### profile likelihood plots ###################################################

if(do.profile){

    system.time(prof <- profile(ml))
    if(FALSE) { ## FIXME (?)
        ## maybe this helps: https://stat.ethz.ch/pipermail/r-help/2005-July/076003.html
        ci <- confint(prof)
        ci
        plot(prof)
    }

    system.time(prof2 <- profile(ml2)) # profiling (time-consuming)
    (ci <- confint(prof2))
    plot(prof2) # => for adjusting stepsize etc., see ?profile.mle2

}


### -log-likelihood plots ######################################################

## build grid (time-consuming)
m <- 20 # number of grid points = number of intervals + 1
th <- seq(I[1,1], I[2,1], length.out=m) # grid points for nu
beta <- seq(I[1,2], I[2,2], length.out=m) # grid points for theta
grid <- expand.grid(theta=th, beta=beta) # grid
system.time(val.grid <- apply(grid, 1, nlogl.GIG., u=U)) # value of the -log-likelihood on the grid

### wireframe

## plot settings
true.theta <- theta
true.val <- c(true.theta, nlogl.GIG.(true.theta, u=U)) # theoretical optimum
opt <- ml@coef # optimizer-optimum
opt.val <- c(opt, nlogl.GIG.(opt, u=U)) # optimizer-optimum and its value
pts <- rbind(true.val, opt.val) # points to add to wireframe plot
title <- "-log-likelihood of an Archimedean GIG copula" # title
sub <- substitute(italic(n) == N ~~~  italic(d)== D ~~~
                  tau == TAU ~~~ "#{eval}:" ~ NIT,
                  list(N=n, D=d, TAU= tau, NIT= ml@details$counts[[1]]))
## lattice bug:
sub <- as.expression(sub)

## wireframe
wireframe(val.grid~grid[,1]*grid[,2], screen=list(z=70, x=-55), zoom=0.95,
          xlab=expression(italic(theta)), ylab=expression(italic(beta)),
          zlab = list(as.expression(-log~L * group("(",list(theta,beta),")")), rot=90),
          main=title, sub=sub, pts=pts, scales=list(col=1, arrows=FALSE),
          par.settings=list(axis.line=list(col="transparent"),
          clip=list(panel="off")), zlim=c(min(val.grid, pts[,3]),
                                   max(val.grid, pts[,3])), aspect=1,
          panel.3d.wireframe = function(x,y,z,xlim,ylim,zlim,xlim.scaled,
          ylim.scaled,zlim.scaled,pts,...){
              panel.3dwire(x=x, y=y, z=z, xlim=xlim, ylim=ylim, zlim=zlim,
                           xlim.scaled=xlim.scaled, ylim.scaled=ylim.scaled,
                           zlim.scaled=zlim.scaled, alpha.regions=0.8, ...)
              panel.3dscatter(x=pts[,1], y=pts[,2], z=pts[,3],
                              xlim=xlim, ylim=ylim, zlim=zlim,
                              xlim.scaled=xlim.scaled, ylim.scaled=ylim.scaled,
                              zlim.scaled=zlim.scaled, type="p", col=c("red","blue"),
                              pch=c(3,4), lex=2, cex=1.4, .scale=TRUE, ...)
          }, key=list(x=0.64, y=1.01, points=list(pch=c(3,4), col=c("red","blue"),
                                      lwd=2, cex=1.4),
             text=list(c("True value", "Optimum of optimizer")), padding.text=3,
             cex=1, align=TRUE, transparent=TRUE))

### levelplot

## plot settings
xlim. <- c(min(grid[,1]),max(grid[,1]))
ylim. <- c(min(grid[,2]),max(grid[,2]))
xeps <- (xlim.[2] - xlim.[1]) * 0.04
yeps <- (ylim.[2] - ylim.[1]) * 0.04
cols <- adjustcolor(colorRampPalette(c("darkgreen", "green", "orange", "yellow"),
                                     space="Lab")(100), 0.8)

## levelplot
levelplot(val.grid~grid[,1]*grid[,2],
          par.settings=list(layout.heights=list(main=3, sub=2),
          regions=list(col=cols)),
          xlim=c(xlim.[1]-xeps, xlim.[2]+xeps),
          ylim=c(ylim.[1]-yeps, ylim.[2]+yeps),
          xlab=expression(italic(theta)), ylab=expression(italic(beta)),
          main=title, sub=sub, pts=pts, aspect=1,
          scales=list(alternating=c(1,1), tck=c(1,0)), contour=TRUE,
          panel=function(x, y, z, pts, ...){
              panel.levelplot(x=x, y=y, z=z, ...)
              grid.points(x=pts[1,1], y=pts[1,2], pch=3,
                          gp=gpar(lwd=2, col="red")) # + true value
              grid.points(x=pts[2,1], y=pts[2,2], pch=4,
                          gp=gpar(lwd=2, col="blue")) # x optimum
          }, key=list(x=0.18, y=1.08, points=list(pch=c(3,4), col=c("red","blue"),
                                      lwd=2, cex=1.4),
             columns=2, text=list(c("True value", "Optimum of optimizer")),
             align=TRUE, transparent=TRUE))
