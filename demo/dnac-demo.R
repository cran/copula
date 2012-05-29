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


### Densities of two-level nested Archimedean copulas ##########################

## Note: see http://arxiv.org/abs/1204.2410 for more details


### Setup ######################################################################

(srcFn <- system.file("Rsource", "dnac.R", package="copula", mustWork = TRUE))
source(srcFn)#--> a.coeff(),  b.coeff()  nacLL()

require(copula)
require(grid)
require(lattice)
## NB:  nacLL() will use package  partitions  and  polynom

### Example 1: ((1,2), (3,4,5))-Gumbel #########################################

## setup
n <- 250
family <- "Gumbel"
cop. <- getAcop(family)
tau <- c(0.2, 0.4, 0.6)
th <- cop.@tauInv(tau)
cop <- onacopulaL(family, list(th[1], NULL, list(list(th[2], 1:2),
                                                 list(th[3], 3:5))))

## sample and compute log-likelihood
U <- rnacopula(n, cop)
nacLL(cop, u=U) # log-likelihood at correct parameters


### Example 2: (1, (2,3), 4, (5,6,7)) Gumbel ###################################

## setup
n <- 250
family <- "Gumbel"
cop. <- getAcop(family)
tau <- c(0.2, 0.4, 0.6)
th <- cop.@tauInv(tau)
cop <- onacopulaL(family, list(th[1], c(1,4), list(list(th[2], 2:3), list(th[3],
                                                                          5:7))))

## Sample and compute log-likelihood
U <- rnacopula(n, cop)
nacLL(cop, u=U) # log-likelihood at correct parameters


### Example 3: (1, (2,3)) Gumbel ###############################################

## setup
n <- 250
family <- "Gumbel"
tau <- c(0.25, 0.5)
th <- getAcop(family)@tauInv(tau)
copTrue <- onacopulaL(family, list(th[1], 1, list(list(th[2], 2:3))))

## Sample and compute log-likelihood
U <- rnacopula(n, copTrue)
nacLL(copTrue, u=U) # log-likelihood at correct parameters

rm(list = ls(patt="^..?.?$"))# cleanup

### -log-Likelihood plots ######################################################

### (1) (1, (2,...)) structure #################################################

## setup: use this to adjust the family and scomp
family <- "Gumbel" # "Clayton" or "Gumbel"
## "*Tr" = true {nesting structure, model, parameters} :
 compTr  <- 1   ## non-sectorial indices
scompTr <- 2:3  ## sectorial indices   --- for plotting, need  '2:d'
## The following restriction is "just for plot labeling"
stopifnot(compTr==1) # otherwise, sub (for plotting, see below) is wrong

## determine values
n <- 100
cop. <- getAcop(family)
tau <- c(0.25, 0.5)
(thTr <- cop.@tauInv(tau))
cop <- onacopulaL(family, list(thTr[1], compTr,
                               list(list(thTr[2], scompTr)))) # copula
U <- rnacopula(n, cop) # sample
h <- 0.2 # delta{tau} for defining a range of theta's
(th0 <- cop.@tauInv(c(tau[1]-h, tau[1]+h)))
(th1 <- cop.@tauInv(c(tau[2]-h, tau[2]+h)))
m <- 20 # number of grid points
grid <- expand.grid(th0= seq(th0[1], th0[2], length.out=m),
                    th1= seq(th1[1], th1[2], length.out=m))
##' negative Log Likelihood for two-parameter case
##'  C_0({u_j}, C_1( {u_k})) where j in 'comp';  k in 'scomp'
nLL2 <- function(th, u, family, comp, scomp) {
    stopifnot(length(th) == 2)
    if(th[1] > th[2]) ## not fulfilling sufficient nesting condition
	return(Inf)# for minimization better than NA
    cop <- onacopulaL(family, list(th[1], comp, list(list(th[2], scomp))))
    -nacLL(cop, u=u)
}
val.grid <- apply(grid, 1, nLL2,
		  u=U, family=family, comp=compTr, scomp=scompTr)

## determine plot supplements
true.val <- c(th0=thTr[1], th1=thTr[2],
              nLL=nLL2(thTr, u=U, family=family, comp=compTr, scomp=scompTr)) # true value
ind <- which.min(val.grid)
opt.val <- c(grid[ind,], nLL=val.grid[ind]) # optimum on the grid
pts <- rbind(true.val, opt.val) # points to add to wireframe plot
title <- paste("-log-likelihood of a nested", family, "copula") # title
mysec <- { if(length(scompTr)==2) bquote(italic(u[3]))
           else substitute(list(...,italic(u[j])), list(j=max(scompTr))) }
sub <- substitute(italic(C(bolditalic(u)))==italic(C[0](u[1],C[1](u[2],MSEC))) ~~~~~~
                  italic(n)==N ~~~~~~ tau(theta[0])==TAU0 ~~~~~~ tau(theta[1])==TAU1,
                  list(MSEC=mysec, N=n, TAU0=tau[1], TAU1=tau[2]))
sub <- as.expression(sub) # lattice "bug" (only needed by lattice)
xlab <- expression(italic(theta[0]))
ylab <- expression(italic(theta[1]))
zlab <- list(as.expression(-log~L *
    group("(",italic(theta[0])*"," ~ italic(theta[1])*";"~bolditalic(u),")")),
             rot = 90)
sTit <- list(c(expression(group("(",list(theta[0],theta[1]),")")^T),
	       expression(group("(",list(hat(theta)["0,n"],hat(theta)["1,n"]),")")^T)))

## wireframe plot
(wf <- wireframe(val.grid~grid[,1]*grid[,2], aspect=1, zoom=1.02, xlim=th0, ylim=th1,
                 zlim= range(val.grid, as.numeric(pts[,3]), finite=TRUE),
                 xlab=xlab, ylab=ylab, zlab=zlab, main=title, sub=sub, pts=pts,
                 par.settings=list(standard.theme(color=FALSE),
                 layout.heights=list(sub=2.4), background=list(col="#ffffff00"),
                 axis.line=list(col="transparent"), clip=list(panel="off")),
                 alpha.regions=0.5, scales=list(col=1, arrows=FALSE),
                 ## add wire/points
                 panel.3d.wireframe = function(x, y, z, xlim, ylim, zlim, xlim.scaled,
                 ylim.scaled, zlim.scaled, pts, ...) {
                     panel.3dwire(x=x, y=y, z=z, xlim=xlim, ylim=ylim, zlim=zlim,
                                  xlim.scaled=xlim.scaled, ylim.scaled=ylim.scaled,
                                  zlim.scaled=zlim.scaled, ...)
                     panel.3dscatter(x=as.numeric(pts[,1]), y=as.numeric(pts[,2]),
                                     z=as.numeric(pts[,3]), xlim=xlim, ylim=ylim, zlim=zlim,
                                     xlim.scaled=xlim.scaled, ylim.scaled=ylim.scaled,
                                     zlim.scaled=zlim.scaled, type="p", col=1,
                                     pch=c(3,4), lex=2, cex=1.4, .scale=TRUE, ...)
                 },
		 key =
		 list(x=-0.01, y=1, points=list(pch=c(3,4), col=1, lwd=2, cex=1.4),
		      text = sTit, padding.text=3, cex=1, align=TRUE, transparent=TRUE)) )

## levelplot
(lp <- levelplot(val.grid~grid[,1]*grid[,2], aspect=1, xlab=xlab, ylab=ylab,
                 par.settings=list(layout.heights=list(main=3, sub=2),
                 regions=list(col=gray(140:400/400))),
                 xlim= extendrange(grid[,1], f = 0.04),
                 ylim= extendrange(grid[,2], f = 0.04),
                 main=title, sub=sub, pts=pts,
                 scales=list(alternating=c(1,1), tck=c(1,0)), contour=TRUE,
                 panel=function(x, y, z, pts, ...){
                     panel.levelplot(x=x, y=y, z=z, ...)
                     grid.points(x=pts[1,1], y=pts[1,2], pch=3,
                                 gp=gpar(lwd=2, col="black")) # + true value
                     grid.points(x=pts[2,1], y=pts[2,2], pch=4,
                                 gp=gpar(lwd=2, col="black")) # x optimum
                 },
		 key =
		 list(x=0.18, y=1.09, points=list(pch=c(3,4), col=1, lwd=2, cex=1.4),
		      columns = 2, text = sTit, align=TRUE, transparent=TRUE)) )


### ---------- Maximum Likelihood estimation -------------------------------------
### ---------- via optimization instead of grid ----------------------------------
## starting not too closely
ropt <- optim(c(1, 3), nLL2,
              u=U, family=family, comp=compTr, scomp=scompTr#, control=list(trace=1)
              )
## compare with previous "grid result" -- it is a bit better
rbind(pts, optim= c(ropt$par, ropt$value))
