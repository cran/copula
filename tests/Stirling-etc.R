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
(isLinux <- identical("Linux", Sys.info()[["sysname"]]))

### Stirling numbers of the 1st kind ###########################################

S1.10 <- c(0, -362880, 1026576, -1172700, 723680,
           -269325, 63273, -9450, 870, -45, 1)
stopifnot(sapply(0:10, Stirling1, n=10) == S1.10,
          Stirling1.all(10) == S1.10[-1])

options(str = strOptions(vec.len = 10, digits.d = 20)) # for ls.str() below

ls.str(copula:::.nacopEnv)
system.time(S  <- Stirling1(30, 7))# updating table -> typically not zero
system.time(S. <- Stirling1(30, 7))# lookup  -->  should be zero
stopifnot(identical(S, S.))

ls.str(copula:::.nacopEnv)

system.time(s1c <- Stirling1(100,10))
s1c
(s1 <- system.time(for(i in 1:20) S. <- Stirling1(100, 10))[[1]])
stopifnot(identical(S., s1c), !isLinux || s1 <= 0.020)
system.time(s2c <- Stirling1(200,190)); s2c
(s2 <- system.time(for(i in 1:20) S. <- Stirling1(200,190))[[1]])
stopifnot(identical(S., s2c), !isLinux || s2 <= 0.020)
## 0.010 occasionally barely fails (prints "0.010") on Martin's X201


### Stirling numbers of the 2nd kind ###########################################

S2.10 <- c(0, 1, 511, 9330, 34105, 42525, 22827, 5880, 750, 45, 1)
stopifnot(sapply(0:10, Stirling2, n=10, method="direct") == S2.10,
          sapply(0:10, Stirling2, n=10, method="lookup") == S2.10,
          Stirling2.all(10) == S2.10[-1])

ls.str(copula:::.nacopEnv)
system.time(S  <- Stirling2(30, 7))# updating table -> typically not zero
system.time(S. <- Stirling2(30, 7))# lookup  -->  should be zero
stopifnot(identical(S, S.),
          all.equal(S, Stirling2(30,7, method="direct"), tol=1e-15))

ls.str(copula:::.nacopEnv)

rbind(C.direct = system.time(Sd <- Stirling2(100,10, method="direct")),
      C.lookup = system.time(Sl <- Stirling2(100,10, method="lookup")))
## should be equal; and lookup time should be "zero" when called again:
(s3 <- system.time(for(i in 1:20) S. <- Stirling2(100, 10))[[1]])
stopifnot(all.equal(Sd, Sl, tol = 1e-15), !isLinux || s3 <= 0.020)
## 0.010 fails on good ole' Solaris when that is busy..
## Here, the direct method already overflows, but the "lookup" still works
rbind(C.direct = system.time(Sd <- Stirling2(200,190, method="direct")),
      C.lookup = system.time(Sl <- Stirling2(200,190, method="lookup")))
Sd ; Sl
(s4 <- system.time(for(i in 1:20) S. <- Stirling2(200,190))[[1]])
stopifnot(!isLinux || s4 <= 0.025)
# 0.010 occasionally barely fails (prints "0.010") on Martin's X201


### Eulerian Numbers ###########################################################

##' cheap "direct" version of Eulerian.all():
Euleri.A <- function(n)
    sapply(0:max(0,n-1), Eulerian, n=n, method="direct")
stopifnot(identical(Euler.l5 <- lapply(0:5, Euleri.A),
		    list(1,
			 1,
			 c(1, 1),
			 c(1, 4, 1),
			 c(1, 11, 11, 1),
			 c(1, 26, 66, 26, 1))))

p.Eul <- function(n) {
    plot(E1 <- Eulerian.all(n), log="y", yaxt="n",
         xlab = "k", ylab = bquote(A(.(n), k)),
         main = bquote("Eulerian numbers "* A(.(n), k)))
    if(require("sfsmisc"))
	eaxis(2, quantile(axTicks(2), (0:16)/16, type=3), at.small=numeric())
    else axis(2)
    lines(E2 <- Euleri.A(n), col="green3", type="o")
    invisible(cbind(E1=E1, E2=E2))
}

if(!dev.interactive(orNone=TRUE)) pdf("Eulerian-ex.pdf")

e60 <- p.Eul(60); all.equal(e60[,2],e60[,1], tol=0) ## 3.82e-09
e70 <- p.Eul(70); all.equal(e70[,2],e70[,1])      ## 2.97e-6
e90 <- p.Eul(90); all.equal(e90[,2],e90[,1])      ## 0.032
e100 <- p.Eul(100); all.equal(e100[,2],e100[,1])  ## 0.80028 --- visible in center
e110 <- p.Eul(110); all.equal(e110[,2],e110[,1])  ## 0.992   --- visible in center
e120 <- p.Eul(120); all.equal(e120[,2],e120[,1])  ## 1 -- problem in center
e150 <- p.Eul(150) ## clear problem in center -- close to overflow though
e170 <- p.Eul(170) ## clear problem in center -- close to overflow though
max(e170[,"E1"]) # 7.5964e+305 -- almost maximum


### Polylogarithm Function #####################################################

EQ <- function(x,y, tol = 1e-15) all.equal(x,y, tol=tol)

x <- (0:127)/128 # < 1
stopifnot(EQ(polylog(s =  1,  x, n.sum=10000), -log(1-x)),
	  EQ(polylog(s = -1, .1, n.sum=	 100), 10/81),
	  EQ(polylog(s = -1, .1, "negI-s-Stirling"),       10/81),
	  EQ(polylog(x, -1, "negI-s-Stirling"), x /(1-x)^2),
	  EQ(polylog(x, -2, "negI-s-Stirling"), x*(1+x)/(1-x)^3),
	  EQ(polylog(x, -4, "negI-s-Stirling"), x*(1+x)*(1+x*(10+x)) / (1-x)^5),
	  identical(	      polylog	  (x, -4, "negI-s-Stirling"),
		    Vectorize(polylog,"z")(x, -4, "negI-s-Stirling")),
	  identical(	      polylog	  (x, -4, "sum", n.sum=10000),
		    Vectorize(polylog,"z")(x, -4, "sum", n.sum=10000)),
	  EQ(polylog(x, -1, "negI-s-Eulerian"), x /(1-x)^2),
	  EQ(polylog(x, -2, "negI-s-Eulerian"), x*(1+x)/(1-x)^3),
	  EQ(polylog(x, -4, "negI-s-Eulerian"), x*(1+x)*(1+x*(10+x)) / (1-x)^5),
          TRUE)

##--> now do plots etc in  ../man/polylog.Rd :
##                         ~~~~~~~~~~~~~~~~~
