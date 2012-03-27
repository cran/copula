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


strictify <- function(val,status)
  {
    val[status>0] <- NaN
    return(val)
  }


debye1 <- function(x, give=FALSE, strict=TRUE){
  attr <- attributes(x)
  x.vec <- as.vector(x)
  jj <- .C(debye_1_C,
           as.double(abs(x.vec)),  ## added abs by JY
           as.integer(length(x.vec)),
           val=as.double(x.vec),
           err=as.double(x.vec),
           status=as.integer(0*x.vec))

  val <- ifelse(x.vec >=0, jj$val, jj$val - x.vec / 2) ## k = 1, Frees & Valdez 1998, p.9
  err <- jj$err
  status <- jj$status
  attributes(val) <- attr
  attributes(err) <- attr
  attributes(status) <- attr

  if(strict){
    val <- strictify(val,status)
  }

  if(give){
      return(list(val=val,err=err,status=status))
  } else {
    return(val)
  }
}

debye2 <- function(x, give=FALSE, strict=TRUE){
  attr <- attributes(x)
  x.vec <- as.vector(x)
  jj <- .C(debye_2,
           as.double(abs(x.vec)),  ## added abs by JY
           as.integer(length(x.vec)),
           val=as.double(x.vec),
           err=as.double(x.vec),
           status=as.integer(0*x.vec))
  val <- ifelse(x.vec >= 0, jj$val, jj$val - x.vec * 2/ 3) ## k = 2
  err <- jj$err
  status <- jj$status
  attributes(val) <- attr
  attributes(err) <- attr
  attributes(status) <- attr


  if(strict){
    val <- strictify(val,status)
  }

  if(give){
      return(list(val=val,err=err,status=status))
  } else {
    return(val)
  }
}

## debye3 <- function(x, give=FALSE, strict=TRUE){
##   attr <- attributes(x)
##   x.vec <- as.vector(x)
##   jj <- .C(debye_3,
##            as.double(x.vec),
##            as.integer(length(x.vec)),
##            val=as.double(x.vec),
##            err=as.double(x.vec),
##            status=as.integer(0*x.vec))
##   val <- jj$val
##   err <- jj$err
##   status <- jj$status
##   attributes(val) <- attr
##   attributes(err) <- attr
##   attributes(status) <- attr

##   if(strict){
##     val <- strictify(val,status)
##   }

##   if(give){
##       return(list(val=val,err=err,status=status))
##   } else {
##     return(val)
##   }
## }

## debye4 <- function(x, give=FALSE, strict=TRUE){
##   attr <- attributes(x)
##   x.vec <- as.vector(x)
##   jj <- .C("debye_4",
##            as.double(x.vec),
##            as.integer(length(x.vec)),
##            val=as.double(x.vec),
##            err=as.double(x.vec),
##            status=as.integer(0*x.vec),
##            PACKAGE="copula"
##            )
##   val <- jj$val
##   err <- jj$err
##   status <- jj$status
##   attributes(val) <- attr
##   attributes(err) <- attr
##   attributes(status) <- attr

##   if(strict){
##     val <- strictify(val,status)
##   }

##   if(give){
##       return(list(val=val,err=err,status=status))
##   } else {
##     return(val)
##   }
## }

## ## debye function is used for compute the assoc measure of frankCopula
## debye <- function(x, k, ...) {
##   Dk.integrand <- function(t) t^k / (exp(t) - 1)
##   Dk.int <- function(x, k, ...) {
##     integrate(Dk.integrand, 0, x, ...)$value
##   }
##   y <- abs(x)
##   Dk <- k / y^k * sapply(y, Dk.int, k = k, ...)
##   ifelse(x < 0, Dk <- Dk + k * y / (k + 1), Dk)
## }
