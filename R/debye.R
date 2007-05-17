strictify <- function(val,status)
  {
    val[status>0] <- NaN
    return(val)
  }


debye1 <- function(x, give=FALSE, strict=TRUE){
  attr <- attributes(x)
  x.vec <- as.vector(x)
  jj <- .C("debye_1",
           as.double(abs(x.vec)),  ## added abs by JY
           as.integer(length(x.vec)),
           val=as.double(x.vec),
           err=as.double(x.vec),
           status=as.integer(0*x.vec),
           PACKAGE="copula"
           )
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
  jj <- .C("debye_2",
           as.double(abs(x.vec)),  ## added abs by JY
           as.integer(length(x.vec)),
           val=as.double(x.vec),
           err=as.double(x.vec),
           status=as.integer(0*x.vec),
           PACKAGE="copula"
           )
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
##   jj <- .C("debye_3",
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
