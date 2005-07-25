setClass("claytonCopula",
         representation = representation("archmCopula"),
#         validity =  validCopula,
         contains = list("copula", "archmCopula")
         )


genFunClayton <- function(copula, u) {
  alpha <- copula@parameters[1]
  u^(-alpha) - 1
}

genInvClayton <- function(copula, s) {
  alpha <- copula@parameters[1]
  (1 + s)^(-1/alpha)
}

claytonCopula <- function(param, dim = 2) {
  ## get expressions of cdf and pdf
  cdfExpr <- function(n) {
    expr <- "u1^(-alpha) - 1"
    for (i in 2:n) {
      cur <- paste( "u", i, "^(-alpha) - 1", sep = "")
      expr <- paste(expr, cur, sep=" + ")
    }
    expr <- paste("(1 + (", expr, "))^ (-1/alpha)")
    parse(text = expr)
  }
  
  pdfExpr <- function(cdf, n) {
    val <- cdf
    for (i in 1:n) {
      val <- D(val, paste("u", i, sep=""))
    }
    val
  }

  if (dim > 2 && param[1] < 0)
    stop("param can be negative only for dim = 2")
  cdf <- cdfExpr(dim)
  pdf <- pdfExpr(cdf, dim)
  val <- new("claytonCopula",
             dimension = dim,
             parameters = param[1],
             exprdist = c(cdf = cdf, pdf = pdf),
             param.names = "param",
             param.lowbnd = -1,
             param.upbnd = Inf,
             message = "Clayton copula family; Archimedean copula")
  val
}

rclaytonBivCopula <- function(copula, n) {
  val <- cbind(runif(n), runif(n))
  alpha <- copula@parameters[1]
  val[,2] <- (val[,1]^(-alpha) * (val[,2]^(-alpha/(alpha + 1)) - 1) + 1)^(-1/alpha)
  val
}


rclaytonCopula <- function(copula, n) {
  dim <- copula@dimension
  if (dim == 2) return (rclaytonBivCopula(copula, n))
  ## gamma frailty
  alpha <- copula@parameters[1]
  val <- matrix(runif(n * dim), nrow = n)
  if (abs(alpha) <= 100 * .Machine$double.eps)
    return (val)  ## the limit is independence
  gam <- rgamma(n, shape = 1/alpha , rate = 1)
  gam <- matrix(gam, nrow = n, ncol = dim)
  genInv(copula, - log(val) / gam)
}


# pclaytonCopula <- function(copula, u) {
#   dim <- copula@dimension
#   if (is.vector(u)) u <- matrix(u, ncol = dim)
#   alpha <- copula@parameters[1]
#   if (abs(alpha) <= 100 * .Machine$double.eps) return (apply(u, 1, prod))
#   myfun.vector <- function(x) {
#     sum(x^(-alpha) - 1)
#   }
#   (1 + apply(u, 1, myfun.vector))^(-1/alpha)
# }


pclaytonCopula <- function(copula, u) {
  dim <- copula@dimension
  if (is.vector(u)) u <- matrix(u, ncol = dim)
  alpha <- copula@parameters[1]
  if (abs(alpha) <= 100 * .Machine$double.eps) return (apply(u, 1, prod))
  cdf <- copula@exprdist$cdf
  dim <- copula@dimension
  for (i in 1:dim) assign(paste("u", i, sep=""), u[,i])
  val <- eval(cdf)
  pmax(val, 0)
}


dclaytonCopula <- function(copula, u) {
  dim <- copula@dimension
  if (is.vector(u)) u <- matrix(u, ncol = dim)
  alpha <- copula@parameters[1]
  if (abs(alpha) <= 100 * .Machine$double.eps) return (rep(1, length(u)))
  pdf <- copula@exprdist$pdf
  dim <- copula@dimension
  for (i in 1:dim) assign(paste("u", i, sep=""), u[,i])
  val <- eval(pdf)
  val[apply(u, 1, function(v) any(v <= 0))] <- 0
  val[apply(u, 1, function(v) any(v >= 1))] <- 0
  if (alpha < 0) {
    cdf <- pcopula(copula, u)
    bad <- cdf == 0
    val[bad] <- 0
  }
  val
}

setMethod("rcopula", signature("claytonCopula"), rclaytonCopula)
setMethod("pcopula", signature("claytonCopula"), pclaytonCopula)
setMethod("dcopula", signature("claytonCopula"), dclaytonCopula)

setMethod("genFun", signature("claytonCopula"), genFunClayton)
setMethod("genInv", signature("claytonCopula"), genInvClayton)

