perspCopula <- function(x, fun, n = 51, theta = -30, phi = 30, expand = 0.618, ...) {
  eps <- (.Machine$double.eps)^(1/4)
  eps <- 0
  xis <- yis <- seq(0 + eps, 1 - eps, len = n)
  grids <- as.matrix(expand.grid(xis, yis))
  zmat <- matrix(fun(x, grids), n, n)
  persp(xis, yis, zmat, theta = theta, phi = phi, expand = expand, ...)
  val <- list(x = xis, y = yis, z = zmat)
  invisible(val)
}



contourCopula <- function(x, fun, n = 51,...) {
  eps <- (.Machine$double.eps)^(1/4)
  eps <- 0
  xis <- yis <- seq(0 + eps, 1 - eps, len = n)
  grids <- as.matrix(expand.grid(xis, yis))
  zmat <- matrix(fun(x, grids), n, n)
  contour(xis, yis, zmat, ...)
  val <- list(x = xis, y = yis, z = zmat)
  invisible(val)
}


perspMvdc <- function(x, fun,
                      xis, yis,
                      theta = -30, phi = 30, expand = 0.618, ...) {
  n <- length(xis)
  grids <- as.matrix(expand.grid(xis, yis))
  zmat <- matrix(fun(x, grids), n, n)
  persp(xis, yis, zmat, theta = theta, phi = phi, expand = expand, ...)
  val <- list(x = xis, y = yis, z = zmat)
  invisible(val)
}



contourMvdc <- function(x, fun,
                       xis, yis,...) {
  n <- length(xis)
  grids <- as.matrix(expand.grid(xis, yis))
  zmat <- matrix(fun(x, grids), n, n)
  contour(xis, yis, zmat, ...)
  val <- list(x = xis, y = yis, z = zmat)
  invisible(val)
}

setMethod("persp", signature("copula"), perspCopula)
setMethod("contour", signature("copula"), contourCopula)

setMethod("persp", signature("mvdc"), perspMvdc)
setMethod("contour", signature("mvdc"), contourMvdc)
