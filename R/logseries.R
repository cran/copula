##### Kemp 1981, Applied Statistics 30(3), pp349--253 for RNG
##### The cf is Log[1 - a exp(i t)] / Log[1 - a]
dlogseries <- function(x, alpha, log = FALSE) {
  val <-  - alpha^x / x / log(1 - alpha)
  if (log) log(val) else val
}


plogseries <- function(q, alpha, lower.tail = TRUE, log.p = FALSE) {
  if (!lower.tail) q <- 1 - q
  val <- 1 - alpha^(1 + q) * hyperg2F1(1 + q, 1, 2 + q, alpha) / (1 + q)
  if (log.p) log(val) else val
}

rlogseries <- function(n, alpha) {
  val <- integer(n)
  alpha <- rep(alpha, len = n)
  .C("rlogseries_R", as.integer(n), as.double(alpha), val = as.integer(val),
     PACKAGE = "copula")$val
}
