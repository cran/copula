# .First.lib <- function(lib, pkg) {
#   require(mvtnorm)
#   require(sn)
#   library.dynam("copula", pkg, lib)
# }

.onLoad <- function(lib, pkg) {
  require(mvtnorm)
  require(sn)
  ##require(scatterplot3d)
  require(methods)
  ##library.dynam("copula", pkg, lib)
}
