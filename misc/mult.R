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


multcomp <- function(u, g, N = 1000, der = 1, multi = 1)
{
  p <- ncol(u)
  n <- nrow(u)
  m <- nrow(g)

  out <- .C("mult",
            as.double(u),
            as.integer(n),
            as.integer(p),
            as.double(g),
            as.integer(m),
            as.integer(N),
            s0 = double(N * m),
            s1 = double(N * m),
            as.integer(der),
            as.integer(multi),
            PACKAGE="copula")

  s0 <- matrix(out$s0, ncol = m, byrow = TRUE)
  s1 <- matrix(out$s1, ncol = m, byrow = TRUE)
  list(s0=s0,s1=s1)
}


