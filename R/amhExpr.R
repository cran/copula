#################################################################################
##
##   R package Copula by Jun Yan Copyright (C) 2008
##
##   This file is part of the R package copula.
##
##   The R package copula is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package copula is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
##   You should have received a copy of the GNU General Public License
##   along with the R package copula. If not, see <http://www.gnu.org/licenses/>.
##
#################################################################################


`amhCopula.pdf.expr` <-
expression(1, -((1 + alpha^2 * (-1 + u1) * (-1 + u2) + alpha * 
    (-2 + u1 + u2 + u1 * u2))/(-1 + alpha * (-1 + u1) * (-1 + 
    u2))^3))
`amhCopula.pdf.algr` <-
expression({
    .value <- 1
    .grad <- array(0, c(length(.value), 1), list(NULL, c("s")))
    .grad[, "s"] <- 0
    attr(.value, "gradient") <- .grad
    .value
}, {
    .expr2 <- -1
    .expr3 <- .expr2 + u1
    .expr5 <- .expr2 + u2
    .value <- -((1 + alpha^2 * .expr3 * .expr5 + alpha * (-2 + 
        u1 + u2 + u1 * u2))/(.expr2 + alpha * .expr3 * .expr5)^3)
    .grad <- array(0, c(length(.value), 1), list(NULL, c("s")))
    .grad[, "s"] <- 0
    attr(.value, "gradient") <- .grad
    .value
})
`amhCopula.genfun.expr` <-
expression((-1 + alpha)/((1 + alpha * (-1 + u)) * u), -(((-1 + 
    alpha) * (1 + alpha * (-1 + 2 * u)))/((1 + alpha * (-1 + 
    u))^2 * u^2)))
`amhCopula.genfun.algr` <-
expression({
    .expr1 <- -1
    .expr2 <- .expr1 + alpha
    .expr5 <- 1 + alpha * (.expr1 + u)
    .expr6 <- .expr5 * u
    .value <- .expr2/.expr6
    .grad <- array(0, c(length(.value), 1), list(NULL, c("u")))
    .grad[, "u"] <- -(.expr2 * (alpha * u + .expr5)/.expr6^2)
    attr(.value, "gradient") <- .grad
    .value
}, {
    .expr1 <- -1
    .expr2 <- .expr1 + alpha
    .expr3 <- 2 * u
    .expr7 <- .expr2 * (1 + alpha * (.expr1 + .expr3))
    .expr10 <- 1 + alpha * (.expr1 + u)
    .expr11 <- .expr10^2
    .expr12 <- u^2
    .expr13 <- .expr11 * .expr12
    .value <- -(.expr7/.expr13)
    .grad <- array(0, c(length(.value), 1), list(NULL, c("u")))
    .grad[, "u"] <- -(.expr2 * (alpha * 2)/.expr13 - .expr7 * 
        (2 * (alpha * .expr10) * .expr12 + .expr11 * .expr3)/.expr13^2)
    attr(.value, "gradient") <- .grad
    .value
})
