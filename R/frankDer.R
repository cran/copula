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


`genFunDerFrank.expr` <-
expression({
    .expr1 <- -1
    .value <- -log((.expr1 + E^-(alpha * u))/(.expr1 + E^-alpha))
    .grad <- array(0, c(length(.value), 1), list(NULL, c("s")))
    .grad[, "s"] <- 0
    attr(.value, "gradient") <- .grad
    .value
}, {
    .value <- alpha/(1 - E^(alpha * u))
    .grad <- array(0, c(length(.value), 1), list(NULL, c("s")))
    .grad[, "s"] <- 0
    attr(.value, "gradient") <- .grad
    .value
}, {
    .expr3 <- E^(alpha * u)
    .value <- alpha^2 * .expr3/(-1 + .expr3)^2
    .grad <- array(0, c(length(.value), 1), list(NULL, c("s")))
    .grad[, "s"] <- 0
    attr(.value, "gradient") <- .grad
    .value
})
`genInvDerFrank.expr` <-
expression({
    .expr4 <- -1 + E^-alpha
    .expr5 <- E^s
    .expr7 <- 1 + .expr4/.expr5
    .value <- -(log(.expr7)/alpha)
    .grad <- array(0, c(length(.value), 1), list(NULL, c("s")))
    .grad[, "s"] <- .expr4 * (.expr5 * log(E))/.expr5^2/.expr7/alpha
    attr(.value, "gradient") <- .grad
    .value
}, {
    .expr1 <- E^alpha
    .expr2 <- 1 - .expr1
    .expr6 <- E^(alpha + s)
    .expr8 <- alpha - alpha * .expr1 + alpha * .expr6
    .value <- .expr2/.expr8
    .grad <- array(0, c(length(.value), 1), list(NULL, c("s")))
    .grad[, "s"] <- -(.expr2 * (alpha * (.expr6 * log(E)))/.expr8^2)
    attr(.value, "gradient") <- .grad
    .value
}, {
    .expr2 <- E^(alpha + s)
    .expr4 <- E^alpha
    .expr5 <- -1 + .expr4
    .expr6 <- .expr2 * .expr5
    .expr8 <- 1 - .expr4 + .expr2
    .expr10 <- alpha * .expr8^2
    .expr13 <- .expr2 * log(E)
    .value <- .expr6/.expr10
    .grad <- array(0, c(length(.value), 1), list(NULL, c("s")))
    .grad[, "s"] <- .expr13 * .expr5/.expr10 - .expr6 * (alpha * 
        (2 * (.expr13 * .expr8)))/.expr10^2
    attr(.value, "gradient") <- .grad
    .value
}, {
    .expr2 <- E^(alpha + s)
    .expr4 <- E^alpha
    .expr5 <- -1 + .expr4
    .expr6 <- .expr2 * .expr5
    .expr7 <- .expr5 + .expr2
    .expr8 <- .expr6 * .expr7
    .expr10 <- 1 - .expr4 + .expr2
    .expr12 <- alpha * .expr10^3
    .expr16 <- .expr2 * log(E)
    .value <- -(.expr8/.expr12)
    .grad <- array(0, c(length(.value), 1), list(NULL, c("s")))
    .grad[, "s"] <- -((.expr16 * .expr5 * .expr7 + .expr6 * .expr16)/.expr12 - 
        .expr8 * (alpha * (3 * (.expr16 * .expr10^2)))/.expr12^2)
    attr(.value, "gradient") <- .grad
    .value
}, {
    .expr1 <- alpha + s
    .expr2 <- E^.expr1
    .expr4 <- E^alpha
    .expr5 <- -1 + .expr4
    .expr6 <- .expr2 * .expr5
    .expr9 <- 2 * alpha
    .expr15 <- E^(2 * .expr1)
    .expr18 <- E^(.expr9 + s)
    .expr20 <- 1 - 2 * .expr4 + E^.expr9 - 4 * .expr2 + .expr15 + 
        4 * .expr18
    .expr21 <- .expr6 * .expr20
    .expr23 <- 1 - .expr4 + .expr2
    .expr25 <- alpha * .expr23^4
    .expr27 <- log(E)
    .expr28 <- .expr2 * .expr27
    .value <- .expr21/.expr25
    .grad <- array(0, c(length(.value), 1), list(NULL, c("s")))
    .grad[, "s"] <- (.expr28 * .expr5 * .expr20 + .expr6 * (.expr15 * 
        (.expr27 * 2) - 4 * .expr28 + 4 * (.expr18 * .expr27)))/.expr25 - 
        .expr21 * (alpha * (4 * (.expr28 * .expr23^3)))/.expr25^2
    attr(.value, "gradient") <- .grad
    .value
}, {
    .expr1 <- alpha + s
    .expr2 <- E^.expr1
    .expr3 <- -1
    .expr4 <- E^alpha
    .expr5 <- .expr3 + .expr4
    .expr6 <- .expr2 * .expr5
    .expr9 <- 2 * alpha
    .expr13 <- 3 * alpha
    .expr19 <- E^(2 * .expr1)
    .expr23 <- E^(3 * .expr1)
    .expr26 <- E^(.expr9 + s)
    .expr30 <- E^(.expr13 + s)
    .expr35 <- E^(.expr13 + 2 * s)
    .expr37 <- .expr3 + 3 * .expr4 - 3 * E^.expr9 + E^.expr13 + 
        11 * .expr2 - 11 * .expr19 + .expr23 - 22 * .expr26 + 
        11 * .expr30 + 11 * .expr35
    .expr38 <- .expr6 * .expr37
    .expr40 <- 1 - .expr4 + .expr2
    .expr42 <- alpha * .expr40^5
    .expr45 <- log(E)
    .expr46 <- .expr2 * .expr45
    .expr50 <- .expr45 * 2
    .value <- -(.expr38/.expr42)
    .grad <- array(0, c(length(.value), 1), list(NULL, c("s")))
    .grad[, "s"] <- -((.expr46 * .expr5 * .expr37 + .expr6 * 
        (11 * .expr46 - 11 * (.expr19 * .expr50) + .expr23 * 
            (.expr45 * 3) - 22 * (.expr26 * .expr45) + 11 * (.expr30 * 
            .expr45) + 11 * (.expr35 * .expr50)))/.expr42 - .expr38 * 
        (alpha * (5 * (.expr46 * .expr40^4)))/.expr42^2)
    attr(.value, "gradient") <- .grad
    .value
})
