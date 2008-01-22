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


`gumbelCopula.pdf.expr` <-
expression(-(((-log(u1))^alpha)^(1/alpha)/(E^((-log(u1))^alpha)^(1/alpha) * 
    u1 * log(u1))), ((-log(u1))^(-1 + alpha) * (-1 + alpha + 
    ((-log(u1))^alpha + (-log(u2))^alpha)^(1/alpha)) * ((-log(u1))^alpha + 
    (-log(u2))^alpha)^(-2 + 1/alpha) * (-log(u2))^(-1 + alpha))/(E^((-log(u1))^alpha + 
    (-log(u2))^alpha)^(1/alpha) * u1 * u2), ((-log(u1))^(-1 + 
    alpha) * (-log(u2))^(-1 + alpha) * (1 + 2 * alpha^2 + 3 * 
    alpha * (-1 + ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha)^(1/alpha)) - 
    3 * ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha)^(1/alpha) + 
    ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha)^(2/alpha)) * 
    ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha)^(-3 + 
        1/alpha) * (-log(u3))^(-1 + alpha))/(E^((-log(u1))^alpha + 
    (-log(u2))^alpha + (-log(u3))^alpha)^(1/alpha) * u1 * u2 * 
    u3), ((-log(u1))^(-1 + alpha) * (-log(u2))^(-1 + alpha) * 
    (-log(u3))^(-1 + alpha) * (-1 + 6 * alpha^3 + 11 * alpha^2 * 
    (-1 + ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
        (-log(u4))^alpha)^(1/alpha)) + 6 * alpha * (1 - 3 * ((-log(u1))^alpha + 
    (-log(u2))^alpha + (-log(u3))^alpha + (-log(u4))^alpha)^(1/alpha) + 
    ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
        (-log(u4))^alpha)^(2/alpha)) + 7 * ((-log(u1))^alpha + 
    (-log(u2))^alpha + (-log(u3))^alpha + (-log(u4))^alpha)^(1/alpha) - 
    6 * ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
        (-log(u4))^alpha)^(2/alpha) + ((-log(u1))^alpha + (-log(u2))^alpha + 
    (-log(u3))^alpha + (-log(u4))^alpha)^(3/alpha)) * ((-log(u1))^alpha + 
    (-log(u2))^alpha + (-log(u3))^alpha + (-log(u4))^alpha)^(-4 + 
    1/alpha) * (-log(u4))^(-1 + alpha))/(E^((-log(u1))^alpha + 
    (-log(u2))^alpha + (-log(u3))^alpha + (-log(u4))^alpha)^(1/alpha) * 
    u1 * u2 * u3 * u4), -(((-log(u1))^(-1 + alpha) * (-log(u2))^(-1 + 
    alpha) * (-log(u3))^(-1 + alpha) * (-log(u4))^(-1 + alpha) * 
    (-1 - 24 * alpha^4 - 50 * alpha^3 * (-1 + ((-log(u1))^alpha + 
        (-log(u2))^alpha + (-log(u3))^alpha + (-log(u4))^alpha + 
        (-log(u5))^alpha)^(1/alpha)) - 35 * alpha^2 * (1 - 3 * 
        ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
            (-log(u4))^alpha + (-log(u5))^alpha)^(1/alpha) + 
        ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
            (-log(u4))^alpha + (-log(u5))^alpha)^(2/alpha)) - 
        10 * alpha * (-1 + 7 * ((-log(u1))^alpha + (-log(u2))^alpha + 
            (-log(u3))^alpha + (-log(u4))^alpha + (-log(u5))^alpha)^(1/alpha) - 
            6 * ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
                (-log(u4))^alpha + (-log(u5))^alpha)^(2/alpha) + 
            ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
                (-log(u4))^alpha + (-log(u5))^alpha)^(3/alpha)) + 
        15 * ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
            (-log(u4))^alpha + (-log(u5))^alpha)^(1/alpha) - 
        25 * ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
            (-log(u4))^alpha + (-log(u5))^alpha)^(2/alpha) + 
        10 * ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
            (-log(u4))^alpha + (-log(u5))^alpha)^(3/alpha) - 
        ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
            (-log(u4))^alpha + (-log(u5))^alpha)^(4/alpha)) * 
    ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
        (-log(u4))^alpha + (-log(u5))^alpha)^(-5 + 1/alpha) * 
    (-log(u5))^(-1 + alpha))/(E^((-log(u1))^alpha + (-log(u2))^alpha + 
    (-log(u3))^alpha + (-log(u4))^alpha + (-log(u5))^alpha)^(1/alpha) * 
    u1 * u2 * u3 * u4 * u5)), ((-log(u1))^(-1 + alpha) * (-log(u2))^(-1 + 
    alpha) * (-log(u3))^(-1 + alpha) * (-log(u4))^(-1 + alpha) * 
    (-log(u5))^(-1 + alpha) * (-1 + 120 * alpha^5 + 274 * alpha^4 * 
    (-1 + ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
        (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha)^(1/alpha)) + 
    225 * alpha^3 * (1 - 3 * ((-log(u1))^alpha + (-log(u2))^alpha + 
        (-log(u3))^alpha + (-log(u4))^alpha + (-log(u5))^alpha + 
        (-log(u6))^alpha)^(1/alpha) + ((-log(u1))^alpha + (-log(u2))^alpha + 
        (-log(u3))^alpha + (-log(u4))^alpha + (-log(u5))^alpha + 
        (-log(u6))^alpha)^(2/alpha)) + 85 * alpha^2 * (-1 + 7 * 
    ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
        (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha)^(1/alpha) - 
    6 * ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
        (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha)^(2/alpha) + 
    ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
        (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha)^(3/alpha)) + 
    15 * alpha * (1 - 15 * ((-log(u1))^alpha + (-log(u2))^alpha + 
        (-log(u3))^alpha + (-log(u4))^alpha + (-log(u5))^alpha + 
        (-log(u6))^alpha)^(1/alpha) + 25 * ((-log(u1))^alpha + 
        (-log(u2))^alpha + (-log(u3))^alpha + (-log(u4))^alpha + 
        (-log(u5))^alpha + (-log(u6))^alpha)^(2/alpha) - 10 * 
        ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
            (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha)^(3/alpha) + 
        ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
            (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha)^(4/alpha)) + 
    31 * ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
        (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha)^(1/alpha) - 
    90 * ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
        (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha)^(2/alpha) + 
    65 * ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
        (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha)^(3/alpha) - 
    15 * ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
        (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha)^(4/alpha) + 
    ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
        (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha)^(5/alpha)) * 
    ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
        (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha)^(-6 + 
        1/alpha) * (-log(u6))^(-1 + alpha))/(E^((-log(u1))^alpha + 
    (-log(u2))^alpha + (-log(u3))^alpha + (-log(u4))^alpha + 
    (-log(u5))^alpha + (-log(u6))^alpha)^(1/alpha) * u1 * u2 * 
    u3 * u4 * u5 * u6), -(((-log(u1))^(-1 + alpha) * (-log(u2))^(-1 + 
    alpha) * (-log(u3))^(-1 + alpha) * (-log(u4))^(-1 + alpha) * 
    (-log(u5))^(-1 + alpha) * (-log(u6))^(-1 + alpha) * (-1 - 
    720 * alpha^6 - 1764 * alpha^5 * (-1 + ((-log(u1))^alpha + 
    (-log(u2))^alpha + (-log(u3))^alpha + (-log(u4))^alpha + 
    (-log(u5))^alpha + (-log(u6))^alpha + (-log(u7))^alpha)^(1/alpha)) - 
    1624 * alpha^4 * (1 - 3 * ((-log(u1))^alpha + (-log(u2))^alpha + 
        (-log(u3))^alpha + (-log(u4))^alpha + (-log(u5))^alpha + 
        (-log(u6))^alpha + (-log(u7))^alpha)^(1/alpha) + ((-log(u1))^alpha + 
        (-log(u2))^alpha + (-log(u3))^alpha + (-log(u4))^alpha + 
        (-log(u5))^alpha + (-log(u6))^alpha + (-log(u7))^alpha)^(2/alpha)) - 
    735 * alpha^3 * (-1 + 7 * ((-log(u1))^alpha + (-log(u2))^alpha + 
        (-log(u3))^alpha + (-log(u4))^alpha + (-log(u5))^alpha + 
        (-log(u6))^alpha + (-log(u7))^alpha)^(1/alpha) - 6 * 
        ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
            (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
            (-log(u7))^alpha)^(2/alpha) + ((-log(u1))^alpha + 
        (-log(u2))^alpha + (-log(u3))^alpha + (-log(u4))^alpha + 
        (-log(u5))^alpha + (-log(u6))^alpha + (-log(u7))^alpha)^(3/alpha)) - 
    175 * alpha^2 * (1 - 15 * ((-log(u1))^alpha + (-log(u2))^alpha + 
        (-log(u3))^alpha + (-log(u4))^alpha + (-log(u5))^alpha + 
        (-log(u6))^alpha + (-log(u7))^alpha)^(1/alpha) + 25 * 
        ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
            (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
            (-log(u7))^alpha)^(2/alpha) - 10 * ((-log(u1))^alpha + 
        (-log(u2))^alpha + (-log(u3))^alpha + (-log(u4))^alpha + 
        (-log(u5))^alpha + (-log(u6))^alpha + (-log(u7))^alpha)^(3/alpha) + 
        ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
            (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
            (-log(u7))^alpha)^(4/alpha)) - 21 * alpha * (-1 + 
    31 * ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
        (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
        (-log(u7))^alpha)^(1/alpha) - 90 * ((-log(u1))^alpha + 
    (-log(u2))^alpha + (-log(u3))^alpha + (-log(u4))^alpha + 
    (-log(u5))^alpha + (-log(u6))^alpha + (-log(u7))^alpha)^(2/alpha) + 
    65 * ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
        (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
        (-log(u7))^alpha)^(3/alpha) - 15 * ((-log(u1))^alpha + 
    (-log(u2))^alpha + (-log(u3))^alpha + (-log(u4))^alpha + 
    (-log(u5))^alpha + (-log(u6))^alpha + (-log(u7))^alpha)^(4/alpha) + 
    ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
        (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
        (-log(u7))^alpha)^(5/alpha)) + 63 * ((-log(u1))^alpha + 
    (-log(u2))^alpha + (-log(u3))^alpha + (-log(u4))^alpha + 
    (-log(u5))^alpha + (-log(u6))^alpha + (-log(u7))^alpha)^(1/alpha) - 
    301 * ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
        (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
        (-log(u7))^alpha)^(2/alpha) + 350 * ((-log(u1))^alpha + 
    (-log(u2))^alpha + (-log(u3))^alpha + (-log(u4))^alpha + 
    (-log(u5))^alpha + (-log(u6))^alpha + (-log(u7))^alpha)^(3/alpha) - 
    140 * ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
        (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
        (-log(u7))^alpha)^(4/alpha) + 21 * ((-log(u1))^alpha + 
    (-log(u2))^alpha + (-log(u3))^alpha + (-log(u4))^alpha + 
    (-log(u5))^alpha + (-log(u6))^alpha + (-log(u7))^alpha)^(5/alpha) - 
    ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
        (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
        (-log(u7))^alpha)^(6/alpha)) * ((-log(u1))^alpha + (-log(u2))^alpha + 
    (-log(u3))^alpha + (-log(u4))^alpha + (-log(u5))^alpha + 
    (-log(u6))^alpha + (-log(u7))^alpha)^(-7 + 1/alpha) * (-log(u7))^(-1 + 
    alpha))/(E^((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
    (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
    (-log(u7))^alpha)^(1/alpha) * u1 * u2 * u3 * u4 * u5 * u6 * 
    u7)), ((-log(u1))^(-1 + alpha) * (-log(u2))^(-1 + alpha) * 
    (-log(u3))^(-1 + alpha) * (-log(u4))^(-1 + alpha) * (-log(u5))^(-1 + 
    alpha) * (-log(u6))^(-1 + alpha) * (-log(u7))^(-1 + alpha) * 
    (-1 + 5040 * alpha^7 + 13068 * alpha^6 * (-1 + ((-log(u1))^alpha + 
        (-log(u2))^alpha + (-log(u3))^alpha + (-log(u4))^alpha + 
        (-log(u5))^alpha + (-log(u6))^alpha + (-log(u7))^alpha + 
        (-log(u8))^alpha)^(1/alpha)) + 13132 * alpha^5 * (1 - 
        3 * ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
            (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
            (-log(u7))^alpha + (-log(u8))^alpha)^(1/alpha) + 
        ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
            (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
            (-log(u7))^alpha + (-log(u8))^alpha)^(2/alpha)) + 
        6769 * alpha^4 * (-1 + 7 * ((-log(u1))^alpha + (-log(u2))^alpha + 
            (-log(u3))^alpha + (-log(u4))^alpha + (-log(u5))^alpha + 
            (-log(u6))^alpha + (-log(u7))^alpha + (-log(u8))^alpha)^(1/alpha) - 
            6 * ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
                (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
                (-log(u7))^alpha + (-log(u8))^alpha)^(2/alpha) + 
            ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
                (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
                (-log(u7))^alpha + (-log(u8))^alpha)^(3/alpha)) + 
        1960 * alpha^3 * (1 - 15 * ((-log(u1))^alpha + (-log(u2))^alpha + 
            (-log(u3))^alpha + (-log(u4))^alpha + (-log(u5))^alpha + 
            (-log(u6))^alpha + (-log(u7))^alpha + (-log(u8))^alpha)^(1/alpha) + 
            25 * ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
                (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
                (-log(u7))^alpha + (-log(u8))^alpha)^(2/alpha) - 
            10 * ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
                (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
                (-log(u7))^alpha + (-log(u8))^alpha)^(3/alpha) + 
            ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
                (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
                (-log(u7))^alpha + (-log(u8))^alpha)^(4/alpha)) + 
        322 * alpha^2 * (-1 + 31 * ((-log(u1))^alpha + (-log(u2))^alpha + 
            (-log(u3))^alpha + (-log(u4))^alpha + (-log(u5))^alpha + 
            (-log(u6))^alpha + (-log(u7))^alpha + (-log(u8))^alpha)^(1/alpha) - 
            90 * ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
                (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
                (-log(u7))^alpha + (-log(u8))^alpha)^(2/alpha) + 
            65 * ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
                (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
                (-log(u7))^alpha + (-log(u8))^alpha)^(3/alpha) - 
            15 * ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
                (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
                (-log(u7))^alpha + (-log(u8))^alpha)^(4/alpha) + 
            ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
                (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
                (-log(u7))^alpha + (-log(u8))^alpha)^(5/alpha)) + 
        28 * alpha * (1 - 63 * ((-log(u1))^alpha + (-log(u2))^alpha + 
            (-log(u3))^alpha + (-log(u4))^alpha + (-log(u5))^alpha + 
            (-log(u6))^alpha + (-log(u7))^alpha + (-log(u8))^alpha)^(1/alpha) + 
            301 * ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
                (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
                (-log(u7))^alpha + (-log(u8))^alpha)^(2/alpha) - 
            350 * ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
                (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
                (-log(u7))^alpha + (-log(u8))^alpha)^(3/alpha) + 
            140 * ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
                (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
                (-log(u7))^alpha + (-log(u8))^alpha)^(4/alpha) - 
            21 * ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
                (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
                (-log(u7))^alpha + (-log(u8))^alpha)^(5/alpha) + 
            ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
                (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
                (-log(u7))^alpha + (-log(u8))^alpha)^(6/alpha)) + 
        127 * ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
            (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
            (-log(u7))^alpha + (-log(u8))^alpha)^(1/alpha) - 
        966 * ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
            (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
            (-log(u7))^alpha + (-log(u8))^alpha)^(2/alpha) + 
        1701 * ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
            (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
            (-log(u7))^alpha + (-log(u8))^alpha)^(3/alpha) - 
        1050 * ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
            (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
            (-log(u7))^alpha + (-log(u8))^alpha)^(4/alpha) + 
        266 * ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
            (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
            (-log(u7))^alpha + (-log(u8))^alpha)^(5/alpha) - 
        28 * ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
            (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
            (-log(u7))^alpha + (-log(u8))^alpha)^(6/alpha) + 
        ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
            (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
            (-log(u7))^alpha + (-log(u8))^alpha)^(7/alpha)) * 
    ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
        (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
        (-log(u7))^alpha + (-log(u8))^alpha)^(-8 + 1/alpha) * 
    (-log(u8))^(-1 + alpha))/(E^((-log(u1))^alpha + (-log(u2))^alpha + 
    (-log(u3))^alpha + (-log(u4))^alpha + (-log(u5))^alpha + 
    (-log(u6))^alpha + (-log(u7))^alpha + (-log(u8))^alpha)^(1/alpha) * 
    u1 * u2 * u3 * u4 * u5 * u6 * u7 * u8), ((-log(u1))^(-1 + 
    alpha) * (-log(u2))^(-1 + alpha) * (-log(u3))^(-1 + alpha) * 
    (-log(u4))^(-1 + alpha) * (-log(u5))^(-1 + alpha) * (-log(u6))^(-1 + 
    alpha) * (-log(u7))^(-1 + alpha) * (-log(u8))^(-1 + alpha) * 
    (1 + 40320 * alpha^8 + 109584 * alpha^7 * (-1 + ((-log(u1))^alpha + 
        (-log(u2))^alpha + (-log(u3))^alpha + (-log(u4))^alpha + 
        (-log(u5))^alpha + (-log(u6))^alpha + (-log(u7))^alpha + 
        (-log(u8))^alpha + (-log(u9))^alpha)^(1/alpha)) + 118124 * 
        alpha^6 * (1 - 3 * ((-log(u1))^alpha + (-log(u2))^alpha + 
        (-log(u3))^alpha + (-log(u4))^alpha + (-log(u5))^alpha + 
        (-log(u6))^alpha + (-log(u7))^alpha + (-log(u8))^alpha + 
        (-log(u9))^alpha)^(1/alpha) + ((-log(u1))^alpha + (-log(u2))^alpha + 
        (-log(u3))^alpha + (-log(u4))^alpha + (-log(u5))^alpha + 
        (-log(u6))^alpha + (-log(u7))^alpha + (-log(u8))^alpha + 
        (-log(u9))^alpha)^(2/alpha)) + 67284 * alpha^5 * (-1 + 
        7 * ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
            (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
            (-log(u7))^alpha + (-log(u8))^alpha + (-log(u9))^alpha)^(1/alpha) - 
        6 * ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
            (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
            (-log(u7))^alpha + (-log(u8))^alpha + (-log(u9))^alpha)^(2/alpha) + 
        ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
            (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
            (-log(u7))^alpha + (-log(u8))^alpha + (-log(u9))^alpha)^(3/alpha)) + 
        22449 * alpha^4 * (1 - 15 * ((-log(u1))^alpha + (-log(u2))^alpha + 
            (-log(u3))^alpha + (-log(u4))^alpha + (-log(u5))^alpha + 
            (-log(u6))^alpha + (-log(u7))^alpha + (-log(u8))^alpha + 
            (-log(u9))^alpha)^(1/alpha) + 25 * ((-log(u1))^alpha + 
            (-log(u2))^alpha + (-log(u3))^alpha + (-log(u4))^alpha + 
            (-log(u5))^alpha + (-log(u6))^alpha + (-log(u7))^alpha + 
            (-log(u8))^alpha + (-log(u9))^alpha)^(2/alpha) - 
            10 * ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
                (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
                (-log(u7))^alpha + (-log(u8))^alpha + (-log(u9))^alpha)^(3/alpha) + 
            ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
                (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
                (-log(u7))^alpha + (-log(u8))^alpha + (-log(u9))^alpha)^(4/alpha)) + 
        4536 * alpha^3 * (-1 + 31 * ((-log(u1))^alpha + (-log(u2))^alpha + 
            (-log(u3))^alpha + (-log(u4))^alpha + (-log(u5))^alpha + 
            (-log(u6))^alpha + (-log(u7))^alpha + (-log(u8))^alpha + 
            (-log(u9))^alpha)^(1/alpha) - 90 * ((-log(u1))^alpha + 
            (-log(u2))^alpha + (-log(u3))^alpha + (-log(u4))^alpha + 
            (-log(u5))^alpha + (-log(u6))^alpha + (-log(u7))^alpha + 
            (-log(u8))^alpha + (-log(u9))^alpha)^(2/alpha) + 
            65 * ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
                (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
                (-log(u7))^alpha + (-log(u8))^alpha + (-log(u9))^alpha)^(3/alpha) - 
            15 * ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
                (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
                (-log(u7))^alpha + (-log(u8))^alpha + (-log(u9))^alpha)^(4/alpha) + 
            ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
                (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
                (-log(u7))^alpha + (-log(u8))^alpha + (-log(u9))^alpha)^(5/alpha)) + 
        546 * alpha^2 * (1 - 63 * ((-log(u1))^alpha + (-log(u2))^alpha + 
            (-log(u3))^alpha + (-log(u4))^alpha + (-log(u5))^alpha + 
            (-log(u6))^alpha + (-log(u7))^alpha + (-log(u8))^alpha + 
            (-log(u9))^alpha)^(1/alpha) + 301 * ((-log(u1))^alpha + 
            (-log(u2))^alpha + (-log(u3))^alpha + (-log(u4))^alpha + 
            (-log(u5))^alpha + (-log(u6))^alpha + (-log(u7))^alpha + 
            (-log(u8))^alpha + (-log(u9))^alpha)^(2/alpha) - 
            350 * ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
                (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
                (-log(u7))^alpha + (-log(u8))^alpha + (-log(u9))^alpha)^(3/alpha) + 
            140 * ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
                (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
                (-log(u7))^alpha + (-log(u8))^alpha + (-log(u9))^alpha)^(4/alpha) - 
            21 * ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
                (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
                (-log(u7))^alpha + (-log(u8))^alpha + (-log(u9))^alpha)^(5/alpha) + 
            ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
                (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
                (-log(u7))^alpha + (-log(u8))^alpha + (-log(u9))^alpha)^(6/alpha)) + 
        36 * alpha * (-1 + 127 * ((-log(u1))^alpha + (-log(u2))^alpha + 
            (-log(u3))^alpha + (-log(u4))^alpha + (-log(u5))^alpha + 
            (-log(u6))^alpha + (-log(u7))^alpha + (-log(u8))^alpha + 
            (-log(u9))^alpha)^(1/alpha) - 966 * ((-log(u1))^alpha + 
            (-log(u2))^alpha + (-log(u3))^alpha + (-log(u4))^alpha + 
            (-log(u5))^alpha + (-log(u6))^alpha + (-log(u7))^alpha + 
            (-log(u8))^alpha + (-log(u9))^alpha)^(2/alpha) + 
            1701 * ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
                (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
                (-log(u7))^alpha + (-log(u8))^alpha + (-log(u9))^alpha)^(3/alpha) - 
            1050 * ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
                (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
                (-log(u7))^alpha + (-log(u8))^alpha + (-log(u9))^alpha)^(4/alpha) + 
            266 * ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
                (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
                (-log(u7))^alpha + (-log(u8))^alpha + (-log(u9))^alpha)^(5/alpha) - 
            28 * ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
                (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
                (-log(u7))^alpha + (-log(u8))^alpha + (-log(u9))^alpha)^(6/alpha) + 
            ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
                (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
                (-log(u7))^alpha + (-log(u8))^alpha + (-log(u9))^alpha)^(7/alpha)) - 
        255 * ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
            (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
            (-log(u7))^alpha + (-log(u8))^alpha + (-log(u9))^alpha)^(1/alpha) + 
        3025 * ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
            (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
            (-log(u7))^alpha + (-log(u8))^alpha + (-log(u9))^alpha)^(2/alpha) - 
        7770 * ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
            (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
            (-log(u7))^alpha + (-log(u8))^alpha + (-log(u9))^alpha)^(3/alpha) + 
        6951 * ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
            (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
            (-log(u7))^alpha + (-log(u8))^alpha + (-log(u9))^alpha)^(4/alpha) - 
        2646 * ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
            (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
            (-log(u7))^alpha + (-log(u8))^alpha + (-log(u9))^alpha)^(5/alpha) + 
        462 * ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
            (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
            (-log(u7))^alpha + (-log(u8))^alpha + (-log(u9))^alpha)^(6/alpha) - 
        36 * ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
            (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
            (-log(u7))^alpha + (-log(u8))^alpha + (-log(u9))^alpha)^(7/alpha) + 
        ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
            (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
            (-log(u7))^alpha + (-log(u8))^alpha + (-log(u9))^alpha)^(8/alpha)) * 
    ((-log(u1))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
        (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
        (-log(u7))^alpha + (-log(u8))^alpha + (-log(u9))^alpha)^(-9 + 
        1/alpha) * (-log(u9))^(-1 + alpha))/(E^((-log(u1))^alpha + 
    (-log(u2))^alpha + (-log(u3))^alpha + (-log(u4))^alpha + 
    (-log(u5))^alpha + (-log(u6))^alpha + (-log(u7))^alpha + 
    (-log(u8))^alpha + (-log(u9))^alpha)^(1/alpha) * u1 * u2 * 
    u3 * u4 * u5 * u6 * u7 * u8 * u9), ((-log(u1))^(-1 + alpha) * 
    (-log(u10))^(-1 + alpha) * (-log(u2))^(-1 + alpha) * (-log(u3))^(-1 + 
    alpha) * (-log(u4))^(-1 + alpha) * (-log(u5))^(-1 + alpha) * 
    (-log(u6))^(-1 + alpha) * (-log(u7))^(-1 + alpha) * (-log(u8))^(-1 + 
    alpha) * (-1 + 362880 * alpha^9 + 1026576 * alpha^8 * (-1 + 
    ((-log(u1))^alpha + (-log(u10))^alpha + (-log(u2))^alpha + 
        (-log(u3))^alpha + (-log(u4))^alpha + (-log(u5))^alpha + 
        (-log(u6))^alpha + (-log(u7))^alpha + (-log(u8))^alpha + 
        (-log(u9))^alpha)^(1/alpha)) + 1172700 * alpha^7 * (1 - 
    3 * ((-log(u1))^alpha + (-log(u10))^alpha + (-log(u2))^alpha + 
        (-log(u3))^alpha + (-log(u4))^alpha + (-log(u5))^alpha + 
        (-log(u6))^alpha + (-log(u7))^alpha + (-log(u8))^alpha + 
        (-log(u9))^alpha)^(1/alpha) + ((-log(u1))^alpha + (-log(u10))^alpha + 
    (-log(u2))^alpha + (-log(u3))^alpha + (-log(u4))^alpha + 
    (-log(u5))^alpha + (-log(u6))^alpha + (-log(u7))^alpha + 
    (-log(u8))^alpha + (-log(u9))^alpha)^(2/alpha)) + 723680 * 
    alpha^6 * (-1 + 7 * ((-log(u1))^alpha + (-log(u10))^alpha + 
    (-log(u2))^alpha + (-log(u3))^alpha + (-log(u4))^alpha + 
    (-log(u5))^alpha + (-log(u6))^alpha + (-log(u7))^alpha + 
    (-log(u8))^alpha + (-log(u9))^alpha)^(1/alpha) - 6 * ((-log(u1))^alpha + 
    (-log(u10))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
    (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
    (-log(u7))^alpha + (-log(u8))^alpha + (-log(u9))^alpha)^(2/alpha) + 
    ((-log(u1))^alpha + (-log(u10))^alpha + (-log(u2))^alpha + 
        (-log(u3))^alpha + (-log(u4))^alpha + (-log(u5))^alpha + 
        (-log(u6))^alpha + (-log(u7))^alpha + (-log(u8))^alpha + 
        (-log(u9))^alpha)^(3/alpha)) + 269325 * alpha^5 * (1 - 
    15 * ((-log(u1))^alpha + (-log(u10))^alpha + (-log(u2))^alpha + 
        (-log(u3))^alpha + (-log(u4))^alpha + (-log(u5))^alpha + 
        (-log(u6))^alpha + (-log(u7))^alpha + (-log(u8))^alpha + 
        (-log(u9))^alpha)^(1/alpha) + 25 * ((-log(u1))^alpha + 
    (-log(u10))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
    (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
    (-log(u7))^alpha + (-log(u8))^alpha + (-log(u9))^alpha)^(2/alpha) - 
    10 * ((-log(u1))^alpha + (-log(u10))^alpha + (-log(u2))^alpha + 
        (-log(u3))^alpha + (-log(u4))^alpha + (-log(u5))^alpha + 
        (-log(u6))^alpha + (-log(u7))^alpha + (-log(u8))^alpha + 
        (-log(u9))^alpha)^(3/alpha) + ((-log(u1))^alpha + (-log(u10))^alpha + 
    (-log(u2))^alpha + (-log(u3))^alpha + (-log(u4))^alpha + 
    (-log(u5))^alpha + (-log(u6))^alpha + (-log(u7))^alpha + 
    (-log(u8))^alpha + (-log(u9))^alpha)^(4/alpha)) + 63273 * 
    alpha^4 * (-1 + 31 * ((-log(u1))^alpha + (-log(u10))^alpha + 
    (-log(u2))^alpha + (-log(u3))^alpha + (-log(u4))^alpha + 
    (-log(u5))^alpha + (-log(u6))^alpha + (-log(u7))^alpha + 
    (-log(u8))^alpha + (-log(u9))^alpha)^(1/alpha) - 90 * ((-log(u1))^alpha + 
    (-log(u10))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
    (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
    (-log(u7))^alpha + (-log(u8))^alpha + (-log(u9))^alpha)^(2/alpha) + 
    65 * ((-log(u1))^alpha + (-log(u10))^alpha + (-log(u2))^alpha + 
        (-log(u3))^alpha + (-log(u4))^alpha + (-log(u5))^alpha + 
        (-log(u6))^alpha + (-log(u7))^alpha + (-log(u8))^alpha + 
        (-log(u9))^alpha)^(3/alpha) - 15 * ((-log(u1))^alpha + 
    (-log(u10))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
    (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
    (-log(u7))^alpha + (-log(u8))^alpha + (-log(u9))^alpha)^(4/alpha) + 
    ((-log(u1))^alpha + (-log(u10))^alpha + (-log(u2))^alpha + 
        (-log(u3))^alpha + (-log(u4))^alpha + (-log(u5))^alpha + 
        (-log(u6))^alpha + (-log(u7))^alpha + (-log(u8))^alpha + 
        (-log(u9))^alpha)^(5/alpha)) + 9450 * alpha^3 * (1 - 
    63 * ((-log(u1))^alpha + (-log(u10))^alpha + (-log(u2))^alpha + 
        (-log(u3))^alpha + (-log(u4))^alpha + (-log(u5))^alpha + 
        (-log(u6))^alpha + (-log(u7))^alpha + (-log(u8))^alpha + 
        (-log(u9))^alpha)^(1/alpha) + 301 * ((-log(u1))^alpha + 
    (-log(u10))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
    (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
    (-log(u7))^alpha + (-log(u8))^alpha + (-log(u9))^alpha)^(2/alpha) - 
    350 * ((-log(u1))^alpha + (-log(u10))^alpha + (-log(u2))^alpha + 
        (-log(u3))^alpha + (-log(u4))^alpha + (-log(u5))^alpha + 
        (-log(u6))^alpha + (-log(u7))^alpha + (-log(u8))^alpha + 
        (-log(u9))^alpha)^(3/alpha) + 140 * ((-log(u1))^alpha + 
    (-log(u10))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
    (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
    (-log(u7))^alpha + (-log(u8))^alpha + (-log(u9))^alpha)^(4/alpha) - 
    21 * ((-log(u1))^alpha + (-log(u10))^alpha + (-log(u2))^alpha + 
        (-log(u3))^alpha + (-log(u4))^alpha + (-log(u5))^alpha + 
        (-log(u6))^alpha + (-log(u7))^alpha + (-log(u8))^alpha + 
        (-log(u9))^alpha)^(5/alpha) + ((-log(u1))^alpha + (-log(u10))^alpha + 
    (-log(u2))^alpha + (-log(u3))^alpha + (-log(u4))^alpha + 
    (-log(u5))^alpha + (-log(u6))^alpha + (-log(u7))^alpha + 
    (-log(u8))^alpha + (-log(u9))^alpha)^(6/alpha)) + 870 * alpha^2 * 
    (-1 + 127 * ((-log(u1))^alpha + (-log(u10))^alpha + (-log(u2))^alpha + 
        (-log(u3))^alpha + (-log(u4))^alpha + (-log(u5))^alpha + 
        (-log(u6))^alpha + (-log(u7))^alpha + (-log(u8))^alpha + 
        (-log(u9))^alpha)^(1/alpha) - 966 * ((-log(u1))^alpha + 
        (-log(u10))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
        (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
        (-log(u7))^alpha + (-log(u8))^alpha + (-log(u9))^alpha)^(2/alpha) + 
        1701 * ((-log(u1))^alpha + (-log(u10))^alpha + (-log(u2))^alpha + 
            (-log(u3))^alpha + (-log(u4))^alpha + (-log(u5))^alpha + 
            (-log(u6))^alpha + (-log(u7))^alpha + (-log(u8))^alpha + 
            (-log(u9))^alpha)^(3/alpha) - 1050 * ((-log(u1))^alpha + 
        (-log(u10))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
        (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
        (-log(u7))^alpha + (-log(u8))^alpha + (-log(u9))^alpha)^(4/alpha) + 
        266 * ((-log(u1))^alpha + (-log(u10))^alpha + (-log(u2))^alpha + 
            (-log(u3))^alpha + (-log(u4))^alpha + (-log(u5))^alpha + 
            (-log(u6))^alpha + (-log(u7))^alpha + (-log(u8))^alpha + 
            (-log(u9))^alpha)^(5/alpha) - 28 * ((-log(u1))^alpha + 
        (-log(u10))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
        (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
        (-log(u7))^alpha + (-log(u8))^alpha + (-log(u9))^alpha)^(6/alpha) + 
        ((-log(u1))^alpha + (-log(u10))^alpha + (-log(u2))^alpha + 
            (-log(u3))^alpha + (-log(u4))^alpha + (-log(u5))^alpha + 
            (-log(u6))^alpha + (-log(u7))^alpha + (-log(u8))^alpha + 
            (-log(u9))^alpha)^(7/alpha)) + 45 * alpha * (1 - 
    255 * ((-log(u1))^alpha + (-log(u10))^alpha + (-log(u2))^alpha + 
        (-log(u3))^alpha + (-log(u4))^alpha + (-log(u5))^alpha + 
        (-log(u6))^alpha + (-log(u7))^alpha + (-log(u8))^alpha + 
        (-log(u9))^alpha)^(1/alpha) + 3025 * ((-log(u1))^alpha + 
    (-log(u10))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
    (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
    (-log(u7))^alpha + (-log(u8))^alpha + (-log(u9))^alpha)^(2/alpha) - 
    7770 * ((-log(u1))^alpha + (-log(u10))^alpha + (-log(u2))^alpha + 
        (-log(u3))^alpha + (-log(u4))^alpha + (-log(u5))^alpha + 
        (-log(u6))^alpha + (-log(u7))^alpha + (-log(u8))^alpha + 
        (-log(u9))^alpha)^(3/alpha) + 6951 * ((-log(u1))^alpha + 
    (-log(u10))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
    (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
    (-log(u7))^alpha + (-log(u8))^alpha + (-log(u9))^alpha)^(4/alpha) - 
    2646 * ((-log(u1))^alpha + (-log(u10))^alpha + (-log(u2))^alpha + 
        (-log(u3))^alpha + (-log(u4))^alpha + (-log(u5))^alpha + 
        (-log(u6))^alpha + (-log(u7))^alpha + (-log(u8))^alpha + 
        (-log(u9))^alpha)^(5/alpha) + 462 * ((-log(u1))^alpha + 
    (-log(u10))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
    (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
    (-log(u7))^alpha + (-log(u8))^alpha + (-log(u9))^alpha)^(6/alpha) - 
    36 * ((-log(u1))^alpha + (-log(u10))^alpha + (-log(u2))^alpha + 
        (-log(u3))^alpha + (-log(u4))^alpha + (-log(u5))^alpha + 
        (-log(u6))^alpha + (-log(u7))^alpha + (-log(u8))^alpha + 
        (-log(u9))^alpha)^(7/alpha) + ((-log(u1))^alpha + (-log(u10))^alpha + 
    (-log(u2))^alpha + (-log(u3))^alpha + (-log(u4))^alpha + 
    (-log(u5))^alpha + (-log(u6))^alpha + (-log(u7))^alpha + 
    (-log(u8))^alpha + (-log(u9))^alpha)^(8/alpha)) + 511 * ((-log(u1))^alpha + 
    (-log(u10))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
    (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
    (-log(u7))^alpha + (-log(u8))^alpha + (-log(u9))^alpha)^(1/alpha) - 
    9330 * ((-log(u1))^alpha + (-log(u10))^alpha + (-log(u2))^alpha + 
        (-log(u3))^alpha + (-log(u4))^alpha + (-log(u5))^alpha + 
        (-log(u6))^alpha + (-log(u7))^alpha + (-log(u8))^alpha + 
        (-log(u9))^alpha)^(2/alpha) + 34105 * ((-log(u1))^alpha + 
    (-log(u10))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
    (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
    (-log(u7))^alpha + (-log(u8))^alpha + (-log(u9))^alpha)^(3/alpha) - 
    42525 * ((-log(u1))^alpha + (-log(u10))^alpha + (-log(u2))^alpha + 
        (-log(u3))^alpha + (-log(u4))^alpha + (-log(u5))^alpha + 
        (-log(u6))^alpha + (-log(u7))^alpha + (-log(u8))^alpha + 
        (-log(u9))^alpha)^(4/alpha) + 22827 * ((-log(u1))^alpha + 
    (-log(u10))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
    (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
    (-log(u7))^alpha + (-log(u8))^alpha + (-log(u9))^alpha)^(5/alpha) - 
    5880 * ((-log(u1))^alpha + (-log(u10))^alpha + (-log(u2))^alpha + 
        (-log(u3))^alpha + (-log(u4))^alpha + (-log(u5))^alpha + 
        (-log(u6))^alpha + (-log(u7))^alpha + (-log(u8))^alpha + 
        (-log(u9))^alpha)^(6/alpha) + 750 * ((-log(u1))^alpha + 
    (-log(u10))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
    (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
    (-log(u7))^alpha + (-log(u8))^alpha + (-log(u9))^alpha)^(7/alpha) - 
    45 * ((-log(u1))^alpha + (-log(u10))^alpha + (-log(u2))^alpha + 
        (-log(u3))^alpha + (-log(u4))^alpha + (-log(u5))^alpha + 
        (-log(u6))^alpha + (-log(u7))^alpha + (-log(u8))^alpha + 
        (-log(u9))^alpha)^(8/alpha) + ((-log(u1))^alpha + (-log(u10))^alpha + 
    (-log(u2))^alpha + (-log(u3))^alpha + (-log(u4))^alpha + 
    (-log(u5))^alpha + (-log(u6))^alpha + (-log(u7))^alpha + 
    (-log(u8))^alpha + (-log(u9))^alpha)^(9/alpha)) * ((-log(u1))^alpha + 
    (-log(u10))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
    (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
    (-log(u7))^alpha + (-log(u8))^alpha + (-log(u9))^alpha)^(-10 + 
    1/alpha) * (-log(u9))^(-1 + alpha))/(E^((-log(u1))^alpha + 
    (-log(u10))^alpha + (-log(u2))^alpha + (-log(u3))^alpha + 
    (-log(u4))^alpha + (-log(u5))^alpha + (-log(u6))^alpha + 
    (-log(u7))^alpha + (-log(u8))^alpha + (-log(u9))^alpha)^(1/alpha) * 
    u1 * u10 * u2 * u3 * u4 * u5 * u6 * u7 * u8 * u9))
`gumbelCopula.pdf.algr` <-
expression({
    .expr1 <- log(u1)
    .expr5 <- ((-.expr1)^alpha)^(1/alpha)
    .value <- -(.expr5/(E^.expr5 * u1 * .expr1))
    .grad <- array(0, c(length(.value), 1), list(NULL, c("s")))
    .grad[, "s"] <- 0
    attr(.value, "gradient") <- .grad
    .value
}, {
    .expr2 <- -log(u1)
    .expr4 <- -1 + alpha
    .expr8 <- -log(u2)
    .expr10 <- .expr2^alpha + .expr8^alpha
    .expr11 <- 1/alpha
    .expr12 <- .expr10^.expr11
    .value <- .expr2^.expr4 * (.expr4 + .expr12) * .expr10^(-2 + 
        .expr11) * .expr8^.expr4/(E^.expr12 * u1 * u2)
    .grad <- array(0, c(length(.value), 1), list(NULL, c("s")))
    .grad[, "s"] <- 0
    attr(.value, "gradient") <- .grad
    .value
}, {
    .expr2 <- -log(u1)
    .expr3 <- -1
    .expr4 <- .expr3 + alpha
    .expr7 <- -log(u2)
    .expr18 <- -log(u3)
    .expr20 <- .expr2^alpha + .expr7^alpha + .expr18^alpha
    .expr21 <- 1/alpha
    .expr22 <- .expr20^.expr21
    .value <- .expr2^.expr4 * .expr7^.expr4 * (1 + 2 * alpha^2 + 
        3 * alpha * (.expr3 + .expr22) - 3 * .expr22 + .expr20^(2/alpha)) * 
        .expr20^(-3 + .expr21) * .expr18^.expr4/(E^.expr22 * 
        u1 * u2 * u3)
    .grad <- array(0, c(length(.value), 1), list(NULL, c("s")))
    .grad[, "s"] <- 0
    attr(.value, "gradient") <- .grad
    .value
}, {
    .expr2 <- -log(u1)
    .expr3 <- -1
    .expr4 <- .expr3 + alpha
    .expr7 <- -log(u2)
    .expr11 <- -log(u3)
    .expr25 <- -log(u4)
    .expr27 <- .expr2^alpha + .expr7^alpha + .expr11^alpha + 
        .expr25^alpha
    .expr28 <- 1/alpha
    .expr29 <- .expr27^.expr28
    .expr37 <- .expr27^(2/alpha)
    .value <- .expr2^.expr4 * .expr7^.expr4 * .expr11^.expr4 * 
        (.expr3 + 6 * alpha^3 + 11 * alpha^2 * (.expr3 + .expr29) + 
            6 * alpha * (1 - 3 * .expr29 + .expr37) + 7 * .expr29 - 
            6 * .expr37 + .expr27^(3/alpha)) * .expr27^(-4 + 
        .expr28) * .expr25^.expr4/(E^.expr29 * u1 * u2 * u3 * 
        u4)
    .grad <- array(0, c(length(.value), 1), list(NULL, c("s")))
    .grad[, "s"] <- 0
    attr(.value, "gradient") <- .grad
    .value
}, {
    .expr2 <- -log(u1)
    .expr3 <- -1
    .expr4 <- .expr3 + alpha
    .expr7 <- -log(u2)
    .expr11 <- -log(u3)
    .expr15 <- -log(u4)
    .expr31 <- -log(u5)
    .expr33 <- .expr2^alpha + .expr7^alpha + .expr11^alpha + 
        .expr15^alpha + .expr31^alpha
    .expr34 <- 1/alpha
    .expr35 <- .expr33^.expr34
    .expr44 <- .expr33^(2/alpha)
    .expr54 <- .expr33^(3/alpha)
    .value <- -(.expr2^.expr4 * .expr7^.expr4 * .expr11^.expr4 * 
        .expr15^.expr4 * (.expr3 - 24 * alpha^4 - 50 * alpha^3 * 
        (.expr3 + .expr35) - 35 * alpha^2 * (1 - 3 * .expr35 + 
        .expr44) - 10 * alpha * (.expr3 + 7 * .expr35 - 6 * .expr44 + 
        .expr54) + 15 * .expr35 - 25 * .expr44 + 10 * .expr54 - 
        .expr33^(4/alpha)) * .expr33^(-5 + .expr34) * .expr31^.expr4/(E^.expr35 * 
        u1 * u2 * u3 * u4 * u5))
    .grad <- array(0, c(length(.value), 1), list(NULL, c("s")))
    .grad[, "s"] <- 0
    attr(.value, "gradient") <- .grad
    .value
}, {
    .expr2 <- -log(u1)
    .expr3 <- -1
    .expr4 <- .expr3 + alpha
    .expr7 <- -log(u2)
    .expr11 <- -log(u3)
    .expr15 <- -log(u4)
    .expr19 <- -log(u5)
    .expr37 <- -log(u6)
    .expr39 <- .expr2^alpha + .expr7^alpha + .expr11^alpha + 
        .expr15^alpha + .expr19^alpha + .expr37^alpha
    .expr40 <- 1/alpha
    .expr41 <- .expr39^.expr40
    .expr50 <- .expr39^(2/alpha)
    .expr61 <- .expr39^(3/alpha)
    .expr73 <- .expr39^(4/alpha)
    .value <- .expr2^.expr4 * .expr7^.expr4 * .expr11^.expr4 * 
        .expr15^.expr4 * .expr19^.expr4 * (.expr3 + 120 * alpha^5 + 
        274 * alpha^4 * (.expr3 + .expr41) + 225 * alpha^3 * 
        (1 - 3 * .expr41 + .expr50) + 85 * alpha^2 * (.expr3 + 
        7 * .expr41 - 6 * .expr50 + .expr61) + 15 * alpha * (1 - 
        15 * .expr41 + 25 * .expr50 - 10 * .expr61 + .expr73) + 
        31 * .expr41 - 90 * .expr50 + 65 * .expr61 - 15 * .expr73 + 
        .expr39^(5/alpha)) * .expr39^(-6 + .expr40) * .expr37^.expr4/(E^.expr41 * 
        u1 * u2 * u3 * u4 * u5 * u6)
    .grad <- array(0, c(length(.value), 1), list(NULL, c("s")))
    .grad[, "s"] <- 0
    attr(.value, "gradient") <- .grad
    .value
}, {
    .expr2 <- -log(u1)
    .expr3 <- -1
    .expr4 <- .expr3 + alpha
    .expr7 <- -log(u2)
    .expr11 <- -log(u3)
    .expr15 <- -log(u4)
    .expr19 <- -log(u5)
    .expr23 <- -log(u6)
    .expr43 <- -log(u7)
    .expr45 <- .expr2^alpha + .expr7^alpha + .expr11^alpha + 
        .expr15^alpha + .expr19^alpha + .expr23^alpha + .expr43^alpha
    .expr46 <- 1/alpha
    .expr47 <- .expr45^.expr46
    .expr56 <- .expr45^(2/alpha)
    .expr67 <- .expr45^(3/alpha)
    .expr80 <- .expr45^(4/alpha)
    .expr94 <- .expr45^(5/alpha)
    .value <- -(.expr2^.expr4 * .expr7^.expr4 * .expr11^.expr4 * 
        .expr15^.expr4 * .expr19^.expr4 * .expr23^.expr4 * (.expr3 - 
        720 * alpha^6 - 1764 * alpha^5 * (.expr3 + .expr47) - 
        1624 * alpha^4 * (1 - 3 * .expr47 + .expr56) - 735 * 
        alpha^3 * (.expr3 + 7 * .expr47 - 6 * .expr56 + .expr67) - 
        175 * alpha^2 * (1 - 15 * .expr47 + 25 * .expr56 - 10 * 
            .expr67 + .expr80) - 21 * alpha * (.expr3 + 31 * 
        .expr47 - 90 * .expr56 + 65 * .expr67 - 15 * .expr80 + 
        .expr94) + 63 * .expr47 - 301 * .expr56 + 350 * .expr67 - 
        140 * .expr80 + 21 * .expr94 - .expr45^(6/alpha)) * .expr45^(-7 + 
        .expr46) * .expr43^.expr4/(E^.expr47 * u1 * u2 * u3 * 
        u4 * u5 * u6 * u7))
    .grad <- array(0, c(length(.value), 1), list(NULL, c("s")))
    .grad[, "s"] <- 0
    attr(.value, "gradient") <- .grad
    .value
}, {
    .expr2 <- -log(u1)
    .expr3 <- -1
    .expr4 <- .expr3 + alpha
    .expr7 <- -log(u2)
    .expr11 <- -log(u3)
    .expr15 <- -log(u4)
    .expr19 <- -log(u5)
    .expr23 <- -log(u6)
    .expr27 <- -log(u7)
    .expr49 <- -log(u8)
    .expr51 <- .expr2^alpha + .expr7^alpha + .expr11^alpha + 
        .expr15^alpha + .expr19^alpha + .expr23^alpha + .expr27^alpha + 
        .expr49^alpha
    .expr52 <- 1/alpha
    .expr53 <- .expr51^.expr52
    .expr62 <- .expr51^(2/alpha)
    .expr73 <- .expr51^(3/alpha)
    .expr86 <- .expr51^(4/alpha)
    .expr101 <- .expr51^(5/alpha)
    .expr117 <- .expr51^(6/alpha)
    .value <- .expr2^.expr4 * .expr7^.expr4 * .expr11^.expr4 * 
        .expr15^.expr4 * .expr19^.expr4 * .expr23^.expr4 * .expr27^.expr4 * 
        (.expr3 + 5040 * alpha^7 + 13068 * alpha^6 * (.expr3 + 
            .expr53) + 13132 * alpha^5 * (1 - 3 * .expr53 + .expr62) + 
            6769 * alpha^4 * (.expr3 + 7 * .expr53 - 6 * .expr62 + 
                .expr73) + 1960 * alpha^3 * (1 - 15 * .expr53 + 
            25 * .expr62 - 10 * .expr73 + .expr86) + 322 * alpha^2 * 
            (.expr3 + 31 * .expr53 - 90 * .expr62 + 65 * .expr73 - 
                15 * .expr86 + .expr101) + 28 * alpha * (1 - 
            63 * .expr53 + 301 * .expr62 - 350 * .expr73 + 140 * 
            .expr86 - 21 * .expr101 + .expr117) + 127 * .expr53 - 
            966 * .expr62 + 1701 * .expr73 - 1050 * .expr86 + 
            266 * .expr101 - 28 * .expr117 + .expr51^(7/alpha)) * 
        .expr51^(-8 + .expr52) * .expr49^.expr4/(E^.expr53 * 
        u1 * u2 * u3 * u4 * u5 * u6 * u7 * u8)
    .grad <- array(0, c(length(.value), 1), list(NULL, c("s")))
    .grad[, "s"] <- 0
    attr(.value, "gradient") <- .grad
    .value
}, {
    .expr2 <- -log(u1)
    .expr3 <- -1
    .expr4 <- .expr3 + alpha
    .expr7 <- -log(u2)
    .expr11 <- -log(u3)
    .expr15 <- -log(u4)
    .expr19 <- -log(u5)
    .expr23 <- -log(u6)
    .expr27 <- -log(u7)
    .expr31 <- -log(u8)
    .expr55 <- -log(u9)
    .expr57 <- .expr2^alpha + .expr7^alpha + .expr11^alpha + 
        .expr15^alpha + .expr19^alpha + .expr23^alpha + .expr27^alpha + 
        .expr31^alpha + .expr55^alpha
    .expr58 <- 1/alpha
    .expr59 <- .expr57^.expr58
    .expr68 <- .expr57^(2/alpha)
    .expr79 <- .expr57^(3/alpha)
    .expr92 <- .expr57^(4/alpha)
    .expr107 <- .expr57^(5/alpha)
    .expr124 <- .expr57^(6/alpha)
    .expr142 <- .expr57^(7/alpha)
    .value <- .expr2^.expr4 * .expr7^.expr4 * .expr11^.expr4 * 
        .expr15^.expr4 * .expr19^.expr4 * .expr23^.expr4 * .expr27^.expr4 * 
        .expr31^.expr4 * (1 + 40320 * alpha^8 + 109584 * alpha^7 * 
        (.expr3 + .expr59) + 118124 * alpha^6 * (1 - 3 * .expr59 + 
        .expr68) + 67284 * alpha^5 * (.expr3 + 7 * .expr59 - 
        6 * .expr68 + .expr79) + 22449 * alpha^4 * (1 - 15 * 
        .expr59 + 25 * .expr68 - 10 * .expr79 + .expr92) + 4536 * 
        alpha^3 * (.expr3 + 31 * .expr59 - 90 * .expr68 + 65 * 
        .expr79 - 15 * .expr92 + .expr107) + 546 * alpha^2 * 
        (1 - 63 * .expr59 + 301 * .expr68 - 350 * .expr79 + 140 * 
            .expr92 - 21 * .expr107 + .expr124) + 36 * alpha * 
        (.expr3 + 127 * .expr59 - 966 * .expr68 + 1701 * .expr79 - 
            1050 * .expr92 + 266 * .expr107 - 28 * .expr124 + 
            .expr142) - 255 * .expr59 + 3025 * .expr68 - 7770 * 
        .expr79 + 6951 * .expr92 - 2646 * .expr107 + 462 * .expr124 - 
        36 * .expr142 + .expr57^(8/alpha)) * .expr57^(-9 + .expr58) * 
        .expr55^.expr4/(E^.expr59 * u1 * u2 * u3 * u4 * u5 * 
        u6 * u7 * u8 * u9)
    .grad <- array(0, c(length(.value), 1), list(NULL, c("s")))
    .grad[, "s"] <- 0
    attr(.value, "gradient") <- .grad
    .value
}, {
    .expr2 <- -log(u1)
    .expr3 <- -1
    .expr4 <- .expr3 + alpha
    .expr7 <- -log(u10)
    .expr11 <- -log(u2)
    .expr15 <- -log(u3)
    .expr19 <- -log(u4)
    .expr23 <- -log(u5)
    .expr27 <- -log(u6)
    .expr31 <- -log(u7)
    .expr35 <- -log(u8)
    .expr61 <- -log(u9)
    .expr63 <- .expr2^alpha + .expr7^alpha + .expr11^alpha + 
        .expr15^alpha + .expr19^alpha + .expr23^alpha + .expr27^alpha + 
        .expr31^alpha + .expr35^alpha + .expr61^alpha
    .expr64 <- 1/alpha
    .expr65 <- .expr63^.expr64
    .expr74 <- .expr63^(2/alpha)
    .expr85 <- .expr63^(3/alpha)
    .expr98 <- .expr63^(4/alpha)
    .expr113 <- .expr63^(5/alpha)
    .expr130 <- .expr63^(6/alpha)
    .expr149 <- .expr63^(7/alpha)
    .expr169 <- .expr63^(8/alpha)
    .value <- .expr2^.expr4 * .expr7^.expr4 * .expr11^.expr4 * 
        .expr15^.expr4 * .expr19^.expr4 * .expr23^.expr4 * .expr27^.expr4 * 
        .expr31^.expr4 * .expr35^.expr4 * (.expr3 + 362880 * 
        alpha^9 + 1026576 * alpha^8 * (.expr3 + .expr65) + 1172700 * 
        alpha^7 * (1 - 3 * .expr65 + .expr74) + 723680 * alpha^6 * 
        (.expr3 + 7 * .expr65 - 6 * .expr74 + .expr85) + 269325 * 
        alpha^5 * (1 - 15 * .expr65 + 25 * .expr74 - 10 * .expr85 + 
        .expr98) + 63273 * alpha^4 * (.expr3 + 31 * .expr65 - 
        90 * .expr74 + 65 * .expr85 - 15 * .expr98 + .expr113) + 
        9450 * alpha^3 * (1 - 63 * .expr65 + 301 * .expr74 - 
            350 * .expr85 + 140 * .expr98 - 21 * .expr113 + .expr130) + 
        870 * alpha^2 * (.expr3 + 127 * .expr65 - 966 * .expr74 + 
            1701 * .expr85 - 1050 * .expr98 + 266 * .expr113 - 
            28 * .expr130 + .expr149) + 45 * alpha * (1 - 255 * 
        .expr65 + 3025 * .expr74 - 7770 * .expr85 + 6951 * .expr98 - 
        2646 * .expr113 + 462 * .expr130 - 36 * .expr149 + .expr169) + 
        511 * .expr65 - 9330 * .expr74 + 34105 * .expr85 - 42525 * 
        .expr98 + 22827 * .expr113 - 5880 * .expr130 + 750 * 
        .expr149 - 45 * .expr169 + .expr63^(9/alpha)) * .expr63^(-10 + 
        .expr64) * .expr61^.expr4/(E^.expr65 * u1 * u10 * u2 * 
        u3 * u4 * u5 * u6 * u7 * u8 * u9)
    .grad <- array(0, c(length(.value), 1), list(NULL, c("s")))
    .grad[, "s"] <- 0
    attr(.value, "gradient") <- .grad
    .value
})
`gumbelCopula.genfun.expr` <-
expression((alpha * (-log(u))^alpha)/(u * log(u)), (alpha * (-1 + 
    alpha - log(u)) * (-log(u))^(-2 + alpha))/u^2)
`gumbelCopula.genfun.algr` <-
expression({
    .expr1 <- log(u)
    .expr2 <- -.expr1
    .expr4 <- alpha * .expr2^alpha
    .expr5 <- u * .expr1
    .expr9 <- 1/u
    .value <- .expr4/.expr5
    .grad <- array(0, c(length(.value), 1), list(NULL, c("u")))
    .grad[, "u"] <- -(alpha * (.expr2^(alpha - 1) * (alpha * 
        .expr9))/.expr5 + .expr4 * (.expr1 + u * .expr9)/.expr5^2)
    attr(.value, "gradient") <- .grad
    .value
}, {
    .expr3 <- log(u)
    .expr5 <- alpha * (-1 + alpha - .expr3)
    .expr6 <- -.expr3
    .expr8 <- -2 + alpha
    .expr9 <- .expr6^.expr8
    .expr10 <- .expr5 * .expr9
    .expr11 <- u^2
    .expr15 <- 1/u
    .value <- .expr10/.expr11
    .grad <- array(0, c(length(.value), 1), list(NULL, c("u")))
    .grad[, "u"] <- -((.expr5 * (.expr6^(.expr8 - 1) * (.expr8 * 
        .expr15)) + alpha * .expr15 * .expr9)/.expr11 + .expr10 * 
        (2 * u)/.expr11^2)
    attr(.value, "gradient") <- .grad
    .value
})
