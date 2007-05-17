`frankCopula.pdf.expr` <-
expression(1, (alpha * E^(alpha * (1 + u1 + u2)) * (-1 + E^alpha))/(E^alpha - 
    E^(alpha + alpha * u1) + E^(alpha * (u1 + u2)) - E^(alpha + 
    alpha * u2))^2, (alpha^2 * E^(alpha * (-2 + u1 + u2 + u3)) * 
    (-1 + E^alpha)^2 * (-1 + E^(alpha * u1) + E^(alpha * u2) - 
    E^(alpha * (u1 + u2)) + E^(alpha * u3) - E^(alpha * (u1 + 
    u3)) - E^(alpha * (u2 + u3)) + E^(alpha * (-2 + u1 + u2 + 
    u3)) - 2 * E^(alpha * (-1 + u1 + u2 + u3)) + 2 * E^(alpha * 
    (u1 + u2 + u3))))/(1 - E^(alpha * u1) - E^(alpha * u2) + 
    E^(alpha * (u1 + u2)) - E^(alpha * u3) + E^(alpha * (u1 + 
    u3)) + E^(alpha * (u2 + u3)) + E^(alpha * (-2 + u1 + u2 + 
    u3)) - 2 * E^(alpha * (-1 + u1 + u2 + u3)))^3, (alpha^3 * 
    E^(alpha * (1 + u1 + u2 + u3 + u4)) * (-1 + E^(-alpha))^4 * 
    (-1 + E^alpha) * (1 - 2 * E^alpha + E^(2 * alpha) + (E^(2 * 
    alpha * (-3 + u1 + u2 + u3 + u4)) * (-1 + E^alpha)^8)/((-1 + 
    E^(alpha * u1))^2 * (-1 + E^(alpha * u2))^2 * (-1 + E^(alpha * 
    u3))^2 * (-1 + E^(alpha * u4))^2) - (4 * E^(alpha * (-3 + 
    u1 + u2 + u3 + u4)) * (-1 + E^alpha)^4)/((-1 + E^(alpha * 
    u1)) * (-1 + E^(alpha * u2)) * (-1 + E^(alpha * u3)) * (-1 + 
    E^(alpha * u4))) + (4 * E^(alpha * (-2 + u1 + u2 + u3 + u4)) * 
    (-1 + E^alpha)^4)/((-1 + E^(alpha * u1)) * (-1 + E^(alpha * 
    u2)) * (-1 + E^(alpha * u3)) * (-1 + E^(alpha * u4)))))/((-1 + 
    E^(alpha * u1))^2 * (-1 + E^(alpha * u2))^2 * (-1 + E^(alpha * 
    u3))^2 * (-1 + E^(alpha * u4))^2 * (1 - E^alpha + (E^(alpha * 
    (-3 + u1 + u2 + u3 + u4)) * (-1 + E^alpha)^4)/((-1 + E^(alpha * 
    u1)) * (-1 + E^(alpha * u2)) * (-1 + E^(alpha * u3)) * (-1 + 
    E^(alpha * u4))))^4), -((alpha^4 * E^(alpha * (1 + u1 + u2 + 
    u3 + u4 + u5)) * (-1 + E^(-alpha))^5 * (-1 + E^alpha) * (-1 + 
    3 * E^alpha - 3 * E^(2 * alpha) + E^(3 * alpha) + (11 * E^(3 * 
    alpha) * (-1 + E^(-alpha))^10)/((-1 + E^(-(alpha * u1)))^2 * 
    (-1 + E^(-(alpha * u2)))^2 * (-1 + E^(-(alpha * u3)))^2 * 
    (-1 + E^(-(alpha * u4)))^2 * (-1 + E^(-(alpha * u5)))^2) + 
    (E^(3 * alpha * (-4 + u1 + u2 + u3 + u4 + u5)) * (-1 + E^alpha)^15)/((-1 + 
        E^(alpha * u1))^3 * (-1 + E^(alpha * u2))^3 * (-1 + E^(alpha * 
        u3))^3 * (-1 + E^(alpha * u4))^3 * (-1 + E^(alpha * u5))^3) - 
    (11 * E^(2 * alpha * (-4 + u1 + u2 + u3 + u4 + u5)) * (-1 + 
        E^alpha)^10)/((-1 + E^(alpha * u1))^2 * (-1 + E^(alpha * 
        u2))^2 * (-1 + E^(alpha * u3))^2 * (-1 + E^(alpha * u4))^2 * 
        (-1 + E^(alpha * u5))^2) + (11 * E^(alpha * (-4 + u1 + 
    u2 + u3 + u4 + u5)) * (-1 + E^alpha)^5)/((-1 + E^(alpha * 
    u1)) * (-1 + E^(alpha * u2)) * (-1 + E^(alpha * u3)) * (-1 + 
    E^(alpha * u4)) * (-1 + E^(alpha * u5))) - (22 * E^(alpha * 
    (-3 + u1 + u2 + u3 + u4 + u5)) * (-1 + E^alpha)^5)/((-1 + 
    E^(alpha * u1)) * (-1 + E^(alpha * u2)) * (-1 + E^(alpha * 
    u3)) * (-1 + E^(alpha * u4)) * (-1 + E^(alpha * u5))) + (11 * 
    E^(alpha * (-2 + u1 + u2 + u3 + u4 + u5)) * (-1 + E^alpha)^5)/((-1 + 
    E^(alpha * u1)) * (-1 + E^(alpha * u2)) * (-1 + E^(alpha * 
    u3)) * (-1 + E^(alpha * u4)) * (-1 + E^(alpha * u5)))))/((-1 + 
    E^(alpha * u1))^2 * (-1 + E^(alpha * u2))^2 * (-1 + E^(alpha * 
    u3))^2 * (-1 + E^(alpha * u4))^2 * (-1 + E^(alpha * u5))^2 * 
    (1 - E^alpha + (E^(alpha * (-4 + u1 + u2 + u3 + u4 + u5)) * 
        (-1 + E^alpha)^5)/((-1 + E^(alpha * u1)) * (-1 + E^(alpha * 
        u2)) * (-1 + E^(alpha * u3)) * (-1 + E^(alpha * u4)) * 
        (-1 + E^(alpha * u5))))^5)), (alpha^5 * E^(alpha * (1 + 
    u1 + u2 + u3 + u4 + u5 + u6)) * (-1 + E^(-alpha))^6 * (-1 + 
    E^alpha) * (1 - 4 * E^alpha + 6 * E^(2 * alpha) - 4 * E^(3 * 
    alpha) + E^(4 * alpha) + (26 * E^(4 * alpha) * (-1 + E^(-alpha))^18)/((-1 + 
    E^(-(alpha * u1)))^3 * (-1 + E^(-(alpha * u2)))^3 * (-1 + 
    E^(-(alpha * u3)))^3 * (-1 + E^(-(alpha * u4)))^3 * (-1 + 
    E^(-(alpha * u5)))^3 * (-1 + E^(-(alpha * u6)))^3) - (132 * 
    E^(3 * alpha) * (-1 + E^(-alpha))^12)/((-1 + E^(-(alpha * 
    u1)))^2 * (-1 + E^(-(alpha * u2)))^2 * (-1 + E^(-(alpha * 
    u3)))^2 * (-1 + E^(-(alpha * u4)))^2 * (-1 + E^(-(alpha * 
    u5)))^2 * (-1 + E^(-(alpha * u6)))^2) + (E^(4 * alpha * (-5 + 
    u1 + u2 + u3 + u4 + u5 + u6)) * (-1 + E^alpha)^24)/((-1 + 
    E^(alpha * u1))^4 * (-1 + E^(alpha * u2))^4 * (-1 + E^(alpha * 
    u3))^4 * (-1 + E^(alpha * u4))^4 * (-1 + E^(alpha * u5))^4 * 
    (-1 + E^(alpha * u6))^4) - (26 * E^(3 * alpha * (-5 + u1 + 
    u2 + u3 + u4 + u5 + u6)) * (-1 + E^alpha)^18)/((-1 + E^(alpha * 
    u1))^3 * (-1 + E^(alpha * u2))^3 * (-1 + E^(alpha * u3))^3 * 
    (-1 + E^(alpha * u4))^3 * (-1 + E^(alpha * u5))^3 * (-1 + 
    E^(alpha * u6))^3) + (66 * E^(2 * alpha * (-5 + u1 + u2 + 
    u3 + u4 + u5 + u6)) * (-1 + E^alpha)^12)/((-1 + E^(alpha * 
    u1))^2 * (-1 + E^(alpha * u2))^2 * (-1 + E^(alpha * u3))^2 * 
    (-1 + E^(alpha * u4))^2 * (-1 + E^(alpha * u5))^2 * (-1 + 
    E^(alpha * u6))^2) + (66 * E^(2 * alpha * (-4 + u1 + u2 + 
    u3 + u4 + u5 + u6)) * (-1 + E^alpha)^12)/((-1 + E^(alpha * 
    u1))^2 * (-1 + E^(alpha * u2))^2 * (-1 + E^(alpha * u3))^2 * 
    (-1 + E^(alpha * u4))^2 * (-1 + E^(alpha * u5))^2 * (-1 + 
    E^(alpha * u6))^2) - (26 * E^(alpha * (-5 + u1 + u2 + u3 + 
    u4 + u5 + u6)) * (-1 + E^alpha)^6)/((-1 + E^(alpha * u1)) * 
    (-1 + E^(alpha * u2)) * (-1 + E^(alpha * u3)) * (-1 + E^(alpha * 
    u4)) * (-1 + E^(alpha * u5)) * (-1 + E^(alpha * u6))) + (78 * 
    E^(alpha * (-4 + u1 + u2 + u3 + u4 + u5 + u6)) * (-1 + E^alpha)^6)/((-1 + 
    E^(alpha * u1)) * (-1 + E^(alpha * u2)) * (-1 + E^(alpha * 
    u3)) * (-1 + E^(alpha * u4)) * (-1 + E^(alpha * u5)) * (-1 + 
    E^(alpha * u6))) - (78 * E^(alpha * (-3 + u1 + u2 + u3 + 
    u4 + u5 + u6)) * (-1 + E^alpha)^6)/((-1 + E^(alpha * u1)) * 
    (-1 + E^(alpha * u2)) * (-1 + E^(alpha * u3)) * (-1 + E^(alpha * 
    u4)) * (-1 + E^(alpha * u5)) * (-1 + E^(alpha * u6))) + (26 * 
    E^(alpha * (-2 + u1 + u2 + u3 + u4 + u5 + u6)) * (-1 + E^alpha)^6)/((-1 + 
    E^(alpha * u1)) * (-1 + E^(alpha * u2)) * (-1 + E^(alpha * 
    u3)) * (-1 + E^(alpha * u4)) * (-1 + E^(alpha * u5)) * (-1 + 
    E^(alpha * u6)))))/((-1 + E^(alpha * u1))^2 * (-1 + E^(alpha * 
    u2))^2 * (-1 + E^(alpha * u3))^2 * (-1 + E^(alpha * u4))^2 * 
    (-1 + E^(alpha * u5))^2 * (-1 + E^(alpha * u6))^2 * (1 - 
    E^alpha + (E^(alpha * (-5 + u1 + u2 + u3 + u4 + u5 + u6)) * 
    (-1 + E^alpha)^6)/((-1 + E^(alpha * u1)) * (-1 + E^(alpha * 
    u2)) * (-1 + E^(alpha * u3)) * (-1 + E^(alpha * u4)) * (-1 + 
    E^(alpha * u5)) * (-1 + E^(alpha * u6))))^6))
`frankCopula.pdf.algr` <-
expression({
    .value <- 1
    .grad <- array(0, c(length(.value), 1), list(NULL, c("s")))
    .grad[, "s"] <- 0
    attr(.value, "gradient") <- .grad
    .value
}, {
    .expr7 <- E^alpha
    .value <- alpha * E^(alpha * (1 + u1 + u2)) * (-1 + .expr7)/(.expr7 - 
        E^(alpha + alpha * u1) + E^(alpha * (u1 + u2)) - E^(alpha + 
        alpha * u2))^2
    .grad <- array(0, c(length(.value), 1), list(NULL, c("s")))
    .grad[, "s"] <- 0
    attr(.value, "gradient") <- .grad
    .value
}, {
    .expr7 <- E^(alpha * (-2 + u1 + u2 + u3))
    .expr9 <- -1
    .expr15 <- E^(alpha * u1)
    .expr18 <- E^(alpha * u2)
    .expr20 <- u1 + u2
    .expr22 <- E^(alpha * .expr20)
    .expr25 <- E^(alpha * u3)
    .expr29 <- E^(alpha * (u1 + u3))
    .expr33 <- E^(alpha * (u2 + u3))
    .expr41 <- 2 * E^(alpha * (.expr9 + u1 + u2 + u3))
    .value <- alpha^2 * .expr7 * (.expr9 + E^alpha)^2 * (.expr9 + 
        .expr15 + .expr18 - .expr22 + .expr25 - .expr29 - .expr33 + 
        .expr7 - .expr41 + 2 * E^(alpha * (.expr20 + u3)))/(1 - 
        .expr15 - .expr18 + .expr22 - .expr25 + .expr29 + .expr33 + 
        .expr7 - .expr41)^3
    .grad <- array(0, c(length(.value), 1), list(NULL, c("s")))
    .grad[, "s"] <- 0
    attr(.value, "gradient") <- .grad
    .value
}, {
    .expr9 <- -1
    .expr15 <- E^alpha
    .expr16 <- .expr9 + .expr15
    .expr20 <- 2 * alpha
    .expr27 <- -3 + u1 + u2 + u3 + u4
    .expr34 <- .expr9 + E^(alpha * u1)
    .expr38 <- .expr9 + E^(alpha * u2)
    .expr43 <- .expr9 + E^(alpha * u3)
    .expr48 <- .expr9 + E^(alpha * u4)
    .expr50 <- .expr34^2 * .expr38^2 * .expr43^2 * .expr48^2
    .expr54 <- E^(alpha * .expr27)
    .expr56 <- .expr16^4
    .expr60 <- .expr34 * .expr38 * .expr43 * .expr48
    .value <- alpha^3 * E^(alpha * (1 + u1 + u2 + u3 + u4)) * 
        (.expr9 + E^-alpha)^4 * .expr16 * (1 - 2 * .expr15 + 
        E^.expr20 + E^(.expr20 * .expr27) * .expr16^8/.expr50 - 
        4 * .expr54 * .expr56/.expr60 + 4 * E^(alpha * (-2 + 
        u1 + u2 + u3 + u4)) * .expr56/.expr60)/(.expr50 * (1 - 
        .expr15 + .expr54 * .expr56/.expr60)^4)
    .grad <- array(0, c(length(.value), 1), list(NULL, c("s")))
    .grad[, "s"] <- 0
    attr(.value, "gradient") <- .grad
    .value
}, {
    .expr10 <- -1
    .expr13 <- .expr10 + E^-alpha
    .expr16 <- E^alpha
    .expr17 <- .expr10 + .expr16
    .expr21 <- 2 * alpha
    .expr25 <- 3 * alpha
    .expr26 <- E^.expr25
    .expr31 <- alpha * u1
    .expr36 <- alpha * u2
    .expr42 <- alpha * u3
    .expr48 <- alpha * u4
    .expr54 <- alpha * u5
    .expr67 <- -4 + u1 + u2 + u3 + u4 + u5
    .expr73 <- .expr10 + E^.expr31
    .expr76 <- .expr10 + E^.expr36
    .expr80 <- .expr10 + E^.expr42
    .expr84 <- .expr10 + E^.expr48
    .expr88 <- .expr10 + E^.expr54
    .expr106 <- .expr73^2 * .expr76^2 * .expr80^2 * .expr84^2 * 
        .expr88^2
    .expr110 <- E^(alpha * .expr67)
    .expr112 <- .expr17^5
    .expr117 <- .expr73 * .expr76 * .expr80 * .expr84 * .expr88
    .value <- -(alpha^4 * E^(alpha * (1 + u1 + u2 + u3 + u4 + 
        u5)) * .expr13^5 * .expr17 * (.expr10 + 3 * .expr16 - 
        3 * E^.expr21 + .expr26 + 11 * .expr26 * .expr13^10/((.expr10 + 
        E^-.expr31)^2 * (.expr10 + E^-.expr36)^2 * (.expr10 + 
        E^-.expr42)^2 * (.expr10 + E^-.expr48)^2 * (.expr10 + 
        E^-.expr54)^2) + E^(.expr25 * .expr67) * .expr17^15/(.expr73^3 * 
        .expr76^3 * .expr80^3 * .expr84^3 * .expr88^3) - 11 * 
        E^(.expr21 * .expr67) * .expr17^10/.expr106 + 11 * .expr110 * 
        .expr112/.expr117 - 22 * E^(alpha * (-3 + u1 + u2 + u3 + 
        u4 + u5)) * .expr112/.expr117 + 11 * E^(alpha * (-2 + 
        u1 + u2 + u3 + u4 + u5)) * .expr112/.expr117)/(.expr106 * 
        (1 - .expr16 + .expr110 * .expr112/.expr117)^5))
    .grad <- array(0, c(length(.value), 1), list(NULL, c("s")))
    .grad[, "s"] <- 0
    attr(.value, "gradient") <- .grad
    .value
}, {
    .expr11 <- -1
    .expr14 <- .expr11 + E^-alpha
    .expr17 <- E^alpha
    .expr18 <- .expr11 + .expr17
    .expr22 <- 2 * alpha
    .expr26 <- 3 * alpha
    .expr27 <- E^.expr26
    .expr30 <- 4 * alpha
    .expr31 <- E^.expr30
    .expr36 <- alpha * u1
    .expr39 <- .expr11 + E^-.expr36
    .expr41 <- alpha * u2
    .expr44 <- .expr11 + E^-.expr41
    .expr47 <- alpha * u3
    .expr50 <- .expr11 + E^-.expr47
    .expr53 <- alpha * u4
    .expr56 <- .expr11 + E^-.expr53
    .expr59 <- alpha * u5
    .expr62 <- .expr11 + E^-.expr59
    .expr65 <- alpha * u6
    .expr68 <- .expr11 + E^-.expr65
    .expr95 <- -5 + u1 + u2 + u3 + u4 + u5 + u6
    .expr101 <- .expr11 + E^.expr36
    .expr104 <- .expr11 + E^.expr41
    .expr108 <- .expr11 + E^.expr47
    .expr112 <- .expr11 + E^.expr53
    .expr116 <- .expr11 + E^.expr59
    .expr120 <- .expr11 + E^.expr65
    .expr146 <- .expr18^12
    .expr158 <- .expr101^2 * .expr104^2 * .expr108^2 * .expr112^2 * 
        .expr116^2 * .expr120^2
    .expr167 <- -4 + u1 + u2 + u3 + u4 + u5 + u6
    .expr175 <- E^(alpha * .expr95)
    .expr177 <- .expr18^6
    .expr183 <- .expr101 * .expr104 * .expr108 * .expr112 * .expr116 * 
        .expr120
    .value <- alpha^5 * E^(alpha * (1 + u1 + u2 + u3 + u4 + u5 + 
        u6)) * .expr14^6 * .expr18 * (1 - 4 * .expr17 + 6 * E^.expr22 - 
        4 * .expr27 + .expr31 + 26 * .expr31 * .expr14^18/(.expr39^3 * 
        .expr44^3 * .expr50^3 * .expr56^3 * .expr62^3 * .expr68^3) - 
        132 * .expr27 * .expr14^12/(.expr39^2 * .expr44^2 * .expr50^2 * 
            .expr56^2 * .expr62^2 * .expr68^2) + E^(.expr30 * 
        .expr95) * .expr18^24/(.expr101^4 * .expr104^4 * .expr108^4 * 
        .expr112^4 * .expr116^4 * .expr120^4) - 26 * E^(.expr26 * 
        .expr95) * .expr18^18/(.expr101^3 * .expr104^3 * .expr108^3 * 
        .expr112^3 * .expr116^3 * .expr120^3) + 66 * E^(.expr22 * 
        .expr95) * .expr146/.expr158 + 66 * E^(.expr22 * .expr167) * 
        .expr146/.expr158 - 26 * .expr175 * .expr177/.expr183 + 
        78 * E^(alpha * .expr167) * .expr177/.expr183 - 78 * 
        E^(alpha * (-3 + u1 + u2 + u3 + u4 + u5 + u6)) * .expr177/.expr183 + 
        26 * E^(alpha * (-2 + u1 + u2 + u3 + u4 + u5 + u6)) * 
            .expr177/.expr183)/(.expr158 * (1 - .expr17 + .expr175 * 
        .expr177/.expr183)^6)
    .grad <- array(0, c(length(.value), 1), list(NULL, c("s")))
    .grad[, "s"] <- 0
    attr(.value, "gradient") <- .grad
    .value
})
`frankCopula.genfun.expr` <-
expression(alpha/(1 - E^(alpha * u)), (alpha^2 * E^(alpha * u))/(-1 + 
    E^(alpha * u))^2)
`frankCopula.genfun.algr` <-
expression({
    .expr2 <- E^(alpha * u)
    .expr3 <- 1 - .expr2
    .value <- alpha/.expr3
    .grad <- array(0, c(length(.value), 1), list(NULL, c("u")))
    .grad[, "u"] <- alpha * (.expr2 * (log(E) * alpha))/.expr3^2
    attr(.value, "gradient") <- .grad
    .value
}, {
    .expr1 <- alpha^2
    .expr3 <- E^(alpha * u)
    .expr4 <- .expr1 * .expr3
    .expr6 <- -1 + .expr3
    .expr7 <- .expr6^2
    .expr11 <- .expr3 * (log(E) * alpha)
    .value <- .expr4/.expr7
    .grad <- array(0, c(length(.value), 1), list(NULL, c("u")))
    .grad[, "u"] <- .expr1 * .expr11/.expr7 - .expr4 * (2 * (.expr11 * 
        .expr6))/.expr7^2
    attr(.value, "gradient") <- .grad
    .value
})
