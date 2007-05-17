`claytonCopula.pdf.expr` <-
expression(1/(u1 * (u1^(-alpha))^(1/alpha)), (1 + alpha) * u1^(-1 - 
    alpha) * u2^(-1 - alpha) * (-1 + u1^(-alpha) + u2^(-alpha))^(-2 - 
    1/alpha), (1 + alpha) * (1 + 2 * alpha) * u1^(-1 - alpha) * 
    u2^(-1 - alpha) * u3^(-1 - alpha) * (-2 + u1^(-alpha) + u2^(-alpha) + 
    u3^(-alpha))^(-3 - 1/alpha), (1 + alpha) * (1 + 2 * alpha) * 
    (1 + 3 * alpha) * u1^(-1 - alpha) * u2^(-1 - alpha) * u3^(-1 - 
    alpha) * u4^(-1 - alpha) * (-3 + u1^(-alpha) + u2^(-alpha) + 
    u3^(-alpha) + u4^(-alpha))^(-4 - 1/alpha), (1 + alpha) * 
    (1 + 2 * alpha) * (1 + 3 * alpha) * (1 + 4 * alpha) * u1^(-1 - 
    alpha) * u2^(-1 - alpha) * u3^(-1 - alpha) * u4^(-1 - alpha) * 
    u5^(-1 - alpha) * (-4 + u1^(-alpha) + u2^(-alpha) + u3^(-alpha) + 
    u4^(-alpha) + u5^(-alpha))^(-5 - 1/alpha), (1 + alpha) * 
    (1 + 2 * alpha) * (1 + 3 * alpha) * (1 + 4 * alpha) * (1 + 
    5 * alpha) * u1^(-1 - alpha) * u2^(-1 - alpha) * u3^(-1 - 
    alpha) * u4^(-1 - alpha) * u5^(-1 - alpha) * u6^(-1 - alpha) * 
    (-5 + u1^(-alpha) + u2^(-alpha) + u3^(-alpha) + u4^(-alpha) + 
        u5^(-alpha) + u6^(-alpha))^(-6 - 1/alpha), (1 + alpha) * 
    (1 + 2 * alpha) * (1 + 3 * alpha) * (1 + 4 * alpha) * (1 + 
    5 * alpha) * (1 + 6 * alpha) * u1^(-1 - alpha) * u2^(-1 - 
    alpha) * u3^(-1 - alpha) * u4^(-1 - alpha) * u5^(-1 - alpha) * 
    u6^(-1 - alpha) * u7^(-1 - alpha) * (-6 + u1^(-alpha) + u2^(-alpha) + 
    u3^(-alpha) + u4^(-alpha) + u5^(-alpha) + u6^(-alpha) + u7^(-alpha))^(-7 - 
    1/alpha), (1 + alpha) * (1 + 2 * alpha) * (1 + 3 * alpha) * 
    (1 + 4 * alpha) * (1 + 5 * alpha) * (1 + 6 * alpha) * (1 + 
    7 * alpha) * u1^(-1 - alpha) * u2^(-1 - alpha) * u3^(-1 - 
    alpha) * u4^(-1 - alpha) * u5^(-1 - alpha) * u6^(-1 - alpha) * 
    u7^(-1 - alpha) * u8^(-1 - alpha) * (-7 + u1^(-alpha) + u2^(-alpha) + 
    u3^(-alpha) + u4^(-alpha) + u5^(-alpha) + u6^(-alpha) + u7^(-alpha) + 
    u8^(-alpha))^(-8 - 1/alpha), (1 + alpha) * (1 + 2 * alpha) * 
    (1 + 3 * alpha) * (1 + 4 * alpha) * (1 + 5 * alpha) * (1 + 
    6 * alpha) * (1 + 7 * alpha) * (1 + 8 * alpha) * u1^(-1 - 
    alpha) * u2^(-1 - alpha) * u3^(-1 - alpha) * u4^(-1 - alpha) * 
    u5^(-1 - alpha) * u6^(-1 - alpha) * u7^(-1 - alpha) * u8^(-1 - 
    alpha) * u9^(-1 - alpha) * (-8 + u1^(-alpha) + u2^(-alpha) + 
    u3^(-alpha) + u4^(-alpha) + u5^(-alpha) + u6^(-alpha) + u7^(-alpha) + 
    u8^(-alpha) + u9^(-alpha))^(-9 - 1/alpha), (1 + alpha) * 
    (1 + 2 * alpha) * (1 + 3 * alpha) * (1 + 4 * alpha) * (1 + 
    5 * alpha) * (1 + 6 * alpha) * (1 + 7 * alpha) * (1 + 8 * 
    alpha) * (1 + 9 * alpha) * u1^(-1 - alpha) * u10^(-1 - alpha) * 
    u2^(-1 - alpha) * u3^(-1 - alpha) * u4^(-1 - alpha) * u5^(-1 - 
    alpha) * u6^(-1 - alpha) * u7^(-1 - alpha) * u8^(-1 - alpha) * 
    u9^(-1 - alpha) * (-9 + u1^(-alpha) + u10^(-alpha) + u2^(-alpha) + 
    u3^(-alpha) + u4^(-alpha) + u5^(-alpha) + u6^(-alpha) + u7^(-alpha) + 
    u8^(-alpha) + u9^(-alpha))^(-10 - 1/alpha))
`claytonCopula.pdf.algr` <-
expression({
    .value <- 1/(u1 * (u1^-alpha)^(1/alpha))
    .grad <- array(0, c(length(.value), 1), list(NULL, c("s")))
    .grad[, "s"] <- 0
    attr(.value, "gradient") <- .grad
    .value
}, {
    .expr2 <- -1
    .expr3 <- .expr2 - alpha
    .expr8 <- -alpha
    .value <- (1 + alpha) * u1^.expr3 * u2^.expr3 * (.expr2 + 
        u1^.expr8 + u2^.expr8)^(-2 - 1/alpha)
    .grad <- array(0, c(length(.value), 1), list(NULL, c("s")))
    .grad[, "s"] <- 0
    attr(.value, "gradient") <- .grad
    .value
}, {
    .expr6 <- -1 - alpha
    .expr14 <- -alpha
    .value <- (1 + alpha) * (1 + 2 * alpha) * u1^.expr6 * u2^.expr6 * 
        u3^.expr6 * (-2 + u1^.expr14 + u2^.expr14 + u3^.expr14)^(-3 - 
        1/alpha)
    .grad <- array(0, c(length(.value), 1), list(NULL, c("s")))
    .grad[, "s"] <- 0
    attr(.value, "gradient") <- .grad
    .value
}, {
    .expr9 <- -1 - alpha
    .expr19 <- -alpha
    .value <- (1 + alpha) * (1 + 2 * alpha) * (1 + 3 * alpha) * 
        u1^.expr9 * u2^.expr9 * u3^.expr9 * u4^.expr9 * (-3 + 
        u1^.expr19 + u2^.expr19 + u3^.expr19 + u4^.expr19)^(-4 - 
        1/alpha)
    .grad <- array(0, c(length(.value), 1), list(NULL, c("s")))
    .grad[, "s"] <- 0
    attr(.value, "gradient") <- .grad
    .value
}, {
    .expr12 <- -1 - alpha
    .expr24 <- -alpha
    .value <- (1 + alpha) * (1 + 2 * alpha) * (1 + 3 * alpha) * 
        (1 + 4 * alpha) * u1^.expr12 * u2^.expr12 * u3^.expr12 * 
        u4^.expr12 * u5^.expr12 * (-4 + u1^.expr24 + u2^.expr24 + 
        u3^.expr24 + u4^.expr24 + u5^.expr24)^(-5 - 1/alpha)
    .grad <- array(0, c(length(.value), 1), list(NULL, c("s")))
    .grad[, "s"] <- 0
    attr(.value, "gradient") <- .grad
    .value
}, {
    .expr15 <- -1 - alpha
    .expr29 <- -alpha
    .value <- (1 + alpha) * (1 + 2 * alpha) * (1 + 3 * alpha) * 
        (1 + 4 * alpha) * (1 + 5 * alpha) * u1^.expr15 * u2^.expr15 * 
        u3^.expr15 * u4^.expr15 * u5^.expr15 * u6^.expr15 * (-5 + 
        u1^.expr29 + u2^.expr29 + u3^.expr29 + u4^.expr29 + u5^.expr29 + 
        u6^.expr29)^(-6 - 1/alpha)
    .grad <- array(0, c(length(.value), 1), list(NULL, c("s")))
    .grad[, "s"] <- 0
    attr(.value, "gradient") <- .grad
    .value
}, {
    .expr18 <- -1 - alpha
    .expr34 <- -alpha
    .value <- (1 + alpha) * (1 + 2 * alpha) * (1 + 3 * alpha) * 
        (1 + 4 * alpha) * (1 + 5 * alpha) * (1 + 6 * alpha) * 
        u1^.expr18 * u2^.expr18 * u3^.expr18 * u4^.expr18 * u5^.expr18 * 
        u6^.expr18 * u7^.expr18 * (-6 + u1^.expr34 + u2^.expr34 + 
        u3^.expr34 + u4^.expr34 + u5^.expr34 + u6^.expr34 + u7^.expr34)^(-7 - 
        1/alpha)
    .grad <- array(0, c(length(.value), 1), list(NULL, c("s")))
    .grad[, "s"] <- 0
    attr(.value, "gradient") <- .grad
    .value
}, {
    .expr21 <- -1 - alpha
    .expr39 <- -alpha
    .value <- (1 + alpha) * (1 + 2 * alpha) * (1 + 3 * alpha) * 
        (1 + 4 * alpha) * (1 + 5 * alpha) * (1 + 6 * alpha) * 
        (1 + 7 * alpha) * u1^.expr21 * u2^.expr21 * u3^.expr21 * 
        u4^.expr21 * u5^.expr21 * u6^.expr21 * u7^.expr21 * u8^.expr21 * 
        (-7 + u1^.expr39 + u2^.expr39 + u3^.expr39 + u4^.expr39 + 
            u5^.expr39 + u6^.expr39 + u7^.expr39 + u8^.expr39)^(-8 - 
            1/alpha)
    .grad <- array(0, c(length(.value), 1), list(NULL, c("s")))
    .grad[, "s"] <- 0
    attr(.value, "gradient") <- .grad
    .value
}, {
    .expr24 <- -1 - alpha
    .expr44 <- -alpha
    .value <- (1 + alpha) * (1 + 2 * alpha) * (1 + 3 * alpha) * 
        (1 + 4 * alpha) * (1 + 5 * alpha) * (1 + 6 * alpha) * 
        (1 + 7 * alpha) * (1 + 8 * alpha) * u1^.expr24 * u2^.expr24 * 
        u3^.expr24 * u4^.expr24 * u5^.expr24 * u6^.expr24 * u7^.expr24 * 
        u8^.expr24 * u9^.expr24 * (-8 + u1^.expr44 + u2^.expr44 + 
        u3^.expr44 + u4^.expr44 + u5^.expr44 + u6^.expr44 + u7^.expr44 + 
        u8^.expr44 + u9^.expr44)^(-9 - 1/alpha)
    .grad <- array(0, c(length(.value), 1), list(NULL, c("s")))
    .grad[, "s"] <- 0
    attr(.value, "gradient") <- .grad
    .value
}, {
    .expr27 <- -1 - alpha
    .expr49 <- -alpha
    .value <- (1 + alpha) * (1 + 2 * alpha) * (1 + 3 * alpha) * 
        (1 + 4 * alpha) * (1 + 5 * alpha) * (1 + 6 * alpha) * 
        (1 + 7 * alpha) * (1 + 8 * alpha) * (1 + 9 * alpha) * 
        u1^.expr27 * u10^.expr27 * u2^.expr27 * u3^.expr27 * 
        u4^.expr27 * u5^.expr27 * u6^.expr27 * u7^.expr27 * u8^.expr27 * 
        u9^.expr27 * (-9 + u1^.expr49 + u10^.expr49 + u2^.expr49 + 
        u3^.expr49 + u4^.expr49 + u5^.expr49 + u6^.expr49 + u7^.expr49 + 
        u8^.expr49 + u9^.expr49)^(-10 - 1/alpha)
    .grad <- array(0, c(length(.value), 1), list(NULL, c("s")))
    .grad[, "s"] <- 0
    attr(.value, "gradient") <- .grad
    .value
})
`claytonCopula.genfun.expr` <-
expression(-(alpha * u^(-1 - alpha)), alpha * (1 + alpha) * u^(-2 - 
    alpha))
`claytonCopula.genfun.algr` <-
expression({
    .expr2 <- -1 - alpha
    .value <- -(alpha * u^.expr2)
    .grad <- array(0, c(length(.value), 1), list(NULL, c("u")))
    .grad[, "u"] <- -(alpha * (u^(.expr2 - 1) * .expr2))
    attr(.value, "gradient") <- .grad
    .value
}, {
    .expr2 <- alpha * (1 + alpha)
    .expr4 <- -2 - alpha
    .value <- .expr2 * u^.expr4
    .grad <- array(0, c(length(.value), 1), list(NULL, c("u")))
    .grad[, "u"] <- .expr2 * (u^(.expr4 - 1) * .expr4)
    attr(.value, "gradient") <- .grad
    .value
})
