getExpr<- function(fun, param) {
  n <- length(param)
  if (n == 0) {
    expr <- parse(text = paste(fun, "(", "x", ")"))
  }
  else {
    nm <- names(param)
    expr <- paste(nm, "=", param, collapse=", ")
    expr <- parse(text = paste(fun, "(", "x", ",", expr, ")"))
  }
  ##eval(expr)
  expr
}

dmvdc <- function(mvdc, x) {
  dim <- mvdc@copula@dimension
  densmarg <- 1
  if (is.vector(x)) x <- matrix(x, nrow = 1)
  u <- x
  for (i in 1:dim) {
    pmarg <- paste("p", mvdc@margins[i], sep = "")
    dmarg <- paste("d", mvdc@margins[i], sep = "")
    cdf.expr <- getExpr(pmarg, mvdc@paramMargins[[i]])
    u[,i] <- eval(cdf.expr, list(x = x[,i]))
    pdf.expr <- getExpr(dmarg, mvdc@paramMargins[[i]])
    densmarg <- densmarg * eval(pdf.expr, list(x = x[,i]))
  }
  dcopula(mvdc@copula, u) * densmarg
}


pmvdc <- function(mvdc, x) {
  dim <- mvdc@copula@dimension
  if (is.vector(x)) x <- matrix(x, nrow = 1)
  u <- x
  for (i in 1:dim) {
    pmarg <- paste("p", mvdc@margins[i], sep = "")
    cdf.expr <- getExpr(pmarg, mvdc@paramMargins[[i]])
    u[,i] <- eval(cdf.expr, list(x = x[,i]))
  }
  pcopula(mvdc@copula, u)
}


rmvdc <- function(mvdc, n) {
  dim <- mvdc@copula@dimension
  u <- rcopula(mvdc@copula, n)
  x <- u
  for (i in 1:dim) {
    qmarg <- paste("q", mvdc@margins[i], sep = "")
    qdf.expr <- getExpr(qmarg, mvdc@paramMargins[[i]])
    x[,i] <- eval(qdf.expr, list(x = u[,i]))
  }
  x
}
