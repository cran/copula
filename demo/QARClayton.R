# MLE for Clayton AR(1) copula model with Student marginal

# Parameters

	alpha <- 10
	mu <- 0
	sigma <- 1
	df <- 3

# Student location scale rpqd family

rtls <- function(n,df,mu,sigma) sigma * rt(n,df) + mu
ptls <- function(x,df,mu,sigma) pt((x - mu)/sigma,df)
qtls <- function(u,df,mu,sigma) sigma * qt(u,df) + mu
dtls <- function(u,df,mu,sigma) dt((x - mu)/sigma,df)/sigma

#Generation of Data

n <- 200
rclayton <- function(n,alpha){
	u <- runif(n+1) #innovations
	v <- u
	for(i in 2:(n+1))
		v[i] <- ((u[i]^(-alpha/(1+alpha)) -1)*v[i-1]^(-alpha) +1)^(-1/alpha)
	return(v[2:(n+1)])
	}
u <- rclayton(n,alpha)
u <- qt(u,3)
y <- u[-n]
x <- u[-1]

plot(x,y)

# Estimation with known marginal
cop <- claytonCopula(5,dim=2)
f <- fitCopula(cbind(ptls(x,df,mu,sigma),ptls(y,df,mu,sigma)),cop)

# Estimation with unknown marginal parameters

M1 <- mvdc(claytonCopula(alpha,dim=2),c("tls","tls"),
	list( list(df=df,mu=mu,sigma=sigma), list(df=df,mu=mu,sigma=sigma)),
	marginsIdentical = TRUE)
M2 <- mvdc(claytonCopula(alpha,dim=2),c("tls","tls"),
	list( list(df=df,mu=mu,sigma=sigma), list(df=df,mu=mu,sigma=sigma)))
g <- fitMvdc(cbind(x,y),M1,c(3,1,1,10)) # constrains marginals to be identical
h <- fitMvdc(cbind(x,y),M2,c(3,1,1,3,1,1,10)) # estimates separate marginals

# Plot some true and estimated conditional quantile functions

z <- seq(min(y),max(y),len = 60)

for(i in 1:5){
	tau <- i/6
	uz <- ((tau^(-alpha/(1+alpha)) -1) * pt(z,3)^(-alpha) + 1)^(-1/alpha)
	yz <- qt(uz,3)
	lines(z,yz)
	b <- g@estimate
	uzhat <-((tau^(-b[4]/(1+b[4])) -1) * ptls(z,b[1],b[2],b[3])^(-b[4]) + 1)^(-1/b[4])
        yzhat <- qtls(uzhat,b[1],b[2],b[3])
	lines(z,yzhat,col="red")
	}

print(f)
print(g)
print(h)
