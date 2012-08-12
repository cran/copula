#### Testing fitCopula
######################################

require(copula)

source(system.file("Rsource", "tstFit-fn.R", package="copula", mustWork=TRUE))

## test code
system.time(
rr <- tstFit1cop(normalCopula(), tau.set=seq(0.2,0.8,by=0.2), n.set=c(25,50,100,200), N=200)
)
## ~ 400 seconds

d <- reshape.tstFit(rr)
plots.tstFit(d)# MM: desirable but ugly tick labelling

plots.tstFit(d, log=FALSE)


## t-copula instead of normal -- minimal set for testing here:
rt <- tstFit1cop(tCopula(), tau.set=c(.4, .8), n.set=c(10, 25), N=64)

##----- Now try fitting a tevCopula() :
## test code
system.time(
rtev <- tstFit1cop(tevCopula(), tau.set=seq(0.2,0.8,by=0.2), n.set=c(25,50,100,200), N=200)
)
## ~ ?? seconds
