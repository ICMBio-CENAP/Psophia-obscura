
library(R2jags)

data <- source(here("bin", "img.data.R"))

## supplied data have logit(psi) = -1 + 2.5*x[i]
## where x[i] is the number of neighboring cells occupied
nG <- data$value$nG
griddim <- data$value$griddim
z <- data$value$z
numnn <- data$value$numnn
N <- data$value$N

Zmat <- matrix(NA, nrow=griddim, ncol=griddim)
Zmat[1:nG]<-z

par(cex.axis=1.8,cex.lab=2.0,mar= c(5, 4, 4, 2)*1.2 + 0.1)
image(1:40,1:40,Zmat,col=rev(terrain.colors(10)),xlab="Easting",ylab="Northing")

sink(here("bin", "model.txt"))
cat("
model{

alpha ~ dnorm(0,.01)
beta  ~ dnorm(0,.01)

for(i in 1:nG){
x[i,1]<-0
for(j in 1:numnn[i]){
x[i,j+1]<-x[i,j]+z[N[i,j]]
}
logit(psi[i])<- alpha + beta*(x[i,numnn[i]+1]/numnn[i])
z[i]~dbern(psi[i])
}
}
",fill=TRUE)
sink()

data <- list ( "z","numnn","N","nG")
inits <- function()
  list ( alpha=rnorm(1),beta=rnorm(1) )

parameters <- c("alpha","beta")

ni=2000
nb=1000
nt=1
nc=3

fit <- jags(data, inits, parameters, here("bin", "model.txt"), n.thin=nt,n.chains=nc, n.burnin=nb,n.iter=ni)

fit
