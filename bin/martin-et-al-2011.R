
# Martin et al 2011. Accounting for non-independent detection when estimating abundance 
# of organisms with a Bayesian approach. MEE

#Appendix S1

############################################################
#Simulation of data; data follow a binomial distribution;
#Note this script was modified from Kery & Schaub (2012)
#this script requires package "R2WinBUGS"
############################################################
n.site <- 200# 200 sites visited
lam <- 4.6 #parameter lambda
N=rpois(n = n.site, lambda = lam)  #true abundance at each site
R <- n.site #Number of sites
T <- 3			# Number of replicate counts at each site
y <- array(dim = c(R, T)) #create array
# Simulate count results of first through last surveys, here p = 0.5
for(i in 1:T){
  y[,i] <- rbinom(n = n.site, size = N, prob = 0.5)
}
C<-c(y) #count in column format

site = 1:R				# 'Short' version of site covariate
site.p <- rep(site, T)		# 'Long' version of site covariate

############################################################
#WinBUGS codes: standard binomial mixture model
############################################################
library("here")
library("R2jags")			# Load R2WinBUGS package
sink(here("bin", "Model.txt"))
cat("
model {
# Priors
p0~dunif(0,1)      #detection probability parameter
lam~dgamma(.01,.01)   #lambda parameter

# Likelihood

# Biological model for true abundance
 for (i in 1:R) {			# Loops over R sites
   N[i] ~ dpois(lambda[i]) #Abundance at each site follows a Poisson distribution
   lambda[i] <- lam
  }

# Observation model for replicated counts 
 for (i in 1:n) {			# Loops over all n observations
   C[i] ~ dbin(p[i], N[site.p[i]]) #the count data follows a binomial distribution
   p[i] <-p0
 }

# Derived quantities
 totalN <- sum(N[])	# Estimate total population size across all sites
}
",fill=TRUE)
sink()

# Package data for WinBUGS
R = dim(y)[1]   # number of sites
n = dim(y)[1] * dim(y)[2]#number of observations (sites*surveys)
win.data <- list(R = R, n = n, C = C, site.p = site.p)

# Initial values
Nst <- apply(y, 1, max) + 1
inits <- function(){list(N = Nst,lam=runif(1,1,8),  p0=runif(1))}

# parameters to monitor
parameters <- c( "totalN", "p0", "lam")

# MCMC settings
nc <- 3 #number of chains
nb <- 4000 # “burn in”
ni <- 14000# “number of iterations”
nt <- 1 # “thinning”

# run model
out <- jags (win.data, inits, parameters, here("bin", "Model.txt"), n.chains=nc, n.iter=ni, n.burn = nb, n.thin=nt)

bias= mean(out$BUGSoutput$sims.list$totalN)-sum(N)
relative_bias= bias/sum(N)#This is the relative bias for 1 simulated data set


#Appendix S2

############################################################
#Simulation of data; data follow a beta-binomial distribution;
#this script requires package "VGAM" and "R2WinBUGS"
############################################################
library(VGAM) #load package "VGAM", necessary to be able to use rbetabin()
n.site <- 200# 200 sites visited
lam <- 4.6 #parameter lambda
N=rpois(n = n.site, lambda = lam)  #true abundance at each site

R <- n.site #number of sites
T <- 3			# Number of replicate counts at each site
y <- array(dim = c(R, T))  #array of dimension RxT

for(r in 1:R){        #loop that fills the y array with counts that follow a beta-binomial distribution
  if (N[r]>0) {y[r,]<- rbetabin(n=T,size=N[r],prob=0.5,rho=0.3)}  #with parameter p =0.5 and rho =0.3
  else
  {y[r,] <- 0}
}


C<-c(y) #count in column format
C <- as.numeric(C)

site = 1:R				# 'Short' version of site covariate
site.p <- rep(site, T)		# 'Long' version of site covariate



############################################################
#WinBUGS codes: beta-binomial mixture model
############################################################
library(R2jags)			# Load R2WinBUGS package
sink(here("bin", "Model.txt"))
cat("
model {

# Priors
lam~dgamma(.01,.01)
alpha.p~dgamma(.01,.01)
beta.p~dgamma(.01,.01)

# Likelihood
# Biological model for true abundance
for (i in 1:R) { # Loops over R sites
   N[i] ~ dpois(lambda[i])  #Abundance at each site follows a Poisson distribution
   lambda[i] <- lam
}


# Observation model for replicated counts
 for (i in 1:n) {			# Loops over all n observations
C[i] ~ dbin(p[i], N[site.p[i]]) #counts follow a beta-binomial distribution
        p[i]~dbeta(alpha.p,beta.p)  #detection probability p follows a beta distribution
 }


# Derived quantities
totalN <- sum(N[])	# Estimate total population size across all sites
p.derived<-alpha.p/(alpha.p+beta.p) #derived detection probability
rho.derived<-1/(alpha.p+beta.p+1)   #derived correlation coefficient
 }

",fill=TRUE)
sink()


# Package data for WinBUGS
R = dim(y)[1]   # number of sites:
n = dim(y)[1] * dim(y)[2]#number of observations (sites*surveys)
win.data <- list(R = R, n = n,C = C, site.p = site.p)

# Initial values
Nst <- apply(y, 1, max,na.rm=T) + 1
inits <- function(){list(N = Nst, lam=runif(1,1,8),alpha.p=runif(1,0.5,1.5), beta.p=runif(1,0.5,1.5))}

# Parameters to monitor
parameters <- c("lam","totalN","alpha.p","beta.p","p.derived","rho.derived")

# MCMC settings
nc <- 3
nb <- 4000
ni <- 14000
nt <- 1

# Output
out <- jags (win.data, inits, parameters, here("bin", "Model.txt"), n.chains=nc, n.iter=ni, n.burn = nb, n.thin=nt)




#Appendix S3. Estimation of the probability of a site being occupied conditional on manatees not being detected at the site (ψcondl).
#The probability that a site was occupied by a manatee given that no manatees were detected after K surveys( ψcondl) was (MacKenzie et al. 2006):
#S1
#ψ and p were estimated with a Bayesian formulation of site occupancy models (MacKenzie et al. 2006, p. 127).
#According to this model the estimate of ψ was 0.25 and the estimate of p was 0.76. Therefore, based on eqn S1, the estimate of ψcondl was 0.0046.


#Appendix S4
############################################################
#Beta-binomial mixture model, with abundance assumed to 
#follow a Poisson distribution.
#To account for the variation in size of the sampling units,
#we used the log of the size of the sampling unit as an 
#additive offset.
############################################################
library(R2jags)
sink(here("bin", "Model.txt"))
cat("
model {

# Priors
alpha.lam~dgamma(.01,.01)
alpha.p~dgamma(.01,.01)
beta.p~dgamma(.01,.01)

# Likelihood
# Ecological model for true abundance
for(i in 1:R){

 N[i] ~ dpois(lambda[i])
   log(lambda[i]) <-alpha.lam + 1* logA[i] # we used the log of the size of the sampling unit (logA) as an additive offset
}
      # Observation model for replicated counts
      for (i in 1:n) {                   # Loop over temporal replications of counts
         C[i] ~ dbin(p[i], N[site.p[i]])   # Counts follow a binomial distribution ; count data should be formatted as shown in ESM1
         p[i]~dbeta(alpha.p,beta.p)        #Detection follows a beta distribution


      } # ends i loop


# Derived and other quantities
   totalN <- sum(N[])	# Estimate of total abundance across all sites
   mean.abundance <- mean(lambda[]) #mean expected abundance per site
   p.derived<-alpha.p/(alpha.p+beta.p) #derived detection probability
   rho.derived<-1/(alpha.p+beta.p+1)  #derived correlation coefficient
}
",fill = TRUE)
sink()


# Bundle data
R = nrow(y)  #number of sites
T = ncol(y)  #number of temporal replicates
n = dim(y)[1] * dim(y)[2]#number of observations (sites*surveys)

win.data <- list(C = C, n = n,R = R, site.p = site.p,logA=log(Area)) #Area: size of the units


# Initial values
Nst <- apply(y, 1, max) + 1
Nst[is.na(Nst)] <- 1
inits <- function(){list(N = Nst, alpha.lam = runif(1, 1, 8),alpha.p=runif(1,0.5,1), beta.p=runif(1,0.5,1))}

# Define parameters to be monitored
params <- c("totalN","mean.abundance", "p.derived", "rho.derived", "alpha.p","beta.p")

# MCMC settings
ni <- 14000
nt <- 1
nb <- 4000
nc <- 3

out <- jags(win.data, inits, params, here("bin", "Model.txt"), n.chains = nc,
            n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())


# Inspect output
print(out, dig = 3)




#Appendix S5
############################################################
#Beta-binomial mixture model, with abundance assumed to 
#follow a zero-truncated Poisson distribution.
#the script also presents a goodness of fit procedure.
############################################################
library(R2jags)
sink(here("bin", "Model.txt"))
cat("
model {

# Priors
lam~dgamma(.01,.01)
alpha.p~dgamma(.01,.01)
beta.p~dgamma(.01,.01)

# Likelihood
#the next four lines are used to model N as a zero-truncated distribution
probs[R+1]<- 1-sum(probs[1:R])
for(i in 1:R){
probs[i]<- exp(-lam)*(pow(lam,x[i]))/(exp(logfact(x[i])) * (1-exp(-lam)))
 N[i] ~ dcat(probs[])}

      # Observation model for replicated counts
      for (i in 1:n) {                   
         C[i] ~ dbin(p[i], N[site.p[i]])
         p[i]~dbeta(alpha.p,beta.p)

         # Assess model fit using Chi-squared discrepancy
         # Compute fit statistic for observed data
         eval[i] <- p[i]*N[site.p[i]]
         E[i] <- pow((C[i] - eval[i]),2) / (eval[i] + 0.5)
         # Generate replicate data and compute fit stats
         C.new[i] ~ dbin(p[i], N[site.p[i]])
         E.new[i] <- pow((C.new[i] - eval[i]),2) / (eval[i] + 0.5)

      } # ends i loop


# Derived and other quantities
   totalN <- sum(N[])	# Estimate abundance across all sites
   mean.abundance <- lam #mean expected abundance per plot
   p.derived<-alpha.p/(alpha.p+beta.p) #derived detection probability
   rho.derived<-1/(alpha.p+beta.p+1)  #correlation coefficient


fit <- sum(E[])
fit.new <- sum(E.new[])

}
",fill = TRUE)
sink()


# Bundle data
R = nrow(y)
T = ncol(y)
n = dim(y)[1] * dim(y)[2]#number of observations (sites*surveys)

win.data <- list(C = C, n=n,R = R, site.p = site.p,x=1:R)

# Initial values
Nst <- apply(y, 1, max) + 1
Nst[is.na(Nst)] <- 1
inits <- function(){list(N = Nst, lam = runif(1, 1, 7),alpha.p=runif(1,0.5,1), beta.p=runif(1,0.5,1))}

# Define parameters to be monitored
params <- c("totalN", "mean.abundance", "lam", "p.derived", "rho.derived", "fit", "fit.new","alpha.p","beta.p")

# MCMC settings
ni <- 14000
nt <- 1
nb <- 4000
nc <- 3

out <- jags(win.data, inits, params, here("bin", "Model.txt"), n.chains = nc,
            n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())


# Inspect output
print(out, dig = 3)

mean(out$BUGSoutput$sims.list$fit.new > out$BUGSoutput$sims.list$fit) #Bayesian P-value (reported in the text)

#Appendix S6. Graphical representation of relative bias for each of the 16 estimators considered in Table 1. The dots correspond to the estimates of relative bias and the error bars correspond to the 95%CI.
#Appendix S7. Posterior distribution for λ, p and ρ obtained by fitting the beta-binomial mixture models (with zero-truncated distribution for λ) to aerial survey data of manatees.