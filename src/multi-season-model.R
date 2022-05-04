# Run multi-season occupancy model in JAGS

# Some chunks were adapted from:
# https://www.vogelwarte.ch/de/projekte/publikationen/bpa/code-for-running-bpa-using-jags
# J. Andrew Royle and Marc Kéry. 2007. A Bayesian state-space formulation of dynamic occupancy models. Ecology 88:1813–1823 (Supplement)

#----- 1 - Load libraries-----
library(here)
library(tidyverse)
library(R2jags)
library(ggplot2)

#----- 2 - Source files-----
# (no source files needed for this part)


#----- 3 - Read and prepare data -----
jags.data <- read_rds(here("data", "psophia_data_ijk.rds"))

# check correlation between predictors
with(jags.data, cor(cbind(block, elevation, range.elevation, distEdge, distWater, basalArea, treeDensity, recovery, bouts)))


#----- 4 - Dynamic occupancy model with covariates -----
# Specify model in BUGS language
sink(here("bin", "Dynocc_covariates_ijk.jags"))
cat("
model {

## Priors

# random site effects (intercepts)
for (i in 1:nsite){
   alpha.psi[i] ~ dnorm(mu.psi, tau.alpha.psi)
   alpha.p[i] ~ dnorm(mu.p, tau.alpha.p)
   }
mu.psi ~ dnorm(0, 0.01)  # hyperparameter
tau.alpha.psi <- 1 / (sd.alpha.psi*sd.alpha.psi)  # hyperparameter
sd.alpha.psi ~ dunif(0, 2)
mu.p ~ dnorm(0, 0.01)  # hyperparameter
tau.alpha.p <- 1 / (sd.alpha.p*sd.alpha.p)  # hyperparameter
sd.alpha.p ~ dunif(0, 2)

# psi coefficients
for (j in 1:5) {  # five predictors for psi: elevation, basalArea, treeDens, distEdge and recovery
  beta.psi[j] ~ dnorm(0, 0.01)
}


# random year effects (for p only)
for (k in 1:nyear){
   eps[k] ~ dnorm(0, tau.year)
} # k
tau.year <- 1 / (sd.year*sd.year)
sd.year ~ dunif(0, 1)

# gamma and phi varying by year
for (k in 1:(nyear-1)) {
  gamma[k] ~ dunif(0, 1)
  phi[k] ~ dunif(0, 1)
}


## Likelihood

# ecological model
for (i in 1:nsite){
  logit(psi1[i]) <- alpha.psi[i] + beta.psi[1]*elevation[i] + beta.psi[2]*distEdge[i] + beta.psi[3]*basalArea[i] + beta.psi[4]*treeDensity[i] + beta.psi[5]*recovery[i]
  z[i,1] ~ dbern(psi1[i])

for (k in 2:nyear){
  muZ[i,k] <- z[i,k-1]*phi[k-1] + (1-z[i,k-1])*gamma[k-1]
  z[i,k] ~ dbern(muZ[i,k])
  } #k
 } #i

## Observation model
for (i in 1:nsite){
  for (j in 1:nocc){
    for (k in 1:nyear){
    
      logit(p[i,j,k]) <- alpha.p[i] + eps[k]
      muy[i,j,k] <- z[i,k]*p[i,j,k] # can only be detected if z=1
      y[i,j,k] ~ dbern(muy[i,j,k])

    }#k
  } #j
} #i

# Derived parameters: Population occupancy, growth rate and turnover
for (i in 1:nsite){
  psi[i,1] <- psi1[i]
 } # i
n.occ[1] <- sum(z[1:nsite,1])
for (i in 1:nsite){
  for (k in 2:nyear){
    psi[i,k] <- psi[i,k-1]*phi[k-1] + (1-psi[i,k-1])*gamma[k-1]
    } #k
  } #i
  
for (k in 2:nyear){
  n.occ[k] <- sum(z[,k])
  growthr[k-1] <- mean(psi[,k])/mean(psi[,k-1]) # originally we had growthr[k]. JAGS seem to dislike vectoring going from 2..K.
  turnover[k-1] <- (1 - mean(psi[,k-1])) * gamma[k-1]/mean(psi[,k])
 } #k

}
",fill = TRUE)
sink()


# Initial values: another alternative
inits <- function(){
  inits_a <- rowSums(jags.data$y[,,1], na.rm=TRUE)
  inits_a[inits_a > 0] <- 1
  inits_b <- rowSums(jags.data$y[,,2], na.rm=TRUE)
  inits_b[inits_b > 0] <- 1
  inits_c <- rowSums(jags.data$y[,,3], na.rm=TRUE)
  inits_c[inits_c > 0] <- 1
  inits_d <- rowSums(jags.data$y[,,4], na.rm=TRUE)
  inits_d[inits_d > 0] <- 1
  inits_e <- rowSums(jags.data$y[,,5], na.rm=TRUE)
  inits_e[inits_e > 0] <- 1
  list(z=as.matrix(cbind(inits_a, inits_b, inits_c, inits_d, inits_e)))
  }

# Parameters monitored
params <- c("z", "psi", "phi", "gamma", "p",
            "alpha.psi", "alpha.p",
            "beta.psi", #"beta.p",
            "eps",
            "growthr", "turnover", "n.occ")


# MCMC settings
ni <- 250000
nt <- 150
nb <- 150000
nc <- 3

# remove unused variables from jags.data
jags.data$block <- NULL
jags.data$nblock <- NULL
jags.data$range.elevation <- NULL
jags.data$distWater <- NULL
jags.data$bouts <- NULL

# Call JAGS from R
out <- jags(jags.data, inits, params, here("bin", "Dynocc_covariates_ijk.jags"), n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

# save model results
saveRDS(out, here("results", "pobscura_model.rds"))

