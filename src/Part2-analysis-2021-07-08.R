# Based on book BPA with WinBUGS
# Codes adapted from site BPA with JAGS

#----- 1 - Load libraries-----
#library(dplyr)
#library(lubridate)
library(here)
library(tidyverse)
library(R2jags)
library(rjags)
library(ggplot2)

#----- 2 - Source files-----
#source(here("bin", "ahumada_codes.R"))
#source(here("bin", "f-matrix-creator-experimental-probably-ok-but-need-check.R"))
#source(here("src", "fix_species_names.R")) # fix some names and remove "false" species
source(here("bin", "figures.R"))

#----- 3 - Read and prepare data -----
jags.data <- read_rds(here("data", "psophia_data.rds"))
jags.data$Ind <- as.numeric(ifelse(jags.data$block == 2, 1, 0)) # recode blocks, B1=0; B2=1)
#attach(jags.data)

#----- 4 - Dynamic occupancy model with covariates -----

# Specify model in BUGS language
sink(here("bin", "Dynocc_covariates.jags"))
cat("
model {

## Priors

# block intercepts
for(i in 1:4) { # one for each parameter psi, phi, gamma, p
  u.B1[i] ~ dnorm(mu.B1[i], tau.B1[i])
  mu.B1[i] <- log(B1.mean[i]) - log(1-B1.mean[i])
  B1.mean[i] ~ dunif(0,1)
  tau.B1[i] ~ dgamma(0.1,0.1)

  u.B2[i] ~ dnorm(mu.B2[i], tau.B2[i])
  mu.B2[i] <- log(B2.mean[i]) - log(1-B2.mean[i])
  B2.mean[i] ~ dunif(0,1)
  tau.B2[i] ~ dgamma(0.1,0.1)
  }

# coefficients
for (j in 1:3) {
  beta.psi[j] ~ dnorm(0, 0.01)
  beta.phi[j] ~ dnorm(0, 0.01)
  beta.gam[j] ~ dnorm(0, 0.01)
  }

# random year effects for p
for (k in 1:nyear){
  eps[k] ~ dnorm(0, 0.01)
  } #k


## Likelihood

# ecological model
for (i in 1:nsite){
  logit(psi[i,1]) <- u.B1[1]*(1-Ind[i]) + u.B2[1]*Ind[i] + beta.psi[1]*elevation[i] + beta.psi[2]*basalArea[i] + beta.psi[3]*recovery[i]
  z[i,1] ~ dbern(psi[i,1])

for (k in 2:nyear){
  logit(phi[i,k-1]) <- u.B1[2]*(1-Ind[i]) + u.B2[2]*Ind[i] + beta.phi[1]*elevation[i] + beta.phi[2]*basalArea[i] + beta.phi[3]*recovery[i]
  logit(gamma[i,k-1]) <- u.B1[3]*(1-Ind[i]) + u.B2[3]*Ind[i] + beta.gam[1]*elevation[i] + beta.gam[2]*basalArea[i] + beta.gam[3]*recovery[i]

  muZ[i,k] <- z[i,k-1]*phi[i,k-1] + (1-z[i,k-1])*gamma[i,k-1]
  z[i,k] ~ dbern(muZ[i,k])
  } #k
 } #i

## Observation model
for (i in 1:nsite){
  for (k in 1:nyear){
    logit(p[i,k]) <- u.B1[4]*(1-Ind[i]) + u.B2[4]*Ind[i] + eps[k]
    muy[i,k] <- z[i,k]*p[i,k] # can only be detected if z=1
    y[i,k] ~ dbin(muy[i,k], nrep[i,k])
    } #k
  } #i

# Derived parameters: Sample and population occupancy, growth rate and turnover
#psi[i,1] <- psi[i,1]
for (i in 1:nsite){
  for (k in 2:nyear){
    psi[i,k] <- psi[i,k-1]*phi[i,k-1] + (1-psi[i,k-1])*gamma[i,k-1]
    growthr[i,k-1] <- psi[i,k]/psi[i,k-1] # originally we had growthr[k]. JAGS seem to dislike vectoring going from 2..K.
    turnover[i,k-1] <- (1 - psi[i,k-1]) * gamma[i,k-1]/psi[i,k]
    } # k
  } #i
}
",fill = TRUE)
sink()


# Initial values
zst <- apply(jags.data$y, c(1, 2), max)	# Observed occurrence as inits for z
inits <- function(){ 
  zst[is.na(zst)] <- 1 # NAs will result in error (node inconsistent with parents)
  zst[zst > 1] <- 1
  list(z = zst)
}
#inits()

# Parameters monitored
params <- c("z", "psi", "phi", "gamma", "p",
            "u.B1", "u.B2",
            #"alpha.psi", "alpha.phi", "alpha.gam", "alpha.p",
            "beta.psi", "beta.phi", "beta.gam", "beta.p",
            #"w.psi", "w.phi", "w.gam", #"w.p",
            "eps",
            "growthr", "turnover")#,
#"fit", "fit.new")


# MCMC settings
#ni <- 100000
#nt <- 100
#nb <- 50000
ni <- 25000
nt <- 10
nb <- 1000
nc <- 3

# remove unused variables from jags.data
jags.data$block <- NULL
jags.data$nblock <- NULL
#jags.data$elevation <- NULL
jags.data$distWater <- NULL
jags.data$distEdge <- NULL

# Call JAGS from R (BRT 3 min)
out <- jags(jags.data, inits, params, here("bin", "Dynocc_covariates.jags"), n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
saveRDS(out, here("results", "pobscura_mod_3predictors.rds"))
