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
#jags.data$Ind <- as.numeric(ifelse(jags.data$block == 2, 1, 0)) # recode blocks, B1=0; B2=1)
#attach(jags.data)

#----- 4 - Dynamic occupancy model with covariates -----

# Specify model in BUGS language
sink(here("bin", "Dynocc_covariates.jags"))
cat("
model {

## Priors

# intercepts
alpha.psi ~ dnorm(0, 0.01)
alpha.p ~ dnorm(0, 0.01)

# coefficients
for (j in 1:3) {
  beta.psi[j] ~ dnorm(0, 0.01)
  beta.p[j] ~ dnorm(0, 0.01)
}

for (k in 1:(nyear-1)) {
  gamma[k] ~ dunif(0, 1)
  phi[k] ~ dunif(0, 1)
}


## Likelihood

# ecological model
for (i in 1:nsite){
  logit(psi[i,1]) <- alpha.psi + beta.psi[1]*elevation[i] + beta.psi[2]*basalArea[i] + beta.psi[3]*recovery[i]
  z[i,1] ~ dbern(psi[i,1])

for (k in 2:nyear){
  muZ[i,k] <- z[i,k-1]*phi[k-1] + (1-z[i,k-1])*gamma[k-1]
  z[i,k] ~ dbern(muZ[i,k])
  } #k
 } #i

## Observation model
for (i in 1:nsite){
  for (k in 1:nyear){
    logit(p[i,k]) <- alpha.p + beta.p[1]*elevation[i] + beta.p[2]*basalArea[i] + beta.p[3]*recovery[i]
    muy[i,k] <- z[i,k]*p[i,k] # can only be detected if z=1
    y[i,k] ~ dbin(muy[i,k], nrep[i,k])
    } #k
  } #i

# Derived parameters: Sample and population occupancy, growth rate and turnover
#psi[i,1] <- psi[i,1] # turned off because it was already defined in the likelihood
for (i in 1:nsite){
  for (k in 2:nyear){
    psi[i,k] <- psi[i,k-1]*phi[k-1] + (1-psi[i,k-1])*gamma[k-1]
    growthr[i,k-1] <- mean(psi[i,k])/mean(psi[i,k-1]) # originally we had growthr[k]. JAGS seem to dislike vectoring going from 2..K.
    #growthr[k-1] <- mean(psi[,k])/mean(psi[,k-1])
    turnover[i,k-1] <- (1 - psi[i,k-1]) * gamma[k-1]/psi[i,k]
    #turnover[k-1] <- (1 - mean(psi[,k-1])) * gamma[k-1]/mean(psi[,k])
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
            "alpha.psi", "alpha.p",
            "beta.psi", "beta.p",
            #"eps.p",
            "growthr", "turnover")#,
#"fit", "fit.new")


# MCMC settings
ni <- 150000
nt <- 100
nb <- 75000
#ni <- 25000
#nt <- 10
#nb <- 1000
nc <- 3

# remove unused variables from jags.data
jags.data$block <- NULL
jags.data$nblock <- NULL
#jags.data$elevation <- NULL
jags.data$distWater <- NULL
jags.data$distEdge <- NULL

# Call JAGS from R (BRT 3 min)
out <- jags(jags.data, inits, params, here("bin", "Dynocc_covariates.jags"), n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
saveRDS(out, here("results", "pobscura_mod_3predictors_simplest.rds"))


# coefficients 
coef.function <- function(x) {
  coefs <- data.frame(x$BUGSoutput$summary[c("beta.psi[1]", "beta.psi[2]", "beta.psi[3]",
                                             "beta.p[1]", "beta.p[2]", "beta.p[3]"),])
  coefs <- tibble(predictor=rep(c("elevation", "basal.area", "recovery"), 2),
                  coeff=row.names(coefs), mean=coefs$mean, lower=coefs$X2.5., upper=coefs$X97.5.,
                  Rhat=coefs$Rhat, n.eff=coefs$n.eff)
  .GlobalEnv$coefs <- coefs
  coefs
}
coef.function(out)

