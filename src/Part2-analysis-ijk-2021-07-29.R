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
jags.data <- read_rds(here("data", "psophia_data_ijk.rds"))
#jags.data$Ind <- as.numeric(ifelse(jags.data$block == 2, 1, 0)) # recode blocks, B1=0; B2=1)
#attach(jags.data)
#with(jags.data, cor(cbind(block, elevation, range.elevation, distEdge, distWater, basalArea, recovery)))
# use elevation, distEdge, distWater, basalArea and recovery

jags.data$y <- as.numeric(jags.data$y)
#jags.data$y <- as.array(list(jags.data$y[,,1], jags.data$y[,,2], jags.data$y[,,3], jags.data$y[,,4]))
#array(y, dim=c(61,14,4))
#str(jags.data$y)
#class(jags.data$y)

#----- 4 - Dynamic occupancy model with covariates -----

# Specify model in BUGS language
sink(here("bin", "Dynocc_covariates.jags"))
cat("
model {

## Priors

# intercepts
alpha.psi ~ dnorm(0, 0.01)
alpha.p ~ dnorm(0, 0.01)

# psi coefficients
for (j in 1:5) {  # five predictors for psi: elevation, distEdge, distWater, basalArea and recovery
  beta.psi[j] ~ dnorm(0, 0.01)
}

# p coefficients
#beta.p ~ dnorm(0, 0.01)

# random year effects (for p)
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
  logit(psi[i,1]) <- alpha.psi + beta.psi[1]*elevation[i] + beta.psi[2]*distWater[i] + beta.psi[3]*basalArea[i] + beta.psi[4]*distEdge[i] + beta.psi[5]*recovery[i]
  z[i,1] ~ dbern(psi[i,1])

for (k in 2:nyear){
  muZ[i,k] <- z[i,k-1]*phi[k-1] + (1-z[i,k-1])*gamma[k-1]
  z[i,k] ~ dbern(muZ[i,k])
  } #k
 } #i

## Observation model
for (i in 1:nsite){
  for (j in 1:14){
    for (k in 1:nyear){
    logit(p[i,j,k]) <- alpha.p + eps[k]
    muy[i,j,k] <- z[i,k]*p[i,j,k] # can only be detected if z=1
    y[i,j,k] ~ dbern(muy[i,j,k])
    } #k
  } #j
 } #i

# Derived parameters: Population occupancy, growth rate and turnover
#psi[i,1] <- psi[i,1] # turned off because it was already defined in the likelihood
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


# Initial values
#zst <- apply(jags.data$y, c(1, 2), max)	# Observed occurrence as inits for z
#inits <- function(){ 
#  zst[is.na(zst)] <- 1 # NAs will result in error (node inconsistent with parents)
#  zst[zst > 1] <- 1
#  list(z = zst)
#}
#inits()

# Parameters monitored
params <- c("z", "psi", "phi", "gamma", "p",
            "alpha.psi", "alpha.p",
            "beta.psi", #"beta.p",
            "eps",
            "growthr", "turnover", "n.occ")#,
#"fit", "fit.new")


# MCMC settings
ni <- 150000
nt <- 100
nb <- 100000
#ni <- 25000
#nt <- 10
#nb <- 1000
nc <- 3

# remove unused variables from jags.data
jags.data$block <- NULL
jags.data$nblock <- NULL
#jags.data$elevation <- NULL
jags.data$range.elevation <- NULL
#jags.data$distWater <- NULL
#jags.data$distEdge <- NULL

# Call JAGS from R (BRT 3 min)
out <- jags(jags.data, inits=NULL, params, here("bin", "Dynocc_covariates.jags"), n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
saveRDS(out, here("results", "pobscura_ijk.rds"))


# coefficients 
coef.function <- function(x) {
  coefs <- data.frame(x$BUGSoutput$summary[c("beta.psi[1]", "beta.psi[2]", "beta.psi[3]", "beta.psi[4]", "beta.psi[5]"),])
  coefs <- tibble(predictor=c("elevation", "distWater", "basalArea", "distEdge", "recovery"),
                  coeff=row.names(coefs), mean=coefs$mean, lower=coefs$X2.5., upper=coefs$X97.5.,
                  Rhat=coefs$Rhat, n.eff=coefs$n.eff)
  .GlobalEnv$coefs <- coefs
  coefs
}
coef.function(out)

