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
jags.data <- read_rds(here("data", "psophia_data_ijk_40d.rds"))
#jags.data$Ind <- as.numeric(ifelse(jags.data$block == 2, 1, 0)) # recode blocks, B1=0; B2=1)
#attach(jags.data)
with(jags.data, cor(cbind(block, elevation, range.elevation, distEdge, distWater, basalArea, treeDensity, recovery, bouts)))
# use elevation, distEdge, distWater, basalArea and recovery
#jags.data$site <- seq(1:61)

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


# Initial values
#zst <- apply(jags.data$y, c(1, 2), max)	# Observed occurrence as inits for z
#inits <- function(){ 
#  zst[is.na(zst)] <- 1 # NAs will result in error (node inconsistent with parents)
#  zst[zst > 1] <- 1
#  list(z = zst)
#}

# Initial values: alternative
#inits <- function(){ list(z = apply(y, c(1, 3), max))}
#inits()

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
#inits()

# Parameters monitored
params <- c("z", "psi", "phi", "gamma", "p",
            "alpha.psi", "alpha.p",
            "beta.psi", #"beta.p",
            "eps",
            "growthr", "turnover", "n.occ")#,
#"fit", "fit.new")


# MCMC settings
ni <- 200000
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
jags.data$distWater <- NULL
#jags.data$distEdge <- NULL
jags.data$bouts <- NULL

# Call JAGS from R (BRT 3 min)
out <- jags(jags.data, inits, params, here("bin", "Dynocc_covariates_ijk.jags"), n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
saveRDS(out, here("results", "pobscura_mod_2021-09-20_ijk_40d.rds"))


# coefficients 
coef.function <- function(x) {
  coefs <- data.frame(x$BUGSoutput$summary[c("beta.psi[1]", "beta.psi[2]", "beta.psi[3]", "beta.psi[4]", "beta.psi[5]"),])
  coefs <- tibble(predictor=c("elevation", "distEdge","basalArea", "treeDensity", "recovery"),
                  coeff=row.names(coefs), mean=coefs$mean, lower=coefs$X2.5., upper=coefs$X97.5.,
                  Rhat=coefs$Rhat, n.eff=coefs$n.eff)
  .GlobalEnv$coefs <- coefs
  coefs
}
coef.function(out)


# Plot effects with uncertainty
predictor.effects <- function(x, original.predictor, coef) {
  
  #dev.off()
  
  ## Predict effect of logging on initial abundance with uncertainty
  mcmc.sample <- x$BUGSoutput$n.sims
  
  original.pred <- round(seq(min(original.predictor), max(original.predictor), length.out = 30),2)
  pred <- round((original.pred - mean(original.pred))/sd(original.pred), 2)
  #pred <- original.pred # it was already standardized
  psi.pred <- plogis( mean(x$BUGSoutput$sims.list$a0) + mean(x$BUGSoutput$sims.list$alpha.psi) +
                        mean(x$BUGSoutput$sims.list$beta[, coef]) * pred )
  
  array.psi.pred <- array(NA, dim = c(length(pred), mcmc.sample))
  for (i in 1:mcmc.sample){
    array.psi.pred[,i] <- plogis( x$BUGSoutput$sims.list$a0[i] + x$BUGSoutput$sims.list$alpha.psi[i] +
                                    x$BUGSoutput$sims.list$beta[i,coef] * pred )
  }
  
  # Plot for a subsample of MCMC draws
  sub.set <- sort(sample(1:mcmc.sample, size = 200))
  
  plot(original.pred, psi.pred, main = "", ylab = expression(psi), xlab = "", 
       ylim=c(0, 1), type = "l", lwd = 2, las=1, frame.plot = FALSE)
  for (i in sub.set){
    lines(original.pred, array.psi.pred[,i], type = "l", lwd = 1, col = "gray")
  }
  lines(original.pred, psi.pred, type = "l", lwd = 2, col = "blue")
  #mtext(expression(psi), side=2, line=3)
}
#predictor.effects(out, original.elevation, "a1")

predictor.effects(out, jags.data$elevation, 1)
#predictor.effects(out, original.distWater, 2)
predictor.effects(out2, jags.data$basalArea, 3)
predictor.effects(out2, jags.data$treeDensity, 4)
predictor.effects(out2, jags.data$recovery, 5)
