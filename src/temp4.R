
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
attach(jags.data)

#----- 4 - Dynamic occupancy model with covariates -----

# Specify model in BUGS language
sink(here("bin", "Dynocc_covariates.jags"))
cat("
model {

  ## Priors
  
  # random site effects
  for (i in 1:nsite){
      alpha.psi[i] ~ dnorm(mu.psi, tau.psi)
      alpha.phi[i] ~ dnorm(mu.phi, tau.phi)
      alpha.gam[i] ~ dnorm(mu.gam, tau.gam)
      alpha.p[i] ~ dnorm(mu.p, tau.p)
  } #i
  
  mu.psi ~ dnorm(0, 0.01)  # Hyperparameter 1
  tau.psi <- 1 / (sd.psi*sd.psi)  # Hyperparameter 2
  sd.psi ~ dunif(0, 2)
  
  mu.phi ~ dnorm(0, 0.01)  # Hyperparameter 1
  tau.phi <- 1 / (sd.phi*sd.phi)  # Hyperparameter 2
  sd.phi ~ dunif(0, 2)
  
  mu.gam ~ dnorm(0, 0.01)  # Hyperparameter 1
  tau.gam <- 1 / (sd.gam*sd.gam)  # Hyperparameter 2
  sd.gam ~ dunif(0, 2)

  mu.p ~ dnorm(0, 0.01)  # Hyperparameter 1
  tau.p <- 1 / (sd.p*sd.p)  # Hyperparameter 2
  sd.p ~ dunif(0, 2)
  
  
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
    logit(psi[i,1]) <-  alpha.psi[i] + beta.psi[1]*distEdge[i] + beta.psi[2]*basalArea[i] + beta.psi[3]*recovery[i]
    z[i,1] ~ dbern(psi[i,1])
    
    for (k in 2:nyear){
      logit(phi[i,k-1]) <- alpha.phi[i] + beta.phi[1]*distEdge[i] + beta.phi[2]*basalArea[i] + beta.phi[3]*recovery[i]
      logit(gamma[i,k-1]) <- alpha.gam[i] + beta.gam[1]*distEdge[i] + beta.gam[2]*basalArea[i] + beta.gam[3]*recovery[i]
      
      muZ[i,k] <- z[i,k-1]*phi[i,k-1] + (1-z[i,k-1])*gamma[i,k-1]
      z[i,k] ~ dbern(muZ[i,k])
      } #k
   } #i
   
   ## Observation model
   for (i in 1:nsite){
       for (k in 1:nyear){
         logit(p[i,k]) <-  alpha.p[i] + eps[k]
         muy[i,k] <- z[i,k]*p[i,k]  # can only be detected if z=1
         y[i,k] ~ dbin(muy[i,k], nrep[i,k])
         
   # Assess model fit using Chi-squared discrepancy
   # Compute fit statistic E for observed data
   eval[i,k] <- p[i,k] * z[i,k]   	# Expected values
   E[i,k] <- pow((y[i,k] - eval[i,k]),2) / (eval[i,k] + 0.5)
   # Generate replicate data and compute fit stats for them
   y.new[i,k] ~ dbin(p[i,k], nrep[i,k])
   E.new[i,k] <- pow((y.new[i,k] - eval[i,k]),2) / (eval[i,k] + 0.5)
         
      } #k
    } #i

   
   # Derived parameters: Sample and population occupancy, growth rate and turnover
   #  psi[i,1] <- psi[i,1]
   for (i in 1:nsite){
     for (k in 2:nyear){
       psi[i,k] <- psi[i,k-1]*phi[i,k-1] + (1-psi[i,k-1])*gamma[i,k-1]
       growthr[i,k-1] <- psi[i,k]/psi[i,k-1]  # originally we had growthr[k]. JAGS seem to dislike vectoring going from 2..K.
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
            "alpha.psi", "alpha.p", "alpha.phi", "alpha.gam",
            "beta.psi", "beta.phi", "beta.gam", "beta.p",
            "growthr", "turnover",
            "fit", "fit.new")


# MCMC settings
ni <- 100000
nt <- 100
nb <- 50000
#ni <- 5000
#nt <- 10
nb <- 2500

nc <- 3

# Call JAGS from R (BRT 3 min)
out <- jags(jags.data, inits, params, here("bin", "Dynocc_covariates.jags"), n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)


#saveRDS(out, here("results", "pobscura_model.rds"))


# check convergence
hist(out$BUGSoutput$summary[,"Rhat"], nclass=8, main="Rhat", xlab="", las=1)
summary(out$BUGSoutput$summary[,"Rhat"])

# which parameters did not converged well?
out$BUGSoutput$summary[which(out$BUGSoutput$summary[,"Rhat"] > 1.1),]

# Evaluate fit
plot(out$BUGSoutput$sims.list$fit, out$BUGSoutput$sims.list$fit.new, main = "", xlab = "Discrepancy actual data", ylab = "Discrepancy replicate data", frame.plot = FALSE)
abline(0, 1, lwd = 2, col = "black")
mean(out0$BUGSoutput$sims.list$fit.new > out0$BUGSoutput$sims.list$fit)
mean(out0$BUGSoutput$mean$fit) / mean(out0$BUGSoutput$mean$fit.new)

# bayesian p-values is the proportion of points above the 1:1 line of equality
# bayesian p-values close to zero or one are suspicious, we've got 0.58 which is a good fit!
# ADD CODE

# coefficients 
coef.function <- function(x) {
  coefs <- data.frame(x$BUGSoutput$summary[c("beta.psi[1]", "beta.psi[2]", "beta.psi[3]",
                                             "beta.phi[1]", "beta.phi[2]", "beta.phi[3]",
                                             "beta.gam[1]", "beta.gam[2]", "beta.gam[3]"),])
                                             #"beta.p[1]", "beta.p[2]", "beta.p[3]", "beta.p[4]"),])
  coefs <- tibble(predictor=c("edges.psi", "basal.area.psi", "recovery.psi",
                              "edges.phi", "basal.area.phi", "recovery.phi",
                              "edges.gam", "basal.area.gam", "recovery.gam"),
                              #"block.p", "edges.p", "basal.area.p", "recovery.p"),
                  coeff=row.names(coefs), mean=coefs$mean, lower=coefs$X2.5., upper=coefs$X97.5.,
                  Rhat=coefs$Rhat, n.eff=coefs$n.eff)
  .GlobalEnv$coefs <- coefs
  coefs
}
coef.function(out)

predictor.effects.psi(out, log(original.distEdge), 1)
mtext("Distance to edge (m)", side=1, line=3)

predictor.effects.psi(out, original.basalArea, 2)
mtext("Basal area (m²/ha)", side=1, line=3)

predictor.effects.psi(out, original.recovery, 3)
mtext("Recovery time (years)", side=1, line=3)

predictor.effects.gam(out, original.basalArea, 2)
mtext("Basal area (m²/ha)", side=1, line=3)


