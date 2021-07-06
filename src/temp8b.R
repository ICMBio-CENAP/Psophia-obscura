
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
jags.data$block <- as.numeric(ifelse(jags.data$block == 2, 1, 0)) # recode block to zeroes and ones
#attach(jags.data)

#----- 4 - Dynamic occupancy model with covariates -----

# Specify model in BUGS language
sink(here("bin", "Dynocc_covariates.jags"))
cat("
model {

  ## Priors
  
  # intercepts
  alpha.psi ~ dnorm(0, 0.01)
  alpha.phi ~ dnorm(0, 0.01)
  alpha.gam ~ dnorm(0, 0.01)
  alpha.p ~ dnorm(0, 0.01)
  
  # coefficients
  for (j in 1:4) {
      beta.psi[j] ~ dnorm(0, 0.01)
      beta.phi[j] ~ dnorm(0, 0.01)
      beta.gam[j] ~ dnorm(0, 0.01)
      beta.p[j] ~ dnorm(0, 0.01)
  }
  
  ## Likelihood
  
  # ecological model
  for (i in 1:nsite){
    logit(psi[i,1]) <-  alpha.psi + beta.psi[1]*block[i] + beta.psi[2]*distEdge[i] + beta.psi[3]*basalArea[i] + beta.psi[4]*recovery[i]
    z[i,1] ~ dbern(psi[i,1])
    
    for (k in 2:nyear){
      logit(phi[i,k-1]) <- alpha.phi + beta.phi[1]*block[i] + beta.phi[2]*distEdge[i] + beta.phi[3]*basalArea[i] + beta.phi[4]*recovery[i]
      logit(gamma[i,k-1]) <- alpha.gam + beta.gam[1]*block[i] + beta.gam[2]*distEdge[i] + beta.gam[3]*basalArea[i] + beta.gam[4]*recovery[i]
      
      muZ[i,k] <- z[i,k-1]*phi[i,k-1] + (1-z[i,k-1])*gamma[i,k-1]
      z[i,k] ~ dbern(muZ[i,k])
      } #k
   } #i
   
   ## Observation model
   for (i in 1:nsite){
       for (k in 1:nyear){
         logit(p[i,k]) <-  alpha.p + beta.p[1]*block[i] + beta.p[2]*distEdge[i] + beta.p[3]*basalArea[i] + beta.p[4]*recovery[i]
         muy[i,k] <- z[i,k]*p[i,k]  # can only be detected if z=1
         y[i,k] ~ dbin(muy[i,k], nrep[i,k])
         
   # Assess model fit using Chi-squared discrepancy
   # Compute fit statistic E for observed data
   #eval[i,k] <- p[i,k]*z[i,k]  # Expected values
   #E[i,k] <- pow((y[i,k] - eval[i,k]),2) / (eval[i,k] + 0.5)
   # Generate replicate data and compute fit stats for them
   #y.new[i,k] ~ dbin(muy[i,k], nrep[i,k])
   #E.new[i,k] <- pow((y.new[i,k] - eval[i,k]),2) / (eval[i,k] + 0.5)
         
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
    #fit <- sum(E[,])
    #fit.new <- sum(E.new[,])
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
            "alpha.psi", "alpha.phi", "alpha.gam", "alpha.p",
            "beta.psi", "beta.phi", "beta.gam", "beta.p",
            #"w.psi", "w.phi", "w.gam", "w.p",
            #"betaT.psi", "betaT.phi", "betaT.gam", "betaT.p",
            #"eps",
            "growthr", "turnover")
#"fit", "fit.new")


# MCMC settings
ni <- 250000
nt <- 100
nb <- 100000
#ni <- 5000
#nt <- 10
#nb <- 2500
nc <- 3

# remove unused variables from jags.data
#jags.data$block <- NULL
jags.data$nblock <- NULL
#jags.data$elevation <- NULL
jags.data$distWater <- NULL

# Call JAGS from R (BRT 3 min)
out <- jags(jags.data, inits, params, here("bin", "Dynocc_covariates.jags"), n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
saveRDS(out, here("results", "pobscura_mod.rds"))

#out <- read_rds(here("results", "pobscura_mod_sel.rds"))


# coefficients 
coef.function <- function(x) {
  coefs <- data.frame(x$BUGSoutput$summary[c("beta.psi[1]", "beta.psi[2]", "beta.psi[3]", "beta.psi[4]",
                                             "beta.phi[1]", "beta.phi[2]", "beta.phi[3]","beta.phi[4]",
                                             "beta.gam[1]", "beta.gam[2]", "beta.gam[3]","beta.gam[4]",
                                             "beta.p[1]", "beta.p[2]", "beta.p[3]", "beta.p[4]"),])
  coefs <- tibble(predictor=c("block", "edges", "basal.area", "recovery",
                              "block", "edges", "basal.area", "recovery",
                              "block", "edges", "basal.area", "recovery",
                              "block", "edges", "basal.area", "recovery"),
                  coeff=row.names(coefs), mean=coefs$mean, lower=coefs$X2.5., upper=coefs$X97.5.,
                  Rhat=coefs$Rhat, n.eff=coefs$n.eff)
  .GlobalEnv$coefs <- coefs
  coefs
}
coef.function(out)

# check model averaged distribution of a coefficient
par(mfrow=c(1,2))
hist(x$BUGSoutput$sims.list$beta.psi[,4], nclass=20, xlab="beta", main="")
hist(x$BUGSoutput$sims.list$beta.psi[,5], nclass=20, xlab="beta", main="")
dev.off()


# inclusion variables 
inc.function <- function(x) {
  coefs <- data.frame(x$BUGSoutput$summary[c("w.psi[1]", "w.psi[2]", "w.psi[3]", "w.psi[4]",
                                             "w.phi[1]", "w.phi[2]", "w.phi[3]","w.phi[4]",
                                             "w.gam[1]", "w.gam[2]", "w.gam[3]","w.gam[4]",
                                             "w.p[1]", "w.p[2]", "w.p[3]", "w.p[4]"),])
  coefs <- tibble(predictor=c("block", "edges", "basal.area", "recovery",
                              "block", "edges", "basal.area", "recovery",
                              "block", "edges", "basal.area", "recovery",
                              "block", "edges", "basal.area", "recovery"),
                  coeff=row.names(coefs), mean=coefs$mean, lower=coefs$X2.5., upper=coefs$X97.5.,
                  Rhat=coefs$Rhat, n.eff=coefs$n.eff)
  .GlobalEnv$coefs <- coefs
  coefs
}
inc.function(out)


# read pobscura file to get original predictor values
pobscura <- read_rds(here("data", "pobscura.rds"))
#SiteCovs <- pobscura[,42:50] 
#cor(SiteCovs)

coef.function(out)
dev.off()
predictor.effects.psi(out, pobscura$basal.area, 3)
mtext("Basal area (m²/ha)", side=1, line=3)

predictor.effects.psi(out, pobscura$recovery, 4)
mtext("Recovery time (years)", side=1, line=3)

predictor.effects.p(out, pobscura$recovery, 4)
mtext("Basal area (m²/ha)", side=1, line=3)

predictor.effects.p(out, pobscura$basal.area, 4)
mtext("Basal area (m²/ha)", side=1, line=3)

par(mfrow = c(2, 1))
hist(plogis(out$BUGSoutput$sims.list$alpha.phi), nclass = 40, col = "gray", main = "Block 1", xlab = "Survival probability", xlim = c(0, 1))
hist(plogis(out$BUGSoutput$sims.list$alpha.phi + out$BUGSoutput$sims.list$beta.phi[,1]), nclass = 40, col = "gray", main = "Block 2", xlab = "Survival probability", xlim = c(0, 1))


# check convergence
hist(out$BUGSoutput$summary[,"Rhat"], nclass=8, main="Rhat", xlab="", las=1)
summary(out$BUGSoutput$summary[,"Rhat"])

# which parameters did not converged well?
out$BUGSoutput$summary[which(out$BUGSoutput$summary[,"Rhat"] > 1.1),]

# Evaluate fit
plot(out$BUGSoutput$sims.list$fit, out$BUGSoutput$sims.list$fit.new, main = "", xlab = "Discrepancy actual data", ylab = "Discrepancy replicate data", frame.plot = FALSE)
abline(0, 1, lwd = 2, col = "black")
mean(out$BUGSoutput$sims.list$fit.new > out$BUGSoutput$sims.list$fit)
mean(out$BUGSoutput$mean$fit) / mean(out$BUGSoutput$mean$fit.new)

# bayesian p-values is the proportion of points above the 1:1 line of equality
# bayesian p-values close to zero or one are suspicious, we've got 0.58 which is a good fit!
# ADD CODE
