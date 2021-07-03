# Based on book BPA with WinBUGS
# Codes from site BPA with JAGS

# J. Andrew Royle and Marc Kéry. 2007. A Bayesian state-space formulation of dynamic occupancy models. Ecology 88:1813–1823.
# Supplement
# WinBUGS model specification for European Crossbill and Cerulean Warbler examples in the paper.
# Ecological Archives E088-108-S1.


# WinBUGS model specification for the Cerulean Warbler example in the paper
# In this example, heterogeneity in model parameters among sites was allowed, in addition to yearly variation
# The data structure is slightly different than in the previous example as they are the number of detections per sample unit per year (i.e., summed over all nrep secondary samples).

# x[i,t] are the observations, referenced by a two-dimensional matrix with indices i (site), and t (year).
# The model assumes that these are binomial counts based on a sample size of 50 (the number of replicate samples, or sub-samples along the BBS route)
# nyear : number of years of data (or primary periods)
# nsite : number of sites or sample locations

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

pobscura <- readRDS(here("data", "pobscura.rds"))
names(pobscura)
y <- array(c(unlist(pobscura[,2:11]), unlist(pobscura[,12:21]), unlist(pobscura[,22:31]), unlist(pobscura[,32:41])), c(61, 10, 4))
str(y)

R <- dim(y)[1]
J <- dim(y)[2]
K <- dim(y)[3]

SiteCovs <- pobscura[,42:50]
names(SiteCovs)
original.landCover <- SiteCovs[,1]
original.distWater <- SiteCovs[,2]#/1000 # convert from metres to km
original.distEdge <- SiteCovs[,3]#/1000 # convert from metres to km
original.point.elevation <- SiteCovs[,4]
original.range.elevation <- SiteCovs[,5]
original.basalArea <- SiteCovs[,6]
original.treeDensity <- SiteCovs[,7]
block <- SiteCovs[,8]
original.recovery <- SiteCovs[,9]

# check corr in SiteCovs
covars_correlations <- data.frame(cor(SiteCovs))
covars_correlations
write.csv(covars_correlations, here("results", "covars_correlations.csv"), row.names=FALSE)

# check distributions
hist(original.landCover, main="", xlab="Proportion of forest in 500 m buffer")
hist(log(original.landCover+0.1), main="", xlab="Log proportion of forest in 500 m buffer")
hist(original.distWater, main="", xlab="Distance to water (km)" )
hist(log(original.distWater + 0.1), main="", xlab="Log distance to water (km)" )
hist(original.distEdge, main="", xlab="Distance to pasture (km)" )
hist(log(original.distEdge + 0.1), main="", xlab="Log distance to pasture (km)" )
hist(original.point.elevation, main="", xlab="Elevation (m)")
hist(log(original.point.elevation), main="", xlab="Log elevation (m)")
hist(original.basalArea, main="", xlab="Basal area (m²/ha)")
hist(original.treeDensity, main="", xlab="Tree density (ind/ha)")
hist(log(original.treeDensity), main="", xlab="Log tree density (ind/ha)")
hist(original.recovery, main="", xlab="Recovery time (years)")
hist(log(original.recovery), main="", xlab="Log recovery time (years)")
# use log for distance to water and to edges

# Standardize covariates
landCover <- original.landCover
mean.landCover <- mean(landCover, na.rm = TRUE)
sd.landCover <- sd(landCover[!is.na(landCover)])
landCover <- (landCover-mean.landCover)/sd.landCover     # Standardise landCover
landCover[is.na(landCover)] <- 0               # Impute zeroes (means)
landCover <- round(landCover, 2)
landCover

distWater <- log(original.distWater)
mean.distWater <- mean(distWater, na.rm = TRUE)
sd.distWater <- sd(distWater[!is.na(distWater)])
distWater <- (distWater-mean.distWater)/sd.distWater     # Standardise distWater
distWater[is.na(distWater)] <- 0               # Impute zeroes (means)
distWater <- round(distWater, 2)
distWater

distEdge <- log(original.distEdge)
mean.distEdge <- mean(distEdge, na.rm = TRUE)
sd.distEdge <- sd(distEdge[!is.na(distEdge)])
distEdge <- (distEdge-mean.distEdge)/sd.distEdge     # Standardise distEdge
distEdge[is.na(distEdge)] <- 0               # Impute zeroes (means)
distEdge <- round(distEdge, 2)
distEdge

point.elevation <- original.point.elevation
mean.point.elevation <- mean(point.elevation, na.rm = TRUE)
sd.point.elevation <- sd(point.elevation[!is.na(point.elevation)])
point.elevation <- (point.elevation-mean.point.elevation)/sd.point.elevation     # Standardise point.elevation
point.elevation[is.na(point.elevation)] <- 0               # Impute zeroes (means)
point.elevation <- round(point.elevation, 2)
point.elevation

range.elevation <- original.range.elevation
mean.range.elevation <- mean(range.elevation, na.rm = TRUE)
sd.range.elevation <- sd(range.elevation[!is.na(range.elevation)])
range.elevation <- (range.elevation-mean.range.elevation)/sd.range.elevation     # Standardise range.elevation
range.elevation[is.na(range.elevation)] <- 0               # Impute zeroes (means)
range.elevation <- round(range.elevation, 2)
range.elevation

basalArea <- original.basalArea
mean.basalArea <- mean(basalArea, na.rm = TRUE)
sd.basalArea <- sd(basalArea[!is.na(basalArea)])
basalArea <- (basalArea-mean.basalArea)/sd.basalArea     # Standardise basalArea
basalArea[is.na(basalArea)] <- 0               # Impute zeroes (means)
basalArea <- round(basalArea, 2)
basalArea

treeDensity <- original.treeDensity
mean.treeDensity <- mean(treeDensity, na.rm = TRUE)
sd.treeDensity <- sd(treeDensity[!is.na(treeDensity)])
treeDensity <- (treeDensity-mean.treeDensity)/sd.treeDensity     # Standardise treeDensity
treeDensity[is.na(treeDensity)] <- 0               # Impute zeroes (means)
treeDensity <- round(treeDensity, 2)
treeDensity

recovery <- original.recovery
mean.recovery <- mean(recovery, na.rm = TRUE)
sd.recovery <- sd(recovery[!is.na(recovery)])
recovery <- (recovery-mean.recovery)/sd.recovery     # Standardise recovery
recovery[is.na(recovery)] <- 0               # Impute zeroes (means)
recovery <- round(recovery, 2)
recovery


#----- 4 - Dynamic occupancy model with covariates -----

# Specify model in BUGS language
sink(here("bin", "Dynocc_covariates.jags"))
cat("
model {

  ## Priors
  
  #psi ~ dunif(0,1)
  tau.gam ~ dgamma(.1,.1)
  tau.phi ~ dgamma(.1,.1)
  #taup ~ dgamma(.1,.1)
  #sigma.phi <- 1/sqrt(tau.phi)
  #sigma.gam <- 1/sqrt(tau.gam)
  #sigma.p <- sqrt(1/taup)
  
  #for(i in 2:nyear){
  #  mup.prob[i] ~ dunif(0,1)
  #  logit(mup[i]) <- mup.prob[i]
  #  muphi.prob[i] ~ dunif(0,1)
  #  logit(muphi[i])<-muphi.prob[i]
  #  mugam.prob[i]~dunif(0,1)
  #  logit(mugam[i])<-mugam.prob[i]
  #}
  
  # random site effects (intercept) for psi
  for (i in 1:nsite){
    alpha.psi[i] ~ dnorm(mu.site, tau.site)
  } #i
  mu.site ~ dnorm(0, 0.5)
  tau.site <- 1/(sd.site*sd.site)
  sd.site ~ dunif(0,5)
  
  # random site effects for p
  for (i in 1:nsite){
    alpha.p[i] ~ dnorm(mu.p.site, tau.p.site)
    } #i
  mu.p.site ~ dnorm(0, 0.5)
  tau.p.site <- 1/(sd.p.site*sd.p.site)
  sd.p.site ~ dunif(0,5)
  
  # predictor coefficients (psi)
  for (j in 1:3) {
      beta.psi[j] ~ dnorm(0, 0.01)
  }
  
  # predictor coefficients (gamma)
  for (j in 1:3) {
      beta.gam[j] ~ dnorm(0, 0.01)
  }
  
  # predictor coefficients (gamma)
  for (j in 1:3) {
      beta.phi[j] ~ dnorm(0, 0.01)
  }


  # Ecological model: state conditional on parameters
  for(i in 1:nsite){
    logit(psi[i,1]) <- alpha.psi[i] + beta.psi[1]*distEdge[i] + beta.psi[2]*basalArea[i] + beta.psi[3]*recovery[i]
    z[i,1] ~ dbern(psi[i,1])
    lphi[i] ~ dnorm(0,tau.phi)#I(-12,12)
    lgam[i] ~ dnorm(0,tau.gam)#I(-12,12)
    
    for(t in 2:nyear){
      #logit(gamma[i,t]) <- mugam[t] + lgam[i] + alpha.gam[i] + beta.gam[1]*distEdge[i] + beta.gam[2]*basalArea[i] + beta.gam[3]*recovery[i]
      #logit(phi[i,t]) <-   muphi[t] + lphi[i] + alpha.phi[i] + beta.phi[1]*distEdge[i] + beta.phi[2]*basalArea[i] + beta.phi[3]*recovery[i]
      logit(gamma[i]) <- lgam[i] + beta.gam[1]*distEdge[i] + beta.gam[2]*basalArea[i] + beta.gam[3]*recovery[i]
      logit(phi[i]) <-   lphi[i] + beta.phi[1]*distEdge[i] + beta.phi[2]*basalArea[i] + beta.phi[3]*recovery[i]
      
      muZ[i,t] <- z[i,t-1]*phi[i] + (1-z[i,t-1])*gamma[i]
      z[i,t] ~ dbern(muZ[i,t])
    }
  }
  
  # Observation model
  for (i in 1:nsite){
    for (j in 1:nrep){
      for (t in 1:nyear){
        logit(p[i,j,t]) <- alpha.p[i]
        muy[i,j,t] <- z[i,t]*p[i,j,t]  # can only be detected if z=1
        y[i,j,t] ~ dbern(muy[i,j,t])
        } #k
      } #j
    } #i

  # Derived parameters
  for(i in 1:nyear){
    psi.year[i] <- sum(z[1:nsite,i])
  }
  #for(i in 2:nyear){
  #  growthr[i-1]<-psi.year[i]/psi.year[i-1]
  #}
  for(i in 1:nsite){
     for(t in 2:nyear){
       psi[i,t] <- psi[i,t-1]*phi[i] + (1-psi[i,t-1])*gamma[i]
       growthr[i,t-1] <- psi[i,t]/psi[i,t-1]  # originally we had growthr[k]. JAGS seem to dislike vectoring going from 2..K.
       turnover[i,t-1] <- (1 - psi[i,t-1]) * gamma[i]/psi[i,t]
     } # k
    } #i

}
",fill = TRUE)
sink()

# Bundle data
jags.data <- list(y = y, nsite = dim(y)[1], nrep = dim(y)[2], nyear = dim(y)[3],
                  elevation=point.elevation, distEdge=distEdge, distWater=distWater, basalArea=basalArea, recovery=recovery)

# Initial values
zst <- apply(y, c(1, 3), max)	# Observed occurrence as inits for z
inits <- function(){ 
  zst[is.na(zst)] <- 1 # NAs will result in error (node inconsistent with parents)
  list(z = zst)
}

# Parameters monitored
params <- c(#"psi", "phi", "gamma", "p", "n.occ", "growthr", "turnover",
  "psi", "phi", "gamma", "p", "growthr", "turnover",
  "alpha.psi", "alpha.p", "beta.psi", "beta.p", "beta.gam", "beta.phi", "a", "beta") 


# MCMC settings
ni <- 25000
nt <- 100
nb <- 250
nc <- 3

# Call JAGS from R (BRT 3 min)
out <- jags(jags.data, inits, params, here("bin", "Dynocc_covariates.jags"), n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)


saveRDS(out, here("results", "pobscura_model.rds"))

# Summarize posteriors
print(out, dig = 2)
# psi now has one more dimension so original code doesnt work
# however the code still works for phi, gamma and p:
#psiall <- paste("psi[", 1:K, "]", sep="")
#print(out$BUGSoutput$summary[psiall, c(1, 2, 3, 7)], dig = 3)
phiall <- paste("phi[", 1:(K-1), "]", sep="")
print(out$BUGSoutput$summary[phiall, c(1, 2, 3, 7)], dig = 3)
gammaall <- paste("gamma[", 1:(K-1), "]", sep="")
print(out$BUGSoutput$summary[gammaall, c(1, 2, 3, 7)], dig = 3)
pall <- paste("p[", 1:K, "]", sep="")
print(out$BUGSoutput$summary[pall, c(1, 2, 3, 7)], dig = 3)

#plot(1:K, out$BUGSoutput$mean$psi, type = "l", xlab = "Year", ylab = "Occupancy probability", col = "red", xlim = c(0,K+1), ylim = c(0,1), lwd = 2, lty = 1, frame.plot = FALSE, las = 1)
#segments(1:K, out$BUGSoutput$summary[psiall,3], 1:K, out$BUGSoutput$summary[psiall,7], col = "blue", lwd = 1)


# check convergence
hist(out$BUGSoutput$summary[,"Rhat"], nclass=8, main="Rhat", xlab="", las=1)
summary(out$BUGSoutput$summary[,"Rhat"])

# which parameters did not converged well?
out$BUGSoutput$summary[which(out$BUGSoutput$summary[,"Rhat"] > 1.1),]

# Evaluate fit
# bayesian p-values is the proportion of points above the 1:1 line of equality
# bayesian p-values close to zero or one are suspicious, we've got 0.58 which is a good fit!
# ADD CODE

# coefficients 
coef.psi <- function(x) {
  coefs <- data.frame(x$BUGSoutput$summary[c("beta.psi[1]", "beta.psi[2]", "beta.psi[3]"),]) 
  coefs <- tibble(predictor=c("edges", "basal.area", "recovery"),
                  coeff=row.names(coefs), mean=coefs$mean, lower=coefs$X2.5., upper=coefs$X97.5.,
                  Rhat=coefs$Rhat, n.eff=coefs$n.eff)
  #inds <- data.frame(x$BUGSoutput$summary[c("w[1]", "w[2]", "w[3]", "w[4]", "w[5]", "w[6]"),]) # w
  #inds <- tibble(w.mean=inds$mean)
  #coefs <- bind_cols(coefs, inds)
  .GlobalEnv$coefs <- coefs
  coefs
}
coef.psi(out)

# coefficients 
coef.p <- function(x) {
  coefs <- data.frame(x$BUGSoutput$summary[c("beta.p[1]", "beta.p[2]", "beta.p[3]"),]) 
  coefs <- tibble(predictor=c("edges", "basal.area", "recovery"),
                  coeff=row.names(coefs), mean=coefs$mean, lower=coefs$X2.5., upper=coefs$X97.5.,
                  Rhat=coefs$Rhat, n.eff=coefs$n.eff)
  #inds <- data.frame(x$BUGSoutput$summary[c("w[1]", "w[2]", "w[3]", "w[4]", "w[5]", "w[6]"),]) # w
  #inds <- tibble(w.mean=inds$mean)
  #coefs <- bind_cols(coefs, inds)
  .GlobalEnv$coefs <- coefs
  coefs
}
coef.p(out)



# table summarizing model parameters
parameters.table <- function(x) {
  parameters_table <- tibble(parameter=c("psi", "p", "phi", "gamma"),
                             mean=c(round(mean(x$BUGSoutput$sims.list$psi), 2),
                                    round(mean(x$BUGSoutput$sims.list$p), 2),
                                    round(mean(x$BUGSoutput$sims.list$phi), 2),
                                    round(mean(x$BUGSoutput$sims.list$gamma), 2) ),
                             lower=c(round(quantile(x$BUGSoutput$sims.list$psi, probs=0.025), 3),
                                     round(quantile(x$BUGSoutput$sims.list$p, probs=0.025), 3),
                                     round(quantile(x$BUGSoutput$sims.list$phi, probs=0.025), 3),
                                     round(quantile(x$BUGSoutput$sims.list$gamma, probs=0.025), 3) ),
                             upper= c(round(quantile(x$BUGSoutput$sims.list$psi, probs=0.975), 3),
                                      round(quantile(x$BUGSoutput$sims.list$p, probs=0.975), 3),
                                      round(quantile(x$BUGSoutput$sims.list$phi, probs=0.975), 3),
                                      round(quantile(x$BUGSoutput$sims.list$gamma, probs=0.975), 3) ) )
  parameters_table <- data.frame(parameters_table)
  parameters_table
}

parameters.table(out)

# Figures

# template save as jpeg:
#jpeg(here("results", "basalArea_effect.jpg"), width = 800, height = 400) # Open jpeg file
#plot.basalArea()
#dev.off()


Fig_effects(out)

dev.off()
coef.function(out)

#predictor.effects(out, original.point.elevation, 1)
#mtext("Elevation (m)", side=1, line=3)

predictor.effects(out, log(original.distEdge), 1)
mtext("Distance to edge (m)", side=1, line=3)

#predictor.effects(out, log(original.distWater), 3)
#mtext("Distance to water (m)", side=1, line=3)

predictor.effects(out, original.basalArea, 2)
mtext("Basal area (m²/ha)", side=1, line=3)

predictor.effects.psi(out, original.recovery, 3)
mtext("Recovery time (years)", side=1, line=3)



# save as jpeg:
jpeg(here("results", "basalArea_effect.jpg"), width = 800, height = 400) # Open jpeg file
predictor.effects(out, original.basalArea, 4)
mtext("Basal area (m²/ha)", side=1, line=3)
dev.off()

