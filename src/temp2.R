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
  
  #for (k in 1:(nyear-1)){
  #  phi[k] ~ dunif(0, 1)
  #  gamma[k] ~ dunif(0, 1)
  #}

  # random block effects
  for (i in 1:nblock){
      mu.psi[i] ~ dnorm(0, tau.psi)
      #mu.p[i] ~ dnorm(0, tau.p)
      mu.phi[i] ~ dnorm(0, tau.phi)
      mu.gam[i] ~ dnorm(0, tau.gam)
    } #i
  # hyperpriors
  tau.psi ~ dgamma(0.001,0.001)
  #tau.p ~ dgamma(0.001,0.001)
  tau.phi ~ dgamma(0.001,0.001)
  tau.gam ~ dgamma(0.001,0.001)
  
  # random block and year effects for p
  for (i in 1:nblock){
    for (k in 1:nyear){
    mu.p[i,k] ~ dnorm(0, tau.p)
    } #k
  } #i
  tau.p ~ dgamma(0.001,0.001)
   

  # effects (betas)
  for (j in 1:3) {
      beta.psi[j] ~ dnorm(0, 0.01)
      #beta.p[j] ~ dnorm(0, 0.01)
      beta.phi[j] ~ dnorm(0, 0.01)
      beta.gam[j] ~ dnorm(0, 0.01)
  }

  ## Ecological model: state conditional on parameters
  for (i in 1:nsite){
    logit(psi[i,1]) <- mu.psi[block[i]] + beta.psi[1]*distEdge[i] + beta.psi[2]*basalArea[i] + beta.psi[3]*recovery[i]
    z[i,1] ~ dbern(psi[i,1])
    
    for (k in 2:nyear){
      logit(phi[i,k-1]) <- mu.phi[block[i]] + beta.phi[1]*distEdge[i] + beta.phi[2]*basalArea[i] + beta.phi[3]*recovery[i]
      logit(gamma[i,k-1]) <- mu.gam[block[i]] + beta.gam[1]*distEdge[i] + beta.gam[2]*basalArea[i] + beta.gam[3]*recovery[i]
      
      muZ[i,k] <- z[i,k-1]*phi[i,k-1] + (1-z[i,k-1])*gamma[i,k-1]
      z[i,k] ~ dbern(muZ[i,k])
      } #k
   } #i
   
   ## Observation model
   for (i in 1:nsite){
     for (j in 1:nrep){
       for (k in 1:nyear){
         logit(p[i,j,k]) <- mu.p[block[i],k]
         muy[i,j,k] <- z[i,k]*p[i,j,k]  # can only be detected if z=1
         y[i,j,k] ~ dbern(muy[i,j,k])
         } #k
     } #j
   } #i
   
   # Derived parameters: Sample and population occupancy, growth rate and turnover
   #n.occ[1] <- sum(z[1:nsite,1])
   #  psi[i,1] <- psi[i,1]
   for (i in 1:nsite){
     for (k in 2:nyear){
       psi[i,k] <- psi[i,k-1]*phi[i,k-1] + (1-psi[i,k-1])*gamma[i,k-1]
       #n.occ[k] <- sum(z[1:nsite,k])
       growthr[i,k-1] <- psi[i,k]/psi[i,k-1]  # originally we had growthr[k]. JAGS seem to dislike vectoring going from 2..K.
       turnover[i,k-1] <- (1 - psi[i,k-1]) * gamma[i,k-1]/psi[i,k]
     } # k
    } #i
}
",fill = TRUE)
sink()


#y <- rowSums(y, na.rm=TRUE)


# Bundle data
jags.data <- list(#y = y, block=block, nblock=length(unique(block)), nsite = length(y), nrep = 10, nyear = 4,
                  y = y, block=block, nblock=length(unique(block)), nsite = dim(y)[1], nrep = dim(y)[2], nyear = dim(y)[3],
                  elevation=point.elevation, distEdge=distEdge, distWater=distWater, basalArea=basalArea, recovery=recovery)

# Initial values
zst <- apply(y, c(1, 3), max)	# Observed occurrence as inits for z
inits <- function(){ 
  zst[is.na(zst)] <- 1 # NAs will result in error (node inconsistent with parents)
  list(z = zst)
}

# Parameters monitored
params <- c(#"psi", "phi", "gamma", "p", "n.occ", "growthr", "turnover",
  "psi", "phi", "gamma", "p", "growthr", "turnover", "n.occ",
  "mu.psi", "mu.p", "mu.phi", "mu.gam", "beta.psi", "beta.phi", "beta.gam") 


# MCMC settings
ni <- 100000
nt <- 100
nb <- 50000
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
coef.function <- function(x) {
  #coefs <- data.frame(x$BUGSoutput$summary[c("beta.psi[1]", "beta.psi[2]", "beta.psi[3]"),])
  #coefs <- data.frame(x$BUGSoutput$summary[c("beta.phi[1]", "beta.phi[2]", "beta.phi[3]"),]) 
  coefs <- data.frame(x$BUGSoutput$summary[c("beta.gam[1]", "beta.gam[2]", "beta.gam[3]"),]) 
  coefs <- tibble(predictor=c("edges", "basal.area", "recovery"),
                  coeff=row.names(coefs), mean=coefs$mean, lower=coefs$X2.5., upper=coefs$X97.5.,
                  Rhat=coefs$Rhat, n.eff=coefs$n.eff)
  #inds <- data.frame(x$BUGSoutput$summary[c("w[1]", "w[2]", "w[3]", "w[4]", "w[5]", "w[6]"),]) # w
  #inds <- tibble(w.mean=inds$mean)
  #coefs <- bind_cols(coefs, inds)
  .GlobalEnv$coefs <- coefs
  coefs
}
coef.function(out)


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

predictor.effects.psi(out, log(original.distEdge), 1)
mtext("Distance to edge (m)", side=1, line=3)

predictor.effects.psi(out, original.basalArea, 2)
mtext("Basal area (m²/ha)", side=1, line=3)

predictor.effects.psi(out, original.recovery, 3)
mtext("Recovery time (years)", side=1, line=3)

predictor.effects.gam(out, original.basalArea, 2)
mtext("Basal area (m²/ha)", side=1, line=3)


# save as jpeg:
jpeg(here("results", "basalArea_effect_psi.jpg"), width = 800, height = 400) # Open jpeg file
predictor.effects.psi(out, original.basalArea, 2)
mtext("Basal area (m²/ha)", side=1, line=3)
dev.off()

# save as jpeg:
jpeg(here("results", "recovery_effect_psi.jpg"), width = 800, height = 400) # Open jpeg file
predictor.effects.psi(out, original.recovery, 3)
mtext("Recovery time (years)", side=1, line=3)
dev.off()

# save as jpeg:
jpeg(here("results", "basalArea_effect_gamma.jpg"), width = 800, height = 400) # Open jpeg file
predictor.effects.gam(out, original.basalArea, 2)
mtext("Basal area (m²/ha)", side=1, line=3)
dev.off()

#####################################
#----- 4 - Dynamic occupancy model with site heterogeneity in all parameters -----

# Specify model in JAGS language
# using cerulean as a template
sink(here("bin", "cerulean.jags"))
cat("

model {
  
  psi ~ dunif(0,1)
  tau.gam ~ dgamma(.1,.1)
  tau.phi ~ dgamma(.1,.1)
  taup ~ dgamma(.1,.1)
  sigma.phi <- 1/sqrt(tau.phi)
  sigma.gam <- 1/sqrt(tau.gam)
  sigma.p <- sqrt(1/taup)
  
  mup.prob[1] ~ dunif(0,1)
  logit(mup[1]) <- mup.prob[1]
  for(i in 2:nyear){
    mup.prob[i] ~ dunif(0,1)
    logit(mup[i]) <- mup.prob[i]
    muphi.prob[i] ~ dunif(0,1)
    logit(muphi[i]) <- muphi.prob[i]
    mugam.prob[i] ~ dunif(0,1)
    logit(mugam[i]) <- mugam.prob[i]
  }
  
  for(i in 1:nsite){
    lp[i] ~ dnorm(0,taup)I(-12,12)
    for(t in 1:nyear){
      logit(p[i,t])<-mup[t]+lp[i]
    }
  }
  
  for(i in 1:nsite){
    z[i,1] ~ dbern(psi)
    lphi[i] ~ dnorm(0,tau.phi)I(-12,12)
    lgam[i] ~ dnorm(0,tau.gam)I(-12,12)
    for(t in 2:nyear){
      logit(gamma[i,t]) <- mugam[t] + lgam[i] 
      logit(phi[i,t]) <- muphi[t] + lphi[i] 
      muZ[i,t] <- z[i,t-1]*phi[i,t] + (1-z[i,t-1])*gamma[i,t]
      z[i,t] ~ dbern(muZ[i,t])
    }
  }
  
  for(i in 1:nsite){
    for (t in 1:nyear){
      Px[i,t] <- z[i,t]*p[i,t]
      x[i,t] ~ dbin(Px[i,t],50)
    }
  }
  
  for(i in 1:nyear){
    psi.year[i] <- sum(z[1:nsite,i])
  }
  for(i in 2:nyear){
    growthr[i-1] <- psi.year[i]/psi.year[i-1]
  }

}
",fill = TRUE)
sink()

# Bundle data for JAGS
jags.data <- list(y = y, nsite = dim(y)[1], nrep = dim(y)[2], nyear = dim(y)[3])


# Initial values
#inits <- function(){ list(z = apply(y, c(1, 3), max)) }
inits <- function(){ 
  z = apply(y, c(1, 3), max)
  z[is.na(z)] <- 1
  z <- list(z=z)
}

# Parameters monitored
params <- c("psi", "phi", "gamma", "p", "n.occ", "growthr", "turnover")  

# MCMC settings
ni <- 5000
nt <- 4
nb <- 1000
nc <- 3

out1 <- jags(jags.data, inits, params, here("bin", "cerulean.jags"), n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

# Summarize posteriors
print(out1, dig = 3)
# Estimated parameters:
# gamma = colonization probability
# growth = growth rate
# n.occ = number of occupied sites in a give primary period (year)
# p = detection probability
# phi = survival probability
# psi = occupancy probability
# turnover

# boxplot
YEAR <- cbind(rep(1, out1$BUGSoutput$n.sims), rep(2, out1$BUGSoutput$n.sims), rep(3, out1$BUGSoutput$n.sims), rep(4, out1$BUGSoutput$n.sims))
boxplot(out1$BUGSoutput$sims.list$psi ~ YEAR, col = "gray", ylab = "Occupancy probability", xlab = "Year", las = 1, frame.plot = FALSE)
apply(apply(y, c(1, 3), max), 2, function(x){sum(!is.na(x))})

# a better plot
occ.trends <- function() {
  # extract mean psi
  mean.psi <- apply(out1$BUGSoutput$sims.list$psi, 2, mean) # subset yearly means
  # lower CIs
  quants <- c(0.025)
  lower.psi <- apply( out1$BUGSoutput$sims.list$psi, 2 , quantile , probs = quants , na.rm = TRUE )
  # upper CIs
  quants2 <- c(0.95)
  upper.psi <- apply( out1$BUGSoutput$sims.list$psi, 2 , quantile , probs = quants2 , na.rm = TRUE )
  # join in a dataframe
  df1 <- data.frame(cbind(mean.psi, lower.psi, upper.psi)) # join them in a dataframe
  df1$year <- c(2016:2019) # as.numeric(rownames(df1))
  
  ggplot(data=df1, aes(x=year,y=mean.psi)) +
    geom_ribbon(data=df1, aes(x=year, ymin=lower.psi, ymax=upper.psi), fill="grey", alpha=0.8) +
    geom_line(size=0.5, alpha=1) +
    geom_point(size=1.5, alpha=0.5) +
    ylim(0,1) + 
    scale_x_continuous(breaks=c(2016,2017,2018,2019)) +
    theme(panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    theme(axis.text=element_text(size=12),axis.title=element_text(size=14)) + 
    #theme(axis.title.x = element_blank()) +
    #ggtitle() +
    theme(plot.title = element_text(hjust = 0.5)) +
    xlab("Years") + 
    ylab("Occupancy")
}

occ.trends()

# save as jpeg
jpeg(here("results", "occ_trends.jpg"), width = 600, height = 400) # Open jpeg file
occ.trends()
dev.off()

survival.trends <- function() {
  # extract mean psi
  mean.phi <- apply(out1$BUGSoutput$sims.list$phi, 2, mean) # subset yearly means
  # lower CIs
  quants <- c(0.025)
  lower.phi <- apply( out1$BUGSoutput$sims.list$phi, 2 , quantile , probs = quants , na.rm = TRUE )
  # upper CIs
  quants2 <- c(0.95)
  upper.phi <- apply( out1$BUGSoutput$sims.list$phi, 2 , quantile , probs = quants2 , na.rm = TRUE )
  # join in a dataframe
  df1 <- data.frame(cbind(mean.phi, lower.phi, upper.phi)) # join them in a dataframe
  df1$year <- as.numeric(rownames(df1))
  
  ggplot(data=df1, aes(x=year,y=mean.phi)) +
    geom_ribbon(data=df1, aes(x=year, ymin=lower.phi, ymax=upper.phi), fill="grey", alpha=0.8) +
    geom_line(size=0.5, alpha=1) +
    geom_point(size=1.5, alpha=0.5) +
    ylim(0,1) + 
    scale_x_continuous(breaks=c(2016,2017,2018,2019)) +
    theme(panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    theme(axis.text=element_text(size=12),axis.title=element_text(size=14)) + 
    #theme(axis.title.x = element_blank()) +
    #ggtitle() +
    theme(plot.title = element_text(hjust = 0.5)) +
    xlab("Years") + 
    ylab("Occupancy")
}

survival.trends()






