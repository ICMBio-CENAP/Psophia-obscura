
# EM CONSTRUCAO!!!!
# FALTA ADAPTAR O MODELO

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

SiteCovs <- pobscura[,42:48]
names(SiteCovs)
landCover <- SiteCovs[,1]
distWater <- SiteCovs[,2]
distEdge <- SiteCovs[,3]
elevation <- SiteCovs[,4]
basalArea <- SiteCovs[,5]
treeDensity <- SiteCovs[,6]
block <- SiteCovs[,7]

# check corr in SiteCovs
covars_correlations <- data.frame(cor(SiteCovs))
covars_correlations
write.csv(covars_correlations, here("results", "covars_correlations.csv"), row.names=FALSE)

# Standardize covariates
mean.landCover <- mean(landCover, na.rm = TRUE)
sd.landCover <- sd(landCover[!is.na(landCover)])
landCover <- (landCover-mean.landCover)/sd.landCover     # Standardise landCover
landCover[is.na(landCover)] <- 0               # Impute zeroes (means)
landCover <- round(landCover, 2)
landCover

mean.distWater <- mean(distWater, na.rm = TRUE)
sd.distWater <- sd(distWater[!is.na(distWater)])
distWater <- (distWater-mean.distWater)/sd.distWater     # Standardise distWater
distWater[is.na(distWater)] <- 0               # Impute zeroes (means)
distWater <- round(distWater, 2)
distWater

mean.distEdge <- mean(distEdge, na.rm = TRUE)
sd.distEdge <- sd(distEdge[!is.na(distEdge)])
distEdge <- (distEdge-mean.distEdge)/sd.distEdge     # Standardise distEdge
distEdge[is.na(distEdge)] <- 0               # Impute zeroes (means)
distEdge <- round(distEdge, 2)
distEdge

mean.elevation <- mean(elevation, na.rm = TRUE)
sd.elevation <- sd(elevation[!is.na(elevation)])
elevation <- (elevation-mean.elevation)/sd.elevation     # Standardise elevation
elevation[is.na(elevation)] <- 0               # Impute zeroes (means)
elevation <- round(elevation, 2)
elevation

mean.basalArea <- mean(basalArea, na.rm = TRUE)
sd.basalArea <- sd(basalArea[!is.na(basalArea)])
basalArea <- (basalArea-mean.basalArea)/sd.basalArea     # Standardise basalArea
basalArea[is.na(basalArea)] <- 0               # Impute zeroes (means)
basalArea <- round(basalArea, 2)
basalArea

mean.treeDensity <- mean(treeDensity, na.rm = TRUE)
sd.treeDensity <- sd(treeDensity[!is.na(treeDensity)])
treeDensity <- (treeDensity-mean.treeDensity)/sd.treeDensity     # Standardise treeDensity
treeDensity[is.na(treeDensity)] <- 0               # Impute zeroes (means)
treeDensity <- round(treeDensity, 2)
treeDensity


#----- 4 - Dynamic occupancy model without covariates -----

# Specify model in BUGS language
sink(here("bin", "Dynocc.jags"))
cat("
model {

# Specify priors
psi1 ~ dunif(0, 1)
for (k in 1:(nyear-1)){
   phi[k] ~ dunif(0, 1)
   gamma[k] ~ dunif(0, 1)
   p[k] ~ dunif(0, 1) 
   }
p[nyear] ~ dunif(0, 1)

# Ecological submodel: Define state conditional on parameters
for (i in 1:nsite){
   z[i,1] ~ dbern(psi1)
   for (k in 2:nyear){
      muZ[i,k]<- z[i,k-1]*phi[k-1] + (1-z[i,k-1])*gamma[k-1]
      z[i,k] ~ dbern(muZ[i,k])
      } #k
   } #i

# Observation model
for (i in 1:nsite){
   for (j in 1:nrep){
      for (k in 1:nyear){
         muy[i,j,k] <- z[i,k]*p[k]
         y[i,j,k] ~ dbern(muy[i,j,k])
         } #k
      } #j
   } #i

# Derived parameters: Sample and population occupancy, growth rate and turnover
psi[1] <- psi1
n.occ[1]<-sum(z[1:nsite,1])
for (k in 2:nyear){
   psi[k] <- psi[k-1]*phi[k-1] + (1-psi[k-1])*gamma[k-1]
   n.occ[k] <- sum(z[1:nsite,k])
   growthr[k-1] <- psi[k]/psi[k-1]                         # originally we had growthr[k]. JAGS seem to dislike vectoring going from 2..K.
   turnover[k-1] <- (1 - psi[k-1]) * gamma[k-1]/psi[k]
   }
}
",fill = TRUE)
sink()

# Bundle data
jags.data <- list(y = y, nsite = dim(y)[1], nrep = dim(y)[2], nyear = dim(y)[3])

# Initial values
zst <- apply(y, c(1, 3), max)	# Observed occurrence as inits for z
inits <- function(){ 
  zst[is.na(zst)] <- 1 # NAs will result in error (node inconsistent with parents)
  list(z = zst)
}

# Parameters monitored
params <- c("psi", "phi", "gamma", "p", "n.occ", "growthr", "turnover") 


# MCMC settings
ni <- 2500
nt <- 4
nb <- 500
nc <- 3

# Call JAGS from R (BRT 3 min)
out <- jags(jags.data, inits, params, here("bin", "Dynocc.jags"), n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)


# Summarize posteriors
print(out, dig = 2)
psiall <- paste("psi[", 1:K, "]", sep="")
print(out$BUGSoutput$summary[psiall, c(1, 2, 3, 7)], dig = 3)
phiall <- paste("phi[", 1:(K-1), "]", sep="")
print(out$BUGSoutput$summary[phiall, c(1, 2, 3, 7)], dig = 3)
gammaall <- paste("gamma[", 1:(K-1), "]", sep="")
print(out$BUGSoutput$summary[gammaall, c(1, 2, 3, 7)], dig = 3)
pall <- paste("p[", 1:K, "]", sep="")
print(out$BUGSoutput$summary[pall, c(1, 2, 3, 7)], dig = 3)

plot(1:K, out$BUGSoutput$mean$psi, type = "l", xlab = "Year", ylab = "Occupancy probability", col = "red", xlim = c(0,K+1), ylim = c(0,1), lwd = 2, lty = 1, frame.plot = FALSE, las = 1)
#lines(1:K, data$psi.app, type = "l", col = "black", lwd = 2)
#points(1:K, out$BUGSoutput$mean$psi, type = "l", col = "blue", lwd = 2)
segments(1:K, out$BUGSoutput$summary[psiall,3], 1:K, out$BUGSoutput$summary[psiall,7], col = "blue", lwd = 1)
