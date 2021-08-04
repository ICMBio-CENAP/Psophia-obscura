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
# if using 5-day occasion:
#y <- array(c(unlist(pobscura[,2:15]), unlist(pobscura[,16:29]), unlist(pobscura[,30:43]), unlist(pobscura[,44:57])), c(61, 14, 4))
# if using 10-day occasion:
y <- array(c(unlist(pobscura[,2:13]), unlist(pobscura[,14:25]), unlist(pobscura[,26:37]), unlist(pobscura[,38:49])), c(61, 12, 4))
str(y)

R <- dim(y)[1]
J <- dim(y)[2]
K <- dim(y)[3]

#j2016 <- as.numeric(rowSums(!is.na(y[,,1])))
#j2017 <- as.numeric(rowSums(!is.na(y[,,2])))
#j2018 <- as.numeric(rowSums(!is.na(y[,,3])))
#j2019 <- as.numeric(rowSums(!is.na(y[,,4])))
#nrep <- data.frame(cbind(j2016, j2017, j2018, j2019))
#head(nrep)

#y2016 <- rowSums(y[,,1], na.rm=TRUE)
#y2017 <- rowSums(y[,,2], na.rm=TRUE)
#y2018 <- rowSums(y[,,3], na.rm=TRUE)
#y2019 <- rowSums(y[,,4], na.rm=TRUE)
#y <- data.frame(cbind(y2016, y2017, y2018, y2019))
#head(y)

SiteCovs <- pobscura[,50:59]
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
original.bouts <- SiteCovs[,10]

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
hist(original.bouts, main="", xlab="Recovery time (years)")
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

bouts <- original.bouts
mean.bouts <- mean(bouts, na.rm = TRUE)
sd.bouts <- sd(bouts[!is.na(bouts)])
bouts <- (bouts-mean.bouts)/sd.bouts     # Standardise bouts
bouts[is.na(bouts)] <- 0               # Impute zeroes (means)
bouts <- round(bouts, 2)
bouts


psophia_data_ijk <- list(y=y, nsite=as.numeric(dim(y)[1]), nocc=as.numeric(dim(y)[2]), nyear=as.numeric(dim(y)[3]),
                     block=block, nblock=as.numeric(length(unique(block))),
                     elevation=point.elevation,
                     range.elevation=range.elevation,
                     distEdge=distEdge,
                     distWater=distWater,
                     basalArea=basalArea,
                     treeDensity=treeDensity,
                     recovery=recovery,
                     bouts=bouts)
str(psophia_data_ijk)

saveRDS(psophia_data_ijk, here("data", "psophia_data_ijk.rds"))


