# Part2-Multi-season Bayesian occupancy model
# Elildo Carvalho Jr @ ICMBio/CENAP 2020-08-17
# code based on templates from the book "Bayesian population analysis using WinBUGS - a hierarchical perspective" by Marc Kéry & Michael Schaub (2012, Academic Press)
# and on templates from the JAGS translation available at:
# https://www.vogelwarte.ch/de/projekte/publikationen/bpa/code-for-running-bpa-using-jags

#----- 1 - Load libraries-----
#library(dplyr)
#library(lubridate)
library(here)
library(R2jags)
library(rjags)
library(ggplot2)


#----- 2 - Source files-----
#source(here("bin", "ahumada_codes.R"))
#source(here("bin", "f-matrix-creator-experimental-probably-ok-but-need-check.R"))
#source(here("src", "fix_species_names.R")) # fix some names and remove "false" species


#----- 3 - Read and prepare data -----
Pobscura <- readRDS(here("data", "Pobscura.rds"))
y <- array(c(unlist(Pobscura[,2:12]), unlist(Pobscura[,13:23]), unlist(Pobscura[,24:34]), unlist(Pobscura[,35:45])), c(61, 11, 4))
str(y)
SiteCovs <- Pobscura[,46:51]

# Look at the number of sites with detection for each day
tmp <- apply(y, c(1,3), max, na.rm = TRUE)
tmp[tmp == "-Inf"] <- NA
apply(tmp, 2, sum, na.rm = TRUE)


#----- 4 - mod0, dynamic occupancy with no covariates -----

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


# Specify model in JAGS language
sink(here("bin", "Dynocc_mod0.jags"))
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



# Call JAGS from R (BRT 1 min)
out1 <- jags(jags.data, inits, params, here("bin", "Dynocc.jags_mod0"), n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory =  getwd())

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



#----- 5 - mod1, dynamic occupancy with no covariates -----
# similar to mod0 but with GLM link functions for covariates
# using only siteCovs for the first year 

# Bundle data for JAGS
jags.data.mod1 <- list(y = y, nsite = dim(y)[1], nrep = dim(y)[2], nyear = dim(y)[3], elevation = SiteCovs$elevation, slope=SiteCovs$slope, dist.water=SiteCovs$dist.water)

# Initial values
#inits <- function(){ list(z = apply(y, c(1, 3), max)) }
inits.mod1 <- function(){ 
   #z = matrix(0, 61, 4)
   #z = apply(z, c(1, 2), function(x) sample(c(0, 1), 1))
   z = apply(y, c(1, 3), max)
   z[is.na(z)] <- 1
   list(z=z, alpha.occ=runif(1,-3,3), beta1.occ=runif(1,-3,3), beta2.occ=runif(1,-3,3), beta3.occ=runif(1,-3,3))
}

# Parameters monitored
params.mod1 <- c("psi", "phi", "gamma", "p", "n.occ", "growthr", "turnover", "alpha.occ", "beta1.occ", "beta2.occ", "beta3.occ")  

# MCMC settings
ni <- 5000
nt <- 4
nb <- 1000
nc <- 3

# Specify model in JAGS language
sink(here("bin", "Dynocc_mod1.jags"))
cat("
model {

# Specify priors
#psi1 ~ dunif(0, 1)

for (k in 1:(nyear-1)){
   phi[k] ~ dunif(0, 1)
   gamma[k] ~ dunif(0, 1)
   p[k] ~ dunif(0, 1) 
   }
p[nyear] ~ dunif(0, 1)

# priors for GLM occupancy 1st year
alpha.occ ~ dunif(-10, 10)
beta1.occ ~ dunif(-10, 10) # elevation
beta2.occ ~ dunif(-10, 10) # slope
beta3.occ ~ dunif(-10, 10) # dist.water

# Ecological submodel: Define state conditional on parameters
for (i in 1:nsite){
   z[i,1] ~ dbern(psi1[i])
   logit(psi1[i]) <- alpha.occ + beta1.occ*elevation[i] + beta2.occ*slope[i] + beta3.occ*dist.water[i]
   
   
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

for (i in 1:nsite){
   psi[i,1] <- psi1[i]
   }

n.occ[1] <- sum(z[1:nsite,1])

for (k in 2:nyear){
   n.occ[k] <- sum(z[1:nsite,k])
for (i in 1:nsite){
   psi[i,k] <- psi[i,k-1]*phi[k-1] + (1-psi[i,k-1])*gamma[k-1]
#   growthr[k-1] <- psi[i,k]/psi[i,k-1]           # originally we had growthr[k]. JAGS seem to dislike vectoring going from 2..K.
#   turnover[k-1] <- (1 - psi[i,k-1]) * gamma[k-1]/psi[i,k]
       }# year
   }# site
}
",fill = TRUE)
sink()


# Call JAGS from R (BRT 1 min)
#out.mod1 <- jags(jags.data.mod1, inits.mod1, params.mod1, here("bin", "Dynocc_mod1.jags"), n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
out.mod1 <- jags(jags.data.mod1, ,params.mod1, here("bin", "Dynocc_mod1.jags"), n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(out.mod1, dig = 3)
