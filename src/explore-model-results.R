# Explore results of multi-season occupancy model

# Codes based on book Bayesian Population Analysis using WinBUGS
# https://www.sciencedirect.com/book/9780123870209/bayesian-population-analysis-using-winbugs
# Codes adapted from site BPA with JAGS
# https://www.vogelwarte.ch/de/projekte/publikationen/bpa/code-for-running-bpa-using-jags

#----- 1 - Load libraries-----
library(here)
library(tidyverse)
library(ggplot2)

#----- 2 - Source files-----
source(here("bin", "figures.R"))

#----- 3 - Read and prepare data -----
out <- read_rds(here("results", "pobscura_model.rds"))


jags.data <- read_rds(here("data", "psophia_data_ijk.rds"))
jags.data$Ind <- as.numeric(ifelse(jags.data$block == 2, 1, 0)) # recode blocks, B1=0; B2=1)
attach(jags.data)
pobscura <- read_rds(here("data", "pobscura.rds"))

#----- 4 - Check fit and model results -----

# Summarize posteriors
#print(out, dig = 2)

# check convergence
hist(out$BUGSoutput$summary[,"Rhat"], nclass=8, main="Rhat", xlab="", las=1)
summary(out$BUGSoutput$summary[,"Rhat"])

# which parameters did not converged well?
out$BUGSoutput$summary[which(out$BUGSoutput$summary[,"Rhat"] > 1.1),]
dim(out$BUGSoutput$summary[which(out$BUGSoutput$summary[,"Rhat"] > 1.1),]) # how many failed to converge?


##----- Initial psi and p

# coefficients (predictor effects)
coef.function <- function(x) {
  coefs <- data.frame(x$BUGSoutput$summary[c("beta.psi[1]", "beta.psi[2]", "beta.psi[3]", "beta.psi[4]", "beta.psi[5]"),])
  coefs <- tibble(predictor=c("elevation", "distEdge","basalArea", "treeDensity", "recovery"),
                  coeff=row.names(coefs), mean=coefs$mean, sd=coefs$sd, lower=coefs$X2.5., upper=coefs$X97.5.,
                  Rhat=coefs$Rhat, n.eff=coefs$n.eff)
  .GlobalEnv$coefs <- coefs
  coefs
}
coef.function(out)


# plot significant effects: elevation on psi
predictor.effects.psi(out, pobscura$point.elevation, 1)
mtext("Elevation (masl)", side=1, line=3)

# plot significant effects: tree density on psi
predictor.effects.psi(out, pobscura$tree.density, 4)
mtext("Tree density (trees/ha)", side=1, line=3)

# plot significant effects: recovery on psi
predictor.effects.psi(out, pobscura$recovery, 5)
mtext("Recovery time (years)", side=1, line=3)


##----- temporal trends

# occupancy: psi
psiall.function <- function() {
  mcmc.sample <- out$BUGSoutput$n.sims
  array.psi <- array(NA, dim = c(nyear, mcmc.sample))
  for (i in 1:mcmc.sample){
    array.psi[,i] <- apply(out$BUGSoutput$sims.list$psi[i,,], 2, mean)
  }
  tibble(year=seq(2016,2020),
         mean=round(apply(array.psi, 1, mean), 2),
         sd=round(apply(array.psi, 1, sd), 2),
         lci=round(apply(array.psi, 1, quantile, prob=0.025), 2),
         uci=round(apply(array.psi, 1, quantile, prob=0.975), 2))
}
psiall.function()


# detection: p
pall.function <- function() {
  mcmc.sample <- out$BUGSoutput$n.sims
  array.p <- array(NA, dim = c(nyear, mcmc.sample))
  for (i in 1:mcmc.sample){
    array.p[,i] <- apply(out$BUGSoutput$sims.list$p[i,,,], 3, mean)
  }
  tibble(year=seq(2016,2020),
         mean=round(apply(array.p, 1, mean), 2),
         sd=round(apply(array.p, 1, sd), 2),
         lci=round(apply(array.p, 1, quantile, prob=0.025), 2),
         uci=round(apply(array.p, 1, quantile, prob=0.975), 2))
}
pall.function()



# yearly psi site means
psi.site <- tibble(array=(jags.data$Ind+1),
                   y.2016=apply(out$BUGSoutput$sims.list$psi[,,1], 2, mean),
                   y.2017=apply(out$BUGSoutput$sims.list$psi[,,2], 2, mean),
                   y.2018=apply(out$BUGSoutput$sims.list$psi[,,3], 2, mean),
                   y.2019=apply(out$BUGSoutput$sims.list$psi[,,4], 2, mean),
                   y.2020=apply(out$BUGSoutput$sims.list$psi[,,5], 2, mean))
psi.site
summary(psi.site)


# use previously created function to plot psi temporal trends
plot.psi.temporal.trends()

# growth rate (lambda):
growthr_table <- tibble(year=c(2016, 2017, 2018, 2019),
                        mean=apply(out$BUGSoutput$sims.list$growthr, 2, mean),
                        median=apply(out$BUGSoutput$sims.list$growthr, 2, median),
                        sd=apply(out$BUGSoutput$sims.list$growthr, 2, sd),
                        LCI=apply(out$BUGSoutput$sims.list$growthr, 2, quantile, prob=.025),
                        UCI=apply(out$BUGSoutput$sims.list$growthr, 2, quantile, prob=.975) )
growthr_table


# survival (phi)
phiall.function <- function() {
  mcmc.sample <- out$BUGSoutput$n.sims
  #array.phi <- array(NA, dim = c(nyear-1, mcmc.sample))
  #for (i in 1:mcmc.sample){
  #  array.phi[,i] <- apply(out$BUGSoutput$sims.list$phi, 2, mean)
  #}
  tibble(year=seq(2016,2019),
         mean=round(apply(out$BUGSoutput$sims.list$phi, 2, mean), 2),
         sd=round(apply(out$BUGSoutput$sims.list$phi, 2, sd), 2),
         lci=round(apply(out$BUGSoutput$sims.list$phi, 2, quantile, prob=0.025), 2),
         uci=round(apply(out$BUGSoutput$sims.list$phi, 2, quantile, prob=0.975), 2))
}
phiall.function()


# colonization (gamma)
gammaall.function <- function() {
  mcmc.sample <- out$BUGSoutput$n.sims
  #array.gamma <- array(NA, dim = c(nyear-1, mcmc.sample))
  #for (i in 1:mcmc.sample){
  #  array.gamma[,i] <- apply(out$BUGSoutput$sims.list$gamma, 2, mean)
  #}
  tibble(year=seq(2016,2019),
         mean=round(apply(out$BUGSoutput$sims.list$gamma, 2, mean), 2),
         sd=round(apply(out$BUGSoutput$sims.list$gamma, 2, sd), 2),
         lci=round(apply(out$BUGSoutput$sims.list$gamma, 2, quantile, prob=0.025), 2),
         uci=round(apply(out$BUGSoutput$sims.list$gamma, 2, quantile, prob=0.975), 2))
}
gammaall.function()



# plot survival trends with uncertainty
plot.phi.temporal.trends <- function() {
  mean_phi <- apply(out$BUGSoutput$sims.list$phi, 2, mean)
  mcmc.sample <- out$BUGSoutput$n.sims
  array.phi <- array(NA, dim = c(nyear-1, mcmc.sample))
  for (i in 1:mcmc.sample){
    array.phi[,i] <- out$BUGSoutput$sims.list$phi[i,]
  }
  # Plot for a subsample of MCMC draws
  sub.set <- sort(sample(1:mcmc.sample, size = 200))
  plot(2017:2020, mean_phi, main = "", ylab = expression(phi), xlab = "", 
       ylim=c(0, 1), type = "b", lwd = 2, las=1, xaxt="n")#, frame.plot = FALSE)
  for (i in sub.set){
    lines(2017:2020, array.phi[,i], type = "l", lwd = 0.1, col = "steelblue")
  }
  lines(2017:2020, mean_phi, type = "l", lwd = 1, col = "black")
  axis(1, at = c(2017, 2018, 2019, 2020), labels = seq(2017,2020))
}
plot.phi.temporal.trends()



# plot colonization trends with uncertainty
plot.gamma.temporal.trends <- function() {
  mean_gamma <- apply(out$BUGSoutput$sims.list$gamma, 2, mean)
  mcmc.sample <- out$BUGSoutput$n.sims
  array.gamma <- array(NA, dim = c(nyear-1, mcmc.sample))
  for (i in 1:mcmc.sample){
    array.gamma[,i] <- out$BUGSoutput$sims.list$gamma[i,]
  }
  # Plot for a subsample of MCMC draws
  sub.set <- sort(sample(1:mcmc.sample, size = 200))
  plot(2017:2020, mean_gamma, main = "", ylab = expression(gamma), xlab = "", 
       ylim=c(0, 1), type = "b", lwd = 2, las=1, xaxt="n")#, frame.plot = FALSE)
  for (i in sub.set){
    lines(2017:2020, array.gamma[,i], type = "l", lwd = 0.1, col = "steelblue")
  }
  lines(2017:2020, mean_gamma, type = "l", lwd = 1, col = "black")
  axis(1, at = c(2017, 2018, 2019, 2020), labels = seq(2017,2020))
}
plot.gamma.temporal.trends()

