# Based on book BPA with WinBUGS
# Codes adapted from site BPA with JAGS

#----- 1 - Load libraries-----
#library(dplyr)
#library(lubridate)
library(here)
library(tidyverse)
#library(R2jags)
#library(rjags)
library(ggplot2)

#----- 2 - Source files-----
source(here("bin", "figures.R"))

#----- 3 - Read and prepare data -----
#out <- read_rds(here("results", "pobscura_mod_3predictors_simplest.rds"))
out <- read_rds(here("results", "pobscura_mod_2021-08-02.rds"))


jags.data <- read_rds(here("data", "psophia_data.rds"))
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

# coefficients 
coef.function <- function(x) {
  coefs <- data.frame(x$BUGSoutput$summary[c("beta.psi[1]", "beta.psi[2]", "beta.psi[3]", "beta.psi[4]", "beta.psi[5]", "beta.psi[6]"),])
  coefs <- tibble(predictor=c("elevation", "distWater", "basalArea", "treeDensity", "distEdge", "recovery"),
                  coeff=row.names(coefs), mean=coefs$mean, lower=coefs$X2.5., upper=coefs$X97.5.,
                  Rhat=coefs$Rhat, n.eff=coefs$n.eff)
  .GlobalEnv$coefs <- coefs
  coefs
}
coef.function(out)

# plot significant effects: basal area on psi
#predictor.effects.psi(out, pobscura$basal.area, 3)
#mtext("Basal area (m²/ha)", side=1, line=3)
# save jpeg
#jpeg(here("results", "basalArea_effect_psi.jpg"), res=120, width = 800, height = 600)
#predictor.effects.psi(out, pobscura$basal.area, 3)
#mtext("Basal area (m²/ha)", side=1, line=3)
#dev.off()

# plot significant effects: distWater on psi
predictor.effects.psi(out, pobscura$dist.water, 2)
mtext("Distance to water (m)", side=1, line=3)
# save jpeg
jpeg(here("results", "distWater_effect_psi.jpg"), res=120, width = 800, height = 500)
predictor.effects.psi(out, pobscura$dist.water, 2)
mtext("Distance to water (m)", side=1, line=3)
dev.off()

# plot significant effects: basalArea on psi
predictor.effects.psi(out, pobscura$basal.area, 3)
mtext("Basal area (m²/ha)", side=1, line=3)
# save jpeg
jpeg(here("results", "basalArea_effect_psi.jpg"), res=120, width = 800, height = 500)
predictor.effects.psi(out, pobscura$basal.area, 3)
mtext("Basal area (m²/ha)", side=1, line=3)
dev.off()

# plot significant effects: recovery on psi
predictor.effects.psi(out, pobscura$recovery, 6)
mtext("Recovery time (years)", side=1, line=3)
# save jpeg
jpeg(here("results", "recovery_effect_psi.jpg"), res=120, width = 800, height = 500)
predictor.effects.psi(out, pobscura$recovery, 6)
mtext("Recovery time (years)", side=1, line=3)
dev.off()


# multipanel
add_label_legend <- function(pos = "topleft", label, ...) {
  legend(pos, label, bty = "n", ...)
}
par(mfrow=c(3,1))
par(mar=c(3,5,1,2))
predictor.effects.psi(out, pobscura$dist.water, 2)
mtext("Distance to water (m)", side=1, line=3)
add_label_legend("topleft", "A")
predictor.effects.psi(out, pobscura$basal.area, 3)
mtext("Basal area (m²/ha)", side=1, line=3)
add_label_legend("topleft", "B")
predictor.effects.psi(out, pobscura$recovery, 6)
mtext("Recovery time (years)", side=1, line=3)
add_label_legend("topleft", "C")
dev.off()

# save jpeg
jpeg(here("results", "Fig_2.jpg"), res=120, width = 600, height = 1200)
par(mfrow=c(3,1))
par(mar=c(3,5,1,2))
predictor.effects.psi(out, pobscura$dist.water, 2)
#mtext("Distance to water (m)", side=1, line=3)
add_label_legend("topleft", "A")
predictor.effects.psi(out, pobscura$basal.area, 3)
#mtext("Basal area (m²/ha)", side=1, line=3)
add_label_legend("topleft", "B")
predictor.effects.psi(out, pobscura$recovery, 6)
#mtext("Recovery time (years)", side=1, line=3)
add_label_legend("topleft", "C")
dev.off()


# alternative: check posterior distribution of coefficients estimate
jpeg(here("results", "all_coefficients.jpg"))
par(mfrow=c(2,2))
hist(out$BUGSoutput$sims.list$beta.psi[,1], xlab="Elevation effect on psi", main="" )
abline(v=0, col="red", lty=2)
hist(out$BUGSoutput$sims.list$beta.psi[,2], xlab="Basal area effect on psi", main="" )
abline(v=0, col="red", lty=2)
hist(out$BUGSoutput$sims.list$beta.psi[,3], xlab="Distance to edge effect on psi", main="" )
abline(v=0, col="red", lty=2)
hist(out$BUGSoutput$sims.list$beta.psi[,4], xlab="Recovery time effect on psi", main="" )
abline(v=0, col="red", lty=2)
dev.off()

##----- temporal trends

psiall <- round(apply(out$BUGSoutput$sims.list$psi, 3, mean), 2)
names(psiall) <- paste("psi[", 1:nyear, "]", sep="")
psiall

# p
pall <- round(apply(out$BUGSoutput$sims.list$p, 3, mean), 2)
names(pall) <- paste("p[", 1:nyear, "]", sep="")
pall

# yearly psi site means
psi.site <- tibble(array=(jags.data$Ind+1),
                   y.2016=apply(out$BUGSoutput$sims.list$psi[,,1], 2, mean),
                   y.2017=apply(out$BUGSoutput$sims.list$psi[,,2], 2, mean),
                   y.2018=apply(out$BUGSoutput$sims.list$psi[,,3], 2, mean),
                   y.2019=apply(out$BUGSoutput$sims.list$psi[,,4], 2, mean) )
psi.site
summary(psi.site)

means <- apply(psi.site[2:5], 2, mean)
lower <- apply(psi.site[2:5], 2, quantile, probs=c(0.025))
upper <- apply(psi.site[2:5], 2, quantile, probs=c(0.975))
plot(seq(2016,2019), means, type="b", ylim=c(0,1), las=1, xaxt = "n") 
axis(1, at = c(2016, 2017, 2018, 2019), labels = seq(2016,2019))
segments(c(2016, 2017, 2018, 2019), lower, c(2016, 2017, 2018, 2019), upper)


# a better version: plot occupancy trends with uncertainty
plot.psi.temporal.trends <- function() {
  mean_psi <- apply(out$BUGSoutput$sims.list$psi[,,], 3, mean)
  mcmc.sample <- out$BUGSoutput$n.sims
  array.psi <- array(NA, dim = c(nyear, mcmc.sample))
  for (i in 1:mcmc.sample){
    array.psi[,i] <- apply(out$BUGSoutput$sims.list$psi[i,,], 2, mean)
  }
  # Plot for a subsample of MCMC draws
  sub.set <- sort(sample(1:mcmc.sample, size = 200))
  plot(2016:2019, mean_psi, main = "", ylab = expression(psi), xlab = "", 
       ylim=c(0, 1), type = "l", lwd = 2, las=1, xaxt="n")#, frame.plot = FALSE)
  for (i in sub.set){
    lines(2016:2019, array.psi[,i], type = "l", lwd = 0.5, col = "gray")
  }
  lines(2016:2019, mean_psi, type = "l", lwd = 1, col = "black")
  axis(1, at = c(2016, 2017, 2018, 2019), labels = seq(2016,2019))
}
plot.psi.temporal.trends()


# save jpeg
jpeg(here("results", "psi_temporal_trends.jpg"), res=120, width = 800, height = 600)
plot.psi.temporal.trends()
dev.off()


# growth rate (lambda):
growthr_table <- tibble(year=c(2016, 2017, 2018),
                        mean=apply(out$BUGSoutput$sims.list$growthr, 2, mean),
                        median=apply(out$BUGSoutput$sims.list$growthr, 2, median),
                        LCI=apply(out$BUGSoutput$sims.list$growthr, 2, quantile, prob=.025),
                        UCI=apply(out$BUGSoutput$sims.list$growthr, 2, quantile, prob=.975) )
growthr_table

with(growthr_table, plot(year, median, type="b", ylim=c(0,4), las=1, xaxt = "n") )
axis(1, at = c(2016, 2017, 2018), labels = seq(2016,2018))
segments(c(2016, 2017, 2018), growthr_table$LCI, c(2016, 2017, 2018), growthr_table$UCI)
abline(h=1, lty=2)


# survival (phi)
phiall_mean <- round(apply(out$BUGSoutput$sims.list$phi, 2, mean), 2)
phiall_LCI <- round(apply(out$BUGSoutput$sims.list$phi, 2, quantile, prob=0.025), 2)
phiall_UCI <- round(apply(out$BUGSoutput$sims.list$phi, 2, quantile, prob=0.975), 2)
phiall <- data.frame(rbind(phiall_mean, phiall_LCI, phiall_UCI))
names(phiall) <- paste("phi[", 1:(nyear-1), "]", sep="")
phiall

# colonization (gamma)
gamma_mean <- round(apply(out$BUGSoutput$sims.list$gamma, 2, mean), 2)
gamma_LCI <- round(apply(out$BUGSoutput$sims.list$gamma, 2, quantile, prob=0.025), 2)
gamma_UCI <- round(apply(out$BUGSoutput$sims.list$gamma, 2, quantile, prob=0.975), 2)
gamma <- data.frame(rbind(gamma_mean, gamma_LCI, gamma_UCI))
names(gamma) <- paste("gamma[", 1:(nyear-1), "]", sep="")
gamma


#--------------------------


# table summarizing model parameters (overall means)
parameters.table <- function(x) {
  parameters_table <- tibble(parameter=c("psi", "p", "phi", "gamma", "turnover", "growthr"),
                             mean=c(round(mean(x$BUGSoutput$sims.list$psi), 2),
                                    round(mean(x$BUGSoutput$sims.list$p), 2),
                                    round(mean(x$BUGSoutput$sims.list$phi), 2),
                                    round(mean(x$BUGSoutput$sims.list$gamma), 2),
                                    round(mean(x$BUGSoutput$sims.list$turnover), 2),
                                    round(median(x$BUGSoutput$sims.list$growthr), 2)),
                             lower=c(round(quantile(x$BUGSoutput$sims.list$psi, probs=0.025), 3),
                                     round(quantile(x$BUGSoutput$sims.list$p, probs=0.025), 3),
                                     round(quantile(x$BUGSoutput$sims.list$phi, probs=0.025), 3),
                                     round(quantile(x$BUGSoutput$sims.list$gamma, probs=0.025), 3),
                                     round(quantile(x$BUGSoutput$sims.list$turnover, probs=0.025), 3),
                                     round(quantile(x$BUGSoutput$sims.list$growthr, probs=0.025), 3)),
                             upper= c(round(quantile(x$BUGSoutput$sims.list$psi, probs=0.975), 3),
                                      round(quantile(x$BUGSoutput$sims.list$p, probs=0.975), 3),
                                      round(quantile(x$BUGSoutput$sims.list$phi, probs=0.975), 3),
                                      round(quantile(x$BUGSoutput$sims.list$gamma, probs=0.975), 3),
                                      round(quantile(x$BUGSoutput$sims.list$turnover, probs=0.975), 3),
                                      round(quantile(x$BUGSoutput$sims.list$growthr, probs=0.975), 3)) )
  parameters_table <- data.frame(parameters_table)
  parameters_table
}
parameters.table(out)



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
  plot(2017:2019, mean_phi, main = "", ylab = expression(phi), xlab = "", 
       ylim=c(0, 1), type = "l", lwd = 2, las=1, xaxt="n", frame.plot = FALSE)
  for (i in sub.set){
    lines(2017:2019, array.phi[,i], type = "l", lwd = 1, col = "gray")
  }
  lines(2017:2019, mean_phi, type = "l", lwd = 2, col = "blue")
  axis(1, at = c(2017, 2018, 2019), labels = seq(2017,2019))
}
plot.phi.temporal.trends()

# save jpeg
jpeg(here("results", "phi_temporal_trends.jpg"), res=120, width = 800, height = 600)
plot.phi.temporal.trends()
dev.off()


# plot survival trends with uncertainty
plot.gamma.temporal.trends <- function() {
  mean_gamma <- apply(out$BUGSoutput$sims.list$gamma, 2, mean)
  mcmc.sample <- out$BUGSoutput$n.sims
  array.gamma <- array(NA, dim = c(nyear-1, mcmc.sample))
  for (i in 1:mcmc.sample){
    array.gamma[,i] <- out$BUGSoutput$sims.list$gamma[i,]
  }
  # Plot for a subsample of MCMC draws
  sub.set <- sort(sample(1:mcmc.sample, size = 200))
  plot(2017:2019, mean_gamma, main = "", ylab = expression(gamma), xlab = "", 
       ylim=c(0, 1), type = "l", lwd = 2, las=1, xaxt="n", frame.plot = FALSE)
  for (i in sub.set){
    lines(2017:2019, array.gamma[,i], type = "l", lwd = 1, col = "gray")
  }
  lines(2017:2019, mean_gamma, type = "l", lwd = 2, col = "blue")
  axis(1, at = c(2017, 2018, 2019), labels = seq(2017,2019))
}
plot.gamma.temporal.trends()

# save jpeg
jpeg(here("results", "gamma_temporal_trends.jpg"), res=120, width = 800, height = 600)
plot.gamma.temporal.trends()
dev.off()


# plot growth rate trends with uncertainty
plot.growthr.temporal.trends <- function() {
  mean_growthr <- apply(out$BUGSoutput$sims.list$growthr, 3, median)
  mcmc.sample <- out$BUGSoutput$n.sims
  array.growthr <- array(NA, dim = c(nyear-1, mcmc.sample))
  for (i in 1:mcmc.sample){
    array.growthr[,i] <- apply(out$BUGSoutput$sims.list$growthr[i,,], 2,median)
  }
  # Plot for a subsample of MCMC draws
  sub.set <- sort(sample(1:mcmc.sample, size = 200))
  plot(2017:2019, mean_growthr, main = "", ylab = expression(lambda), xlab = "",
       type = "l", lwd = 2, las=1, xaxt="n", frame.plot = FALSE)
       #ylim=c(0, 1), type = "l", lwd = 2, las=1, xaxt="n", frame.plot = FALSE)
  for (i in sub.set){
    lines(2017:2019, array.growthr[,i], type = "l", lwd = 1, col = "gray")
  }
  lines(2017:2019, mean_growthr, type = "l", lwd = 2, col = "blue")
  axis(1, at = c(2017, 2018, 2019), labels = seq(2017,2019))
}
plot.growthr.temporal.trends()

# save jpeg
#jpeg(here("results", "growthr_temporal_trends.jpg"), res=120, width = 800, height = 600)
#plot.growthr.temporal.trends()
#dev.off()


coef.function(out)

predictor.effects.psi(out, pobscura$basal.area, 2)
mtext("Basal area (m²/ha)", side=1, line=3)

# save jpeg
jpeg(here("results", "basalArea_effect_psi.jpg"), res=120, width = 800, height = 600)
predictor.effects.psi(out, pobscura$basal.area, 2)
mtext("Basal area (m²/ha)", side=1, line=3)
dev.off()



#predictor.effects.phi(out, pobscura$dist.water, 2)
#mtext("distance to water (km)", side=1, line=3)

#predictor.effects.gam(out, pobscura$basal.area, 2)
#mtext("Basal area (m²/ha)", side=1, line=3)

#predictor.effects.gam(out, pobscura$recovery, 4)
#mtext("Recovery time (years)", side=1, line=3)

#predictor.effects.gam(out, pobscura$recovery, 4)
#mtext("Recovery time (years)", side=1, line=3)

#predictor.effects.psi(out, pobscura$recovery, 4)
#mtext("Recovery time (years)", side=1, line=3)

#predictor.effects.p(out, pobscura$recovery, 4)
#mtext("Basal area (m²/ha)", side=1, line=3)

#par(mfrow = c(2, 1))
#hist(plogis(out$BUGSoutput$sims.list$alpha.phi), nclass = 40, col = "gray", main = "Block 1", xlab = "Survival probability", xlim = c(0, 1))
#hist(plogis(out$BUGSoutput$sims.list$alpha.phi + out$BUGSoutput$sims.list$beta.phi[,1]), nclass = 40, col = "gray", main = "Block 2", xlab = "Survival probability", xlim = c(0, 1))

