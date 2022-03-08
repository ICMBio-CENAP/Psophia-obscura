
#-----
# Plot effects on psi with uncertainty
predictor.effects.psi <- function(x, original.predictor, coef) {
  
  #dev.off()
  
  ## Predict effect of logging on initial abundance with uncertainty
  mcmc.sample <- x$BUGSoutput$n.sims
  
  original.pred <- round(seq(min(original.predictor), max(original.predictor), length.out = 30),2)
  pred <- round((original.pred - mean(original.pred))/sd(original.pred), 2)
  #pred <- original.pred # it was already standardized
  psi.pred <- plogis( mean(x$BUGSoutput$sims.list$alpha.psi) +
                        mean(x$BUGSoutput$sims.list$beta.psi[, coef]) * pred )

  array.psi.pred <- array(NA, dim = c(length(pred), mcmc.sample))
  for (i in 1:mcmc.sample){
    array.psi.pred[,i] <- plogis( mean(x$BUGSoutput$sims.list$alpha.psi[i,]) +
                                    x$BUGSoutput$sims.list$beta.psi[i,coef] * pred )
  }
  
  # Plot for a subsample of MCMC draws
  sub.set <- sort(sample(1:mcmc.sample, size = 200))
  par(cex.axis = 1.5, cex.lab=2)
  #plot(original.pred, psi.pred, main = "", ylab = expression(Occupancy ~ (psi) ), xlab = "", yaxp= c(0, 1, 2),
  plot(original.pred, psi.pred, main = "", ylab = "", xlab = "", yaxp= c(0, 1, 2),
       ylim=c(0, 1), type = "l", lwd = 2, las=1, cex=1.5)#, frame.plot = FALSE)
  for (i in sub.set){
    lines(original.pred, array.psi.pred[,i], type = "l", lwd = 0.5, col = "grey")
  }
  lines(original.pred, psi.pred, type = "l", lwd = 0.8, col = "black")
  #mtext(expression(psi), side=2, line=3)
}
#predictor.effects(out, original.elevation, "a1")


#-----
# Plot effects on p with uncertainty
predictor.effects.p <- function(x, original.predictor, coef) {
  
  #dev.off()
  
  ## Predict effect of logging on initial abundance with uncertainty
  mcmc.sample <- x$BUGSoutput$n.sims
  
  original.pred <- round(seq(min(original.predictor), max(original.predictor), length.out = 30),2)
  pred <- round((original.pred - mean(original.pred))/sd(original.pred), 2)
  #pred <- original.pred # it was already standardized
  p.pred <- plogis( mean(x$BUGSoutput$sims.list$alpha.p) +
                        mean(x$BUGSoutput$sims.list$beta.p[, coef]) * pred )
  
  array.p.pred <- array(NA, dim = c(length(pred), mcmc.sample))
  for (i in 1:mcmc.sample){
    array.p.pred[,i] <- plogis( x$BUGSoutput$sims.list$alpha.p[i] +
                                    x$BUGSoutput$sims.list$beta.p[i,coef] * pred )
  }
  
  # Plot for a subsample of MCMC draws
  sub.set <- sort(sample(1:mcmc.sample, size = 200))
  
  plot(original.pred, p.pred, main = "", ylab = expression(p), xlab = "", 
       ylim=c(0, 1), type = "l", lwd = 2, las=1, frame.plot = FALSE)
  for (i in sub.set){
    lines(original.pred, array.p.pred[,i], type = "l", lwd = 3, col = "gray")
  }
  lines(original.pred, p.pred, type = "l", lwd = 2, col = "blue")
  #mtext(expression(p), side=2, line=3)
}
#predictor.effects(out, original.elevation, "a1")

#-----
# multipanel plot of effects
multipanel.4graphs <- function() {
  
  add_label_legend <- function(pos = "topleft", label, ...) {
    legend(pos, label, bty = "n", ...)
  }
  par(mfrow=c(2,2))
  par(mar=c(3,3,1,1))
  predictor.effects.psi(out, pobscura$point.elevation, 1)
  #mtext("Elevation (MASL)", side=1, line=3)
  add_label_legend("topleft", "A")
  predictor.effects.psi(out, pobscura$basal.area, 3)
  #mtext("Basal area (m²/ha)", side=1, line=3)
  add_label_legend("topleft", "B")
  predictor.effects.psi(out, pobscura$tree.density, 4)
  #mtext("Tree density (ind/ha)", side=1, line=3)
  add_label_legend("topleft", "C")
  predictor.effects.psi(out, pobscura$recovery, 5)
  #mtext("Recovery time (years)", side=1, line=3)
  add_label_legend("topleft", "D")

}


multipanel.3graphs <- function() {
  
  add_label_legend <- function(pos = "topleft", label, ...) {
    legend(pos, label, bty = "n", cex=1.5, ...)
  }
  par(mfrow=c(3,1))
  par(mar=c(5,5,1,2))
  predictor.effects.psi(out, pobscura$point.elevation, 1)
  mtext("Elevation (masl)", side=1, line=3)
  add_label_legend("topleft", "A")
  predictor.effects.psi(out, pobscura$tree.density, 4)
  mtext("Tree density (ind/ha)", side=1, line=3)
  add_label_legend("topleft", "B")
  mtext(expression(Occupancy ~ (psi) ), side=2, line=3)
  predictor.effects.psi(out, pobscura$recovery, 5)
  mtext("Recovery time (years)", side=1, line=3)
  add_label_legend("topleft", "C")

}
#multipanel.3graphs()

# alternative to multipanel: check posterior distribution of coefficients estimate
plot.coefs.posterior <- function() {
  par(mfrow=c(3,2))
  hist(out$BUGSoutput$sims.list$beta.psi[,1], xlab="Elevation (b1)", main="" )
  abline(v=0, col="red", lty=2)
  hist(out$BUGSoutput$sims.list$beta.psi[,2], xlab="Distance to edge (b2)", main="" )
  abline(v=0, col="red", lty=2)
  hist(out$BUGSoutput$sims.list$beta.psi[,3], xlab="Basal area (b3)", main="" )
  abline(v=0, col="red", lty=2)
  hist(out$BUGSoutput$sims.list$beta.psi[,4], xlab="Tree density (b4)", main="" )
  abline(v=0, col="red", lty=2)
  hist(out$BUGSoutput$sims.list$beta.psi[,5], xlab="Recovery time (b5)", main="" )
  abline(v=0, col="red", lty=2)
}




#-----
# Plot temporal trends with uncertainty
plot.psi.temporal.trends <- function() {
  mean_psi <- apply(out$BUGSoutput$sims.list$psi[,,], 3, mean)
  mcmc.sample <- out$BUGSoutput$n.sims
  array.psi <- array(NA, dim = c(nyear, mcmc.sample))
  for (i in 1:mcmc.sample){
    array.psi[,i] <- apply(out$BUGSoutput$sims.list$psi[i,,], 2, mean)
  }
  # Plot for a subsample of MCMC draws
  sub.set <- sort(sample(1:mcmc.sample, size = 200))
  par(cex.axis = 1, cex.lab=1.5)
  par(mar=c(5,5,1,1))
  plot(2016:2020, mean_psi, main = "", ylab = expression(Occupancy ~ (psi) ), xlab = "Year", 
       ylim=c(0, 1), type = "b", lwd = 2, las=1, xaxt="n")#, frame.plot = FALSE)
  for (i in sub.set){
    lines(2016:2020, array.psi[,i], type = "l", lwd = 0.5, col = "grey")
  }
  lines(2016:2020, mean_psi, type = "l", lwd = 0.8, col = "black")
  points(2016:2020, mean_psi, pch = 1, lwd = 1, col = "black")
  axis(1, at = c(2016, 2017, 2018, 2019, 2020), labels = seq(2016,2020))
}
#plot.psi.temporal.trends()
#---------------------------------------------------
# probably everything below is rubbish
#-----
Fig_effects.psi <- function(x) {
  par(mfrow=c(2,2))
  par(mar=c(3, 3, 2, 2))
  
  # effect of elevation on initial psi
  hist(x$BUGSoutput$sims.list$beta[,1], main="", xlab="", las=1)
  abline(v = 0, lty=2, lwd=3, col="red")
  #abline(v = quantile(out$BUGSoutput$sims.listbeta[,1], probs=c(0.025, 0.975)), lty=2)
  mtext("a", side = 3, line = -1.3, adj = 0.05, cex = 1.2, font = 2, col = "black")
  
  # effect of elevation on initial psi
  hist(x$BUGSoutput$sims.list$beta[,2], main="", xlab="", las=1)
  abline(v = 0, lty=2, lwd=3, col="red")
  #abline(v = quantile(out$BUGSoutput$sims.list$beta[,2], probs=c(0.025, 0.975)), lty=2)
  mtext("b", side = 3, line = -1.3, adj = 0.05, cex = 1.2, font = 2, col = "black")
  
  # effect of elevation on initial psi
  hist(x$BUGSoutput$sims.list$beta[,3], main="", xlab="", las=1)
  abline(v = 0, lty=2, lwd=3, col="red")
  #abline(v = quantile(out$BUGSoutput$sims.list$beta[,3], probs=c(0.025, 0.975)), lty=2)
  mtext("c", side = 3, line = -1.3, adj = 0.05, cex = 1.2, font = 2, col = "black")
  
  # effect of elevation on initial psi
  hist(x$BUGSoutput$sims.list$beta[,4], main="", xlab="", las=1)
  abline(v = 0, lty=2, lwd=3, col="red")
  #abline(v = quantile(out$BUGSoutput$sims.list$beta[,4], probs=c(0.025, 0.975)), lty=2)
  mtext("d", side = 3, line = -1.3, adj = 0.05, cex = 1.2, font = 2, col = "black")
  
}
#Fig_effects(out)

#-----
# Plot effects with uncertainty
predictor.effects.gam <- function(x, original.predictor, coef) {
  
  #dev.off()
  
  ## Predict effect of logging on initial abundance with uncertainty
  mcmc.sample <- x$BUGSoutput$n.sims
  
  original.pred <- round(seq(min(original.predictor), max(original.predictor), length.out = 30),2)
  pred <- round((original.pred - mean(original.pred))/sd(original.pred), 2)
  #pred <- original.pred # it was already standardized
  gam.pred <- plogis( mean(x$BUGSoutput$sims.list$alpha.gam) +
                        mean(x$BUGSoutput$sims.list$beta.gam[, coef]) * pred )
  
  array.gam.pred <- array(NA, dim = c(length(pred), mcmc.sample))
  for (i in 1:mcmc.sample){
    array.gam.pred[,i] <- plogis( x$BUGSoutput$sims.list$alpha.gam[i] +
                                    x$BUGSoutput$sims.list$beta.gam[i,coef] * pred )
  }
  
  # Plot for a subsample of MCMC draws
  sub.set <- sort(sample(1:mcmc.sample, size = 200))
  
  plot(original.pred, gam.pred, main = "", ylab = expression(gamma), xlab = "", 
       ylim=c(0, 1), type = "l", lwd = 2, las=1, frame.plot = FALSE)
  for (i in sub.set){
    lines(original.pred, array.gam.pred[,i], type = "l", lwd = 1, col = "gray")
  }
  lines(original.pred, gam.pred, type = "l", lwd = 2, col = "blue")
  #mtext(expression(gam), side=2, line=3)
}
#predictor.effects(out, original.elevation, "a1")

