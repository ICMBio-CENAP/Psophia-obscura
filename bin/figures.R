
Fig_effects <- function(x) {
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


# alternative plot with darkened histogram tails
#my_hist <- hist(out$BUGSoutput$sims.list$c1, nclass=15, main="", xlab="", las=1, ylab="")
#my_color= ifelse(my_hist$breaks < quantile(out$BUGSoutput$sims.list$c1, probs=c(0.025)), "blue" , ifelse (my_hist$breaks >= quantile(out$BUGSoutput$sims.list$c1, probs=c(0.975)), "blue", "grey" ))
#plot(my_hist, col=my_color, main="", las=1, xlab="")
#abline(v = quantile(out$BUGSoutput$sims.list$c1, probs=c(0.025, 0.975)), lty=2)


# Plot effects with uncertainty
predictor.effects <- function(x, original.predictor, coef) {
  
  #dev.off()

  ## Predict effect of logging on initial abundance with uncertainty
  mcmc.sample <- x$BUGSoutput$n.sims

  original.pred <- round(seq(min(original.predictor), max(original.predictor), length.out = 30),2)
  pred <- round((original.pred - mean(original.pred))/sd(original.pred), 2)
  #pred <- original.pred # it was already standardized
  psi.pred <- plogis( mean(x$BUGSoutput$sims.list$alpha.psi) +
                        mean(x$BUGSoutput$sims.list$beta[, coef]) * pred )
  
  array.psi.pred <- array(NA, dim = c(length(pred), mcmc.sample))
  for (i in 1:mcmc.sample){
    array.psi.pred[,i] <- plogis( x$BUGSoutput$sims.list$alpha.psi[i] +
                              x$BUGSoutput$sims.list$beta[i,coef] * pred )
  }
  
  # Plot for a subsample of MCMC draws
  sub.set <- sort(sample(1:mcmc.sample, size = 200))
  
  plot(original.pred, psi.pred, main = "", ylab = expression(psi), xlab = "", 
       ylim=c(0, 1), type = "l", lwd = 2, las=1, frame.plot = FALSE)
  for (i in sub.set){
    lines(original.pred, array.psi.pred[,i], type = "l", lwd = 1, col = "gray")
  }
  lines(original.pred, psi.pred, type = "l", lwd = 2, col = "blue")
  #mtext(expression(psi), side=2, line=3)
 }
#predictor.effects(out, original.elevation, "a1")



#-----
# Plot effects with uncertainty
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
    array.psi.pred[,i] <- plogis( x$BUGSoutput$sims.list$alpha.psi[i] +
                                    x$BUGSoutput$sims.list$beta.psi[i,coef] * pred )
  }
  
  # Plot for a subsample of MCMC draws
  sub.set <- sort(sample(1:mcmc.sample, size = 200))
  
  plot(original.pred, psi.pred, main = "", ylab = expression(psi), xlab = "", 
       ylim=c(0, 1), type = "l", lwd = 2, las=1, frame.plot = FALSE)
  for (i in sub.set){
    lines(original.pred, array.psi.pred[,i], type = "l", lwd = 1, col = "gray")
  }
  lines(original.pred, psi.pred, type = "l", lwd = 2, col = "blue")
  #mtext(expression(psi), side=2, line=3)
}
#predictor.effects(out, original.elevation, "a1")


#-----
# Plot effects with uncertainty
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
    lines(original.pred, array.p.pred[,i], type = "l", lwd = 1, col = "gray")
  }
  lines(original.pred, p.pred, type = "l", lwd = 2, col = "blue")
  #mtext(expression(p), side=2, line=3)
}
#predictor.effects(out, original.elevation, "a1")

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