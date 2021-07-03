
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






