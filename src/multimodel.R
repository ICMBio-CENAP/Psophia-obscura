library(here)
library(R2jags)
library(rjags)



# 13.5.2. Dynamic occupancy modeling in a real data set
# Read in the data and put it into 3D array
bdat <- read.table(file = here("data", "burnet.txt"), header = T)
str(bdat)

unique(bdat$site)

# site, day, count1, count2 etc

y <- array(NA, dim = c(95, 2, 7))	# 95 sites, 2 reps, 7 days

for (i in 1:7){
  sel.rows <- bdat$day == i
  y[,,i] <- as.matrix(bdat)[sel.rows, 3:4]
}
str(y)

# Convert counts to detection/nondetection data
y[y>0] <- 1

# Look at the number of sites with detections for each day
tmp <- apply(y, c(1,3), max, na.rm = TRUE)
tmp[tmp == "-Inf"] <- NA
apply(tmp, 2, sum, na.rm = TRUE)

# Bundle data
win.data <- list(y = y, nsite = dim(y)[1], nrep = dim(y)[2], nyear = dim(y)[3])

# Initial values
inits <- function(){ list(z = apply(y, c(1, 3), max))}

# Parameters monitored
params <- c("psi", "phi", "gamma", "p", "n.occ", "growthr", "turnover")  

# MCMC settings
ni <- 5000
nt <- 4
nb <- 1000
nc <- 3

# Call JAGS from R (BRT 1 min)
out1 <- jags(win.data, inits, params, "Dynocc.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(out1, dig = 3)

# Specify model in BUGS language
sink("Dynocc2.jags")
cat("
model {

# Specify priors
psi1 ~ dunif(0, 1)
for (k in 1:(nyear-1)){
   phi[k] ~ dunif(0, 1)
   gamma[k] ~ dunif(0, 1)
   }
p ~ dunif(0, 1)

# Both models at once
for (i in 1:nsite){
   z[i,1] ~ dbern(psi1)     # State model 1: Initial state
   muy[i,1] <- z[i,1]*p
   y[i,1] ~ dbin(muy[i,1], 2)
   for (k in 2:nyear){      # State model 2: State dynamics
      muZ[i,k] <- z[i,k-1]*phi[k-1] + (1-z[i,k-1])*gamma[k-1]
      z[i,k] ~ dbern(muZ[i,k])

      # Observation model
      muy[i,k] <- z[i,k]*p
      y[i,k] ~ dbin(muy[i,k], 2)
      } #k
   } #i

# Derived parameters: Sample and population occupancy, growth rate and turnover
psi[1] <- psi1
n.occ[1] <- sum(z[1:nsite,1])
for (k in 2:nyear){
   psi[k] <- psi[k-1]*phi[k-1] + (1-psi[k-1])*gamma[k-1]
   n.occ[k] <- sum(z[1:nsite,k])
   growthr[k-1] <- psi[k]/psi[k-1]
   turnover[k-1] <- (1 - psi[k-1]) * gamma[k-1]/psi[k]
   }
} # end of model
",fill = TRUE)
sink()

# Aggregate detections over reps within a day and bundle data
yy <- apply(y, c(1, 3), sum, na.rm = TRUE)
win.data <- list(y = yy, nsite = dim(yy)[1], nyear = dim(yy)[2])

# Initial values
inits <- function(){list(z = apply(y, c(1, 3), max))}

# Parameters monitored
params <- c("psi", "phi", "gamma", "p", "n.occ", "growthr", "turnover")  

# MCMC settings
ni <- 2500
nt <- 2
nb <- 500
nc <- 3

# Call JAGS from R (BRT 1 min)
out2 <- jags(win.data, inits, params, "Dynocc2.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posterior
print(out2, dig = 3)

DAY <- cbind(rep(1, out2$BUGSoutput$n.sims), rep(2, out2$BUGSoutput$n.sims), rep(3, out2$BUGSoutput$n.sims), rep(4, out2$BUGSoutput$n.sims), rep(5, out2$BUGSoutput$n.sims), rep(6, out2$BUGSoutput$n.sims), rep(7, out2$BUGSoutput$n.sims))
boxplot(out2$BUGSoutput$sims.list$psi ~ DAY, col = "gray", ylab = "Occupancy probability", xlab = "Day of survey", las = 1, frame.plot = FALSE)
apply(apply(y, c(1, 3), max), 2, function(x){sum(!is.na(x))})


