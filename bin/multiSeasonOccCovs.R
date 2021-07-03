# https://sites.google.com/site/hierarchicalmodelingcourse/home/winbugs-models

model {

# Priors

psi0 ~ dnorm(0, 0.1)	# normal prior on psi intercept
psiDBH ~ dnorm(0, 0.1)  # slope of dbh for psi

g0 ~ dnorm(0, 0.1) 		# gamma intercept
gDBH ~ dnorm(0, 0.1)

e0 ~ dnorm(0, 0.1)		# epsilon intercept
eDBH ~ dnorm(0, 0.1)

p0 ~ dnorm(0, 0.1)		# p interecept
pWIND ~ dnorm(0, 0.1)



# Likelihood

for(i in 1:nSites) {
	logit(psi[i]) <- psi0 + psiDBH * dbh[i, 1]
	Z[i, 1] ~ dbin(psi[i], 1)
	for(t in 2:nYears) {
		logit(gamma[i, t-1]) <- g0 + gDBH * dbh[i, t]
		logit(epsilon[i, t-1]) <- e0 + eDBH * dbh[i, t]
		muZ[i, t] <- Z[i, t-1] * (1 - epsilon[i, t-1]) + 
			(1 - Z[i, t-1]) * gamma[i, t-1]
		Z[i, t] ~ dbin(muZ[i, t], 1)
		}
	for(t in 1:nYears) {	
		for(j in 1:nOccasions) {
			logit(p[i, t, j]) <- p0 + pWIND * wind[i, t, j]
			muY[i, t, j] <- Z[i, t] * p[i, t, j]
			y[i, t, j] ~ dbin(muY[i, t, j], 1)
			}
		}
	}

# Derived parameters (PAO = proportion of sites occupied each year)

for(t in 1:nYears) {
	PAO[t] <- sum(Z[,t]) / nSites
	}
	
	
}