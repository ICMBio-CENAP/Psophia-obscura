#-------Vegetation structure at TEAM sites, Gurupi Biological Reserve-------------

## ----- Load libraries -----

library(tidyverse)
library(here)


## ----- Read and fix data -----

veget <- read.csv(here("data", "rbg_veget_08mar2018.csv"), header=TRUE, sep=",")
colnames(veget)[1] <- "Camera.Trap.Name"
head(veget)

veget$dist.m <- veget$dist/100
veget$dap <- veget$cap/pi # criando nova coluna com DAP calculado a partir de CAP
veget$dap.m <- veget$dap*0.01
veget$ba <- pi*((veget$dap.m/2)^2) # cria nova coluna com area basal de cada arvore

veget$rep[veget$rep==""]  <- NA 
veget$id.ponto <- paste(veget$trilha, veget$quadrante, veget$ponto, sep=".") # criando id para pontos

dim(veget)
veget <- subset(veget, is.na((veget$rep)))
dim(veget)

veget <- veget[,c("Camera.Trap.Name", "id.ponto", "dist.m", "ba")]
head(veget)

# input mean values to NAs
veget$dist.m[which(is.na(veget$dist.m))] <- mean(veget$dist.m, na.rm=TRUE)
veget$ba[which(is.na(veget$ba))] <- mean(veget$ba, na.rm=TRUE)

# find large outlier trees
#tail(sort(veget$ba), n=20)

# a few exceptionally large trees may bias estimates
# this is not be a problem for the method
# but since our sample size was not so large it might be advisable to reduce their effect
# furthermore some large values may be due to typing error
# find largest trees and impute new values
hist(veget$ba)
hist(veget$ba, nclass=50, xlim=c(0,2))
hist(veget$ba, nclass=50, xlim=c(0,1))
veget[which(veget$ba > 1),] # only 13 individual trees has basal area > 1 m2
range(veget$ba)
mean(veget$ba)
median(veget$ba)
quantile(veget$ba, probs=c(0.975))
# replace the value of trees above the 95% quantile by the quantile itself 
nrow(veget)
nrow(veget[which(veget$ba > quantile(veget$ba, probs=c(0.975))),]) # 104 largest trees
veget[which(veget$ba > quantile(veget$ba, probs=c(0.975))),"ba"] <- quantile(veget$ba, probs=c(0.975))
hist(veget$ba)


# check if there are exceptionally large distances
# this may be due to typing errors
hist(veget$dist.m) # seems ok
range(veget$dist.m)
mean(veget$dist.m)
median(veget$dist.m)
quantile(veget$dist.m, probs=c(0.025))
quantile(veget$dist.m, probs=c(0.975))
sort(veget$dist.m)
# replace distances = 0 by 0.1
nrow(veget[which(veget$dist.m < 0.1),]) # 199 trees have zero distance
veget[which(veget$dist.m < 0.1),"dist.m"] <- 0.1
# replace distances > 95% quantile by the quantile value
nrow(veget[which(veget$dist.m > quantile(veget$dist.m, probs=c(0.975))),]) # 94 trees above distance quantile
veget[which(veget$dist.m > quantile(veget$dist.m, probs=c(0.975))), "dist.m"] <- quantile(veget$dist.m, probs=c(0.975))
hist(veget$dist.m)


## ----- Estimate site-level tree density and basal area -----

# Based on: Mitchell (2015) Quantitative Analysis by the Point-Centered Quarter Method

trees <- tibble(Camera.Trap.Name=sort(unique(veget$Camera.Trap.Name)), mean.distance=as.numeric(NA))

# mean.distance is the sum of the nearest point-to-tree distances in
# the quarters surveyed divided by the number of quarters
for(i in 1:nrow(trees)) {
  df1 <- subset(veget, Camera.Trap.Name == trees[[i,1]])
  sum.point.to.tree.distances <- sum(df1$dist.m, na.rm=TRUE)
  number.of.quarters <- length(unique(df1$id.ponto))
  trees[[i,2]] <- sum.point.to.tree.distances/number.of.quarters
}
trees # check

# absolute density is lambda=1/(r^2):
trees$tree.density <- as.numeric(NA)
for(i in 1:nrow(trees)) {
  trees[[i,3]] <- 1/(trees[[i,2]]^2) * 10000 # per ha
}
trees # check


# basal area: determine the total cover or basal area of the trees in the sample
# by species, and then calculate the mean basal area for each species
# (in our case there is no species so all species go together)
trees$basal.area <- as.numeric(NA)
for(i in 1:nrow(trees)) {
  df1 <- subset(veget, Camera.Trap.Name == trees[[i,1]])
  #mean.basal.area <- mean(df1$ba, na.rm=TRUE)
  mean.basal.area <- mean(df1$ba[-which(df1[,"ba"] == max(df1$ba) )], na.rm=TRUE) # basal area without largest tree
  trees[[i,4]] <- trees[[i,3]]*mean.basal.area
}
trees # check

hist(trees$tree.density, xlab="Tree density (trees/ha)", main="")
hist(trees$basal.area, xlab="Basal area (m²2/ha)", main="")

write.csv(trees, here("data", "trees.csv"), row.names=FALSE)
