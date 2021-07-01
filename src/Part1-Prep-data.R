# Part1-Prepare data

# Read in a TEAM data set and create and format so it is ready for wildlife community model
# Written by Jorge Ahumada @ Conservation International
# Adapted by Elildo Carvalho Jr @ ICMBio/CENAP, 2020-04-02

#----- 1 - Load libraries-----
library(dplyr)
library(lubridate)
library(here)


#----- 2 - Source files-----
here <- here::here # to avoid confusion with "here" function from lubridate
#source(here("bin", "camera trap analysis functions-10-06-18.R")) # using package here to build a path to the subdirectory "bin"
source(here("bin", "ahumada_codes.R"))
source(here("bin", "f-matrix-creator-experimental-probably-ok-but-need-check.R"))
source(here("src", "fix_species_names.R")) # fix some names and remove "false" species


## ----Load data-------
#dataRBG <- f.readin.fix.data(here("data", "Wild_ID_RBG_2016to2019.csv"))
dataRBG <- read.csv(here("data", "Wild_ID_RBG_2016to2019.csv"))

# some fixes
dataRBG$Sampling.Unit.Name <- as.factor(dataRBG$Camera.Trap.Name)
colnames(dataRBG)[9] <- "Photo.Time"
dataRBG$bin <- factor(dataRBG$bin)

# fix date formats (only needed if data was read with read.csv instead of f.readin.fix.data)
dataRBG$Photo.Date <- as.Date(dataRBG$Photo.Date)
dataRBG$Camera.Start.Date <- as.Date(dataRBG$Camera.Start.Date)
dataRBG$Camera.End.Date <- as.Date(dataRBG$Camera.End.Date)
dataRBG$Start.Date <- as.Date(dataRBG$Start.Date)
dataRBG$End.Date <- as.Date(dataRBG$End.Date)


# fix species names
f.fix.species.names(dataRBG)
dataRBG <- dataTemp # use new df created by function

#----- 4 - Extract binary presence/absence matrices for each species
species <- unique(dataRBG$bin)
cams <- unique(dataRBG$Camera.Trap.Name)
years <- unique(dataRBG$Sampling.Event)
secondPeriods <- 1:10

# Separate different years - clunky?
dataRBG2016 <- dplyr::filter(dataRBG, Sampling.Event == 2016)
dataRBG2017 <- dplyr::filter(dataRBG, Sampling.Event == 2017)
dataRBG2018 <- dplyr::filter(dataRBG, Sampling.Event == 2018)
dataRBG2019 <- dplyr::filter(dataRBG, Sampling.Event == 2019)


# use only first 50 days of sampling
# since we are using only 1st 50 days, we must reset end dates using max photo date
f.update.end.data <- function(data, duration){
  new.end.date <- min(data$Start.Date)+duration
  df1 <- subset(data, Photo.Date <= new.end.date)
  for(i in 1:nrow(df1)){ 
    if (df1$End.Date[i] > new.end.date) {
      df1$End.Date[i] <- new.end.date
    }
  }
  df1$Camera.Start.Date <- df1$Start.Date
  df1$Camera.End.Date <- df1$End.Date
  assign("df1", df1, envir=.GlobalEnv)
} # End of function

f.update.end.data(dataRBG2016, 50)
dataRBG2016 <- df1
f.update.end.data(dataRBG2017, 50)
dataRBG2017 <- df1
f.update.end.data(dataRBG2018, 50)
dataRBG2018 <- df1
f.update.end.data(dataRBG2019, 50)
dataRBG2019 <- df1


# Create presence/absence matrices for each species each year
# matrix dimensions are all identical accross species and years

# before using f.matrix.creator check sampling duration
duration <- function(data) {
  sampling.days <- max(data$End.Date) - min(data$Start.Date) + 1
  return(sampling.days)
}

duration(dataRBG2016) # 
#round(55/5) # get the number of occasions argument for f.matrix.creator4
duration(dataRBG2017)
duration(dataRBG2018)
duration(dataRBG2019)

paMats2016 <- f.matrix.creator4(dataRBG2016, species, 10)
paMats2017 <- f.matrix.creator4(dataRBG2017, species, 10)
paMats2018 <- f.matrix.creator4(dataRBG2018, species, 10)
paMats2019 <- f.matrix.creator4(dataRBG2019, species, 10)

dim(paMats2016[[1]]) # check
paMats2016[[1]] # check

# check species names
names(paMats2016) # Psophia obscura is the 1st species

# function to create species data
createSppData <- function(x) {
  for(i in 1:length(x)){
    df1 <- as.data.frame(paMats2016[x])
    colnames(df1) <- seq(1:length(colnames(df1))); colnames(df1) <- paste("X2016.", colnames(df1), sep="")
    df2 <- as.data.frame(paMats2017[x]) 
    colnames(df2) <- seq(1:length(colnames(df2))); colnames(df2) <- paste("X2017.", colnames(df2), sep="")
    df3 <- as.data.frame(paMats2018[x]) 
    colnames(df3) <- seq(1:length(colnames(df3))); colnames(df3) <- paste("X2018.", colnames(df3), sep="")
    df4 <- as.data.frame(paMats2019[x]) 
    colnames(df4) <- seq(1:length(colnames(df4))); colnames(df4) <- paste("X2019.", colnames(df4), sep="")
    bla <- cbind(df1, df2, df3, df4)
  }
  assign(paste("dataRBG_species", gsub(" ", "_", x), sep="_"), bla, envir = .GlobalEnv)
}

# check if it works
createSppData("Psophia obscura")
dataRBG_species_Psophia_obscura
dim(dataRBG_species_Psophia_obscura)


#----- 4 - Read covariate data

# Land cover Mapbiomas
cover <- read.csv(here("data", "cover_mapbiomas.csv"))
names(cover)[2] <- "Camera.Trap.Name"
names(cover)[4] <- "cover"
cover$Camera.Trap.Name <- gsub("Ctrbg", "CT-RBG-", cover$Camera.Trap.Name)
cover$Camera.Trap.Name <- gsub("Ctrgb", "CT-RBG-", cover$Camera.Trap.Name)
cover <- cover[,c("Camera.Trap.Name", "cover")]
head(cover)
hist(cover$cover)
sort(cover$cover)
# virtually all sites have 100% forest cover so maybe this variable should not be used
# an alternative would be to use distance to edges

# Distance to water
dist.water <- read.csv(here("data", "dist_agua_conv_trsh6_13fev.csv"))
names(dist.water)[1] <- "Camera.Trap.Name"
dist.water$Camera.Trap.Name <- gsub("Ctrbg", "CT-RBG-", dist.water$Camera.Trap.Name)
dist.water$Camera.Trap.Name <- gsub("Ctrgb", "CT-RBG-", dist.water$Camera.Trap.Name)
names(dist.water)[4] <- "dist.water"
dist.water <- dist.water[,c("Camera.Trap.Name", "dist.water")]
head(dist.water)

# Distance to forest edge
dist.edge <- read.csv(here("data", "dist_to_edge.csv"))
names(dist.edge) <- c("Camera.Trap.Name", "dist.edge")
dist.edge$Camera.Trap.Name <- gsub("Ctrbg", "CT-RBG-", dist.edge$Camera.Trap.Name)
dist.edge$Camera.Trap.Name <- gsub("Ctrgb", "CT-RBG-", dist.edge$Camera.Trap.Name)
head(dist.edge)

# distance to pasture > 10ha is exactly the same as distance to edge so lets keep the former
#dist.pasto <- read.csv(here("data", "dist_pasto10ha.csv"))
#names(dist.pasto)[1] <- "Camera.Trap.Name"
#dist.pasto$Camera.Trap.Name <- gsub("Ctrbg", "CT-RBG-", dist.pasto$Camera.Trap.Name)
#dist.pasto$Camera.Trap.Name <- gsub("Ctrgb", "CT-RBG-", dist.pasto$Camera.Trap.Name)
#names(dist.pasto)[2] <- "dist.pasto"
#dist.pasto <- dist.pasto[,c("Camera.Trap.Name", "dist.pasto")]
#head(dist.pasto)


# elevation
elev <- read.csv(here("data", "slope_elev.csv"))
names(elev)[1] <- "Camera.Trap.Name"
names(elev)[3] <- "elevation"
elev$Camera.Trap.Name <- gsub("Ctrbg", "CT-RBG-", elev$Camera.Trap.Name)
elev$Camera.Trap.Name <- gsub("Ctrgb", "CT-RBG-", elev$Camera.Trap.Name)
elev <- elev[,c("Camera.Trap.Name", "elevation")]
head(elev)

# tree structure
trees <- read.csv(here("data", "trees.csv"))
head(trees)
trees <- trees[,c("Camera.Trap.Name", "basal.area", "density")]

# camera array (block)
block <- read.csv(here("data", "blocos.csv"))
names(block)[1] <- "Camera.Trap.Name"
names(block)[2] <- "block"
block$Camera.Trap.Name <- gsub("Ctrbg", "CT-RBG-", block$Camera.Trap.Name)
block$Camera.Trap.Name <- gsub("Ctrgb", "CT-RBG-", block$Camera.Trap.Name)
block <- block[,c("Camera.Trap.Name", "block")]
head(block)


## create a single covariates dataframe
covars <- merge(cover, dist.water, by="Camera.Trap.Name")
covars <- merge(covars, dist.edge, by="Camera.Trap.Name")
covars <- merge(covars, elev, by="Camera.Trap.Name")
covars <- merge(covars, trees[,c("Camera.Trap.Name", "basal.area", "density")], by="Camera.Trap.Name")
covars <- merge(covars, block, by="Camera.Trap.Name")
head(covars)

# merge 
dataRBG_species_Psophia_obscura$Camera.Trap.Name <- rownames(dataRBG_species_Psophia_obscura)
row.names(dataRBG_species_Psophia_obscura) <- NULL
pobscura <- dataRBG_species_Psophia_obscura
head(pobscura)


# create a separate 2017 dataset for single-season model
#Pobscura2017 <- dataRBG_species_Psophia_obscura[,c(45,12:22)]
#Pobscura2017 <- merge(Pobscura2017, covars, by="Camera.Trap.Name")

pobscura <- merge(pobscura, covars, by="Camera.Trap.Name")
names(pobscura)[1] <- "cams"

# Save to disk
#saveRDS(Pobscura2017, here("data","Pobscura2017.rds"))
saveRDS(pobscura, here("data","pobscura.rds"))
