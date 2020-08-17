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
#source(here("src", "fix_species_names.R")) # fix some names and remove "false" species


## ----Load data-------
#dataRBG <- f.readin.fix.data(here("data", "Wild_ID_RBG_2016to2019.csv"))
dataRBG <- read.csv(here("data", "Wild_ID_RBG_2016to2019.csv"))

# some fixes
dataRBG$Camera.Trap.Name <- as.factor(dataRBG$Camera.Trap.Name)
dataRBG$Sampling.Unit.Name <- as.factor(dataRBG$Camera.Trap.Name)
colnames(dataRBG)[9] <- "Photo.Time"
dataRBG$bin <- factor(dataRBG$bin)

# fix date formats (only needed if data was read with read.csv instead of f.readin.fix.data)
dataRBG$Photo.Date <- as.Date(dataRBG$Photo.Date)
dataRBG$Camera.Start.Date <- as.Date(dataRBG$Camera.Start.Date)
dataRBG$Camera.End.Date <- as.Date(dataRBG$Camera.End.Date)
dataRBG$Start.Date <- as.Date(dataRBG$Start.Date)
dataRBG$End.Date <- as.Date(dataRBG$End.Date)

# use only first 50 days of sampling
#dataRBG <- dataRBG[dataRBG$Photo.Date <= (dataRBG$Start.Date+50),]

# since we are using only 1st 55 days, we must reset end dates using max photo date
f.update.end.data <- function(data, duration){
  for(j in 1: length(unique(data$Sampling.Event))) {
    df <- subset(data, Sampling.Event == data$Sampling.Event[j])
    for(i in 1:nrow(df)){ 
    if (df$End.Date[i] >= min(df$Start.Date)+duration) {
      df$End.Date[i] <- min(df1$Start.Date)+duration 
    }
  }
  df$Camera.End.Date <- df$End.Date
  assign("dataTemp1", df1, envir=.GlobalEnv)
} }# End of function

f.update.end.data(dataRBG, 55) # FUNCAO TEM ERROS CORRIGIR PAREI NO MEIO
dataRBG <- dataTemp1

# fix species names
#f.fix.species.names(dataRBG)
#dataRBG <- dataTemp # use new df created by function

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


# Create presence/absence matrices for each species each year
# matrix dimensions are all identical accross species and years

#paMats2016 <- f.matrix.creator3(dataRBG2016, cams, species)
#paMats2017 <- f.matrix.creator3(dataRBG2017, cams, species)
#paMats2018 <- f.matrix.creator3(dataRBG2018, cams, species)
#paMats2019 <- f.matrix.creator3(dataRBG2019, cams, species)

# before using f.matrix.creator check sampling duration
duration <- function(data) {
  sampling.days <- max(data$End.Date) - min(data$Start.Date) + 1
  return(sampling.days)
}

duration(dataRBG2016) # 55 days, if using 5 days as one occasion detection history will need 11 occasions
duration(dataRBG2017) # 79 days will need 15 occasions
duration(dataRBG2018) # 58 days will need 11 occasions
duration(dataRBG2019) # 56 days will need 11 occasions

paMats2016 <- f.matrix.creator4(dataRBG2016, "Psophia obscura", 12)
paMats2017 <- f.matrix.creator4(dataRBG2017, "Psophia obscura", 15)
paMats2018 <- f.matrix.creator4(dataRBG2018, "Psophia obscura", 11)
paMats2019 <- f.matrix.creator4(dataRBG2019, "Psophia obscura", 11)


#paMats2016 <- f.matrix.creator4(dataRBG2016, species, 12)
#paMats2017 <- f.matrix.creator4(dataRBG2017, species, 15)
#paMats2018 <- f.matrix.creator4(dataRBG2018, species, 11)
#paMats2019 <- f.matrix.creator4(dataRBG2019, species, 11)


# ---- WARNING: ---- 
# colext function requires an equal number of secondary periods (J) per primary periods (T)
# however, in this dataset J varies between years 
# so you must add additional Js with NAs so that all paMats have the same dimensions
# for now just an experiment, ideally this must be included in the createSppData function
newColsNumber16 <- dim(paMats2019[x][[1]])[2] - dim(paMats2016[x][[1]])[2]
newCols16 <- matrix(rep(NA, newColsNumber16*dim(paMats2016[x][[1]])[1]), ncol=newColsNumber16)
newColsNumber17 <- dim(paMats2019[x][[1]])[2] - dim(paMats2017[x][[1]])[2]
newCols17 <- matrix(rep(NA, newColsNumber17*dim(paMats2017[x][[1]])[1]), ncol=newColsNumber17)
newColsNumber18 <- dim(paMats2019[x][[1]])[2] - dim(paMats2018[x][[1]])[2]
newCols18 <- matrix(rep(NA, newColsNumber18*dim(paMats2018[x][[1]])[1]), ncol=newColsNumber18)

# add new cols to paMats
for(i in 1:length(paMats2016)) {
  paMats2016[[i]] <- cbind(paMats2016[[i]], newCols16)
}
for(i in 1:length(paMats2017)) {
  paMats2017[[i]] <- cbind(paMats2017[[i]], newCols17)
}
for(i in 1:length(paMats2018)) {
  paMats2018[[i]] <- cbind(paMats2018[[i]], newCols18)
}

#----- 4 - create species data similar to template from Beaudrot repo

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


#---------------------------------------



#Combine all in one list
paMatsdataRBG <- list(paMats2016, paMats2017, paMats2018, paMats2019)

# The function createArragyG puts the data in a format (list) that will be ready for
# ingestion by a jags model
arrayG_RBG <- createArrayG(data = paMatsdataRBG, species = species, cams = cams, years = years, secondPeriods = secondPeriods, siteName = "dataRBG")

# Save to disk
saveRDS(arrayG_RBG, here("data","arrayG_RBG.rds"))
