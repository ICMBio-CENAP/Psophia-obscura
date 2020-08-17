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

# use only first 50 days of sampling
dataRBG <- dataRBG[dataRBG$Photo.Date <= (dataRBG$Start.Date+50),]

# since we are using only 1st 40 days, we must reset end dates using max photo date
f.update.end.data <- function(data, duration){
  df1 <- data[data$Photo.Date <= (data$Start.Date+50),] # subset using desired duration
  for(i in 1:nrow(df1)){ 
    if (df1$End.Date[i] >= df1$Start.Date[i]+50) {
      df1$End.Date[i] <- df1$Start.Date[i]+50 
    }
  }
  df1$Camera.Start.Date <- df1$Start.Date
  df1$Camera.End.Date <- df1$End.Date
  assign("dataTemp1", df1, envir=.GlobalEnv)
} # End of function

f.update.end.data(dataRBG, 50)
dataRBG <- dataTemp1

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


# Create presence/absence matrices for each species each year
# matrix dimensions are all identical accross species and years

#paMats2016 <- f.matrix.creator3(dataRBG2016, cams, species)
#paMats2017 <- f.matrix.creator3(dataRBG2017, cams, species)
#paMats2018 <- f.matrix.creator3(dataRBG2018, cams, species)
#paMats2019 <- f.matrix.creator3(dataRBG2019, cams, species)

paMats2016 <- f.matrix.creator4(dataRBG2016, species, 10) # 10 is the number of occasions, can be changed
paMats2017 <- f.matrix.creator4(dataRBG2017, species, 10)
paMats2018 <- f.matrix.creator4(dataRBG2018, species, 10)
paMats2019 <- f.matrix.creator4(dataRBG2019, species, 10)


#Combine all in one list
paMatsdataRBG <- list(paMats2016, paMats2017, paMats2018, paMats2019)

# The function createArragyG puts the data in a format (list) that will be ready for
# ingestion by a jags model
arrayG_RBG <- createArrayG(data = paMatsdataRBG, species = species, cams = cams, years = years, secondPeriods = secondPeriods, siteName = "dataRBG")

# Save to disk
saveRDS(arrayG_RBG, here("data","arrayG_RBG.rds"))
