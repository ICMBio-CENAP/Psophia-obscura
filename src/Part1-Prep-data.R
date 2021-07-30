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
source(here("bin", "fix_species_names.R")) # fix some names and remove "false" species


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
dataRBG$Photo.Time <- format(as.POSIXct(dataRBG$Photo.Time, format="%HH %MM %SS"), "%H:%M:%S")
dataRBG$td.photo <- as.POSIXct(paste(dataRBG$Photo.Date, dataRBG$Photo.Time, sep=" "), format="%Y-%m-%d %H:%M:%OS", tz="UTC")


# fix species names
f.fix.species.names(dataRBG)
dataRBG <- dataTemp # use new df created by function

# calculate sampling effort, number of records etc
tempdf <- distinct(dataRBG, Camera.Trap.Name, Sampling.Event, Start.Date, End.Date)
tempdf$effort <- as.numeric(tempdf$End.Date-tempdf$Start.Date)
head(tempdf)
hist(tempdf$effort)

# Sampling was much longer in 2017 than in other years
# let limit maximum sampling to 60 days
# for this we must reset end dates using max photo date
f.update.end.data <- function(data){
  for(i in 1:nrow(data)){
    if (data$End.Date[i] > data$Start.Date[i] + 50) {
      data$End.Date[i] <- data$Start.Date[i] + 50
    }
  }
  data$Camera.Start.Date <- data$Start.Date
  data$Camera.End.Date <- data$End.Date
  .GlobalEnv$temp_data <- data
} # End of function

f.update.end.data(dataRBG)
# check:
tempdf <- distinct(temp_data, Camera.Trap.Name, Sampling.Event, Start.Date, End.Date)
tempdf$effort <- as.numeric(tempdf$End.Date-tempdf$Start.Date)
head(tempdf)
hist(tempdf$effort)

# overwrite dataRBG with updated temp_data
dataRBG <- temp_data

# calculate effort number of records etc
sum(tempdf$effort) # total effort across years
tempdf %>%
  group_by(Sampling.Event) %>%
  dplyr::summarize(mean = mean(effort), total=sum(effort) )

# number of photos
dataRBG %>%
  group_by(Sampling.Event) %>%
  filter(bin == "Psophia obscura") %>%
  summarize(photos = n())

# number of independent records
temp1 <- dataRBG %>%
  filter(bin == "Psophia obscura")
temp1$td.photo <- as.POSIXct(temp1$td.photo) # somehow td.photo class was reverted to factor 
temp1 <- f.separate.events(temp1, 60)
temp1 <- distinct(temp1, Camera.Trap.Name, bin, grp, .keep_all = TRUE)
temp1 %>%
  group_by(Sampling.Event) %>%
  filter(bin == "Psophia obscura") %>%
  summarize(records = n())

temp1 %>%
  group_by(Sampling.Event) %>%
  filter(bin == "Psophia obscura") %>%
  summarize(mean_group = mean(Number.of.Animals), min = min(Number.of.Animals), max = max(Number.of.Animals))

temp2 <- temp1 %>%
  group_by(Sampling.Event) %>%
  filter(bin == "Psophia obscura")
hist(as.numeric(temp2$Number.of.Animals))


#----- 4 - Extract binary presence/absence matrices for each species
species <- unique(dataRBG$bin)
cams <- unique(dataRBG$Camera.Trap.Name)
years <- unique(dataRBG$Sampling.Event)
secondPeriods <- 1:10

# Separate different years
dataRBG2016 <- dplyr::filter(dataRBG, Sampling.Event == 2016)
dataRBG2017 <- dplyr::filter(dataRBG, Sampling.Event == 2017)
dataRBG2018 <- dplyr::filter(dataRBG, Sampling.Event == 2018)
dataRBG2019 <- dplyr::filter(dataRBG, Sampling.Event == 2019)


# use only first 50 days of sampling
# since we are using only 1st 50 days, we must reset end dates using max photo date
# NB! NOT NEEDED AS WE ALREADY UPDATED END DATES
#f.update.end.data <- function(data, duration){
#  new.end.date <- min(data$Start.Date)+duration
#  df1 <- subset(data, Photo.Date <= new.end.date)
#  for(i in 1:nrow(df1)){ 
#    if (df1$End.Date[i] > new.end.date) {
#      df1$End.Date[i] <- new.end.date
#    }
#  }
#  df1$Camera.Start.Date <- df1$Start.Date
#  df1$Camera.End.Date <- df1$End.Date
#  assign("df1", df1, envir=.GlobalEnv)
#} # End of function

#f.update.end.data(dataRBG2016, 50)
#dataRBG2016 <- df1
#f.update.end.data(dataRBG2017, 50)
#dataRBG2017 <- df1
#f.update.end.data(dataRBG2018, 50)
#dataRBG2018 <- df1
#f.update.end.data(dataRBG2019, 50)
#dataRBG2019 <- df1


# Create presence/absence matrices for each species each year
# matrix dimensions are all identical accross species and years

# before using f.matrix.creator check sampling duration
duration <- function(data) {
  sampling.days <- max(data$End.Date) - min(data$Start.Date) + 1
  return(sampling.days)
}

duration(dataRBG2016) #
round(54/5) # get the number of occasions argument for f.matrix.creator4
#round(55/10)
duration(dataRBG2017)
round(60/5)
#round(70/10)
duration(dataRBG2018)
round(57/5)
#round(58/10)
duration(dataRBG2019)
round(56/5)
#round(56/10)

paMats2016 <- f.matrix.creator4(dataRBG2016, species, 11) # 11 if using 5-day occasion
paMats2017 <- f.matrix.creator4(dataRBG2017, species, 12) # 14
paMats2018 <- f.matrix.creator4(dataRBG2018, species, 11) # 12
paMats2019 <- f.matrix.creator4(dataRBG2019, species, 11) # 11

dim(paMats2016[[1]]) # check
paMats2016[[1]] # check

# check species names
names(paMats2016) # Psophia obscura is the 1st species

# matrices from different years should have the same size
# so we have to add NA only columns to some matrices
# check here how many new cols will be needed and add using createSppData function
dim(paMats2016[[1]]) # check
dim(paMats2017[[1]]) # this matrix is the one with more cols, the others must match it 
dim(paMats2018[[1]])
dim(paMats2019[[1]])


# function to create species data
createSppData <- function(x) {
  for(i in 1:length(x)){
    #df1 <- as.data.frame(paMats2016[x])
    df1 <- as.data.frame(cbind(paMats2016[[x]], matrix(NA, 61, 1))) # matrix(NA, 61, 3)) if using 5-day occasion
    colnames(df1) <- seq(1:length(colnames(df1))); colnames(df1) <- paste("X2016.", colnames(df1), sep="")
    df2 <- as.data.frame(paMats2017[x]) 
    colnames(df2) <- seq(1:length(colnames(df2))); colnames(df2) <- paste("X2017.", colnames(df2), sep="")
    #df3 <- as.data.frame(paMats2018[x])
    df3 <- as.data.frame(cbind(paMats2018[[x]], matrix(NA, 61, 1))) # 61,3...
    colnames(df3) <- seq(1:length(colnames(df3))); colnames(df3) <- paste("X2018.", colnames(df3), sep="")
    #df4 <- as.data.frame(paMats2019[x])
    df4 <- as.data.frame(cbind(paMats2019[[x]], matrix(NA, 61, 1))) # 61,3...
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
names(elev)[3] <- "point.elevation"
names(elev)[6] <- "range.elevation"
elev$Camera.Trap.Name <- gsub("Ctrbg", "CT-RBG-", elev$Camera.Trap.Name)
elev$Camera.Trap.Name <- gsub("Ctrgb", "CT-RBG-", elev$Camera.Trap.Name)
elev <- elev[,c("Camera.Trap.Name", "point.elevation", "range.elevation")]
head(elev)

# tree structure
trees <- read.csv(here("data", "trees.csv"))
head(trees)
trees <- trees[,c("Camera.Trap.Name", "basal.area", "tree.density")]

# camera array (block)
block <- read.csv(here("data", "blocos.csv"))
names(block)[1] <- "Camera.Trap.Name"
names(block)[2] <- "block"
block$Camera.Trap.Name <- gsub("Ctrbg", "CT-RBG-", block$Camera.Trap.Name)
block$Camera.Trap.Name <- gsub("Ctrgb", "CT-RBG-", block$Camera.Trap.Name)
block <- block[,c("Camera.Trap.Name", "block")]
head(block)

# recovery time
# use data from the Forest Ecology and Management paper
covars_hmsc <- read.csv(here("data", "covars_hmsc.csv"))
covars_hmsc$recovery <- as.numeric(covars_hmsc$recovery)
recovery <- covars_hmsc[,c("Camera.Trap.Name", "recovery")]
head(recovery)

## create a single covariates dataframe
covars <- merge(cover, dist.water, by="Camera.Trap.Name", all.x=TRUE)
covars <- merge(covars, dist.edge, by="Camera.Trap.Name", all.x=TRUE)
covars <- merge(covars, elev, by="Camera.Trap.Name", all.x=TRUE)
covars <- merge(covars, trees, by="Camera.Trap.Name", all.x=TRUE)
covars <- merge(covars, block, by="Camera.Trap.Name", all.x=TRUE)
covars <- merge(covars, recovery, by="Camera.Trap.Name", all.x=TRUE)
head(covars)
dim(covars)
covars[which(is.na(covars$recovery)),]
# there are two NAs in recovery, the fix is the following
covars[covars$Camera.Trap.Name=="CT-RBG-1-11", "recovery"] <- 6
covars[covars$Camera.Trap.Name=="CT-RBG-2-82", "recovery"] <- 15
covars[which(is.na(covars$recovery)),]

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


bla1 <- data.frame(cams)
names(bla1) <- "Camera.Trap.Name"
bla3 <- merge(bla1, dataRBG[,c("Camera.Trap.Name", "Latitude", "Longitude")], by="Camera.Trap.Name")
bla3
bla4 <- distinct(bla3)
bla4
write.csv(bla4, here("data", "psophia-cams.csv"))
