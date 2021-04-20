# code for revising bumble bee conservation areas for Canada under current and future climate scenarios 

#First section is making SDMs for each NA bumble bee species under current and future climate
#Prioritization analyses will be in the second half.

#SDM packages
library(raster) #for accessing raster data
library(tidyverse) #for combining outputs, data manipulation and cleaning
library(usdm) #vif function for checking collinearity
library(ENMeval) #function for MaxEnt to check feature class settings and regularization settings
library(dismo) #used to add background points to the SDM

## specify temporary directory for rasters
#rasterOptions(tmpdir= "D://Documents//Temp")

#Data for SDMS
bombus<-read.csv("NABumblebeeRecordsClean.csv")#bumble bee occurrences. Code for prep below
coordinates(bombus)<-~longitude+latitude # assigns the occurrences to spatial points data frame
proj4string(bombus)<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" #sets the CRS for occurrence points
bias<-raster("bias.raster.tif")
poly<-shapefile("NApolygon.shp")
current<-stack("Current.tif")

### Current Climate Analyses

## testing new package using one species for now - affinis
## checking for and removing collinear variables for affinis

affinis<-subset(bombus, species=="affinis")

random <- sampleRandom(current, 7070, xy=T) #adding random background points for comparison. 10 x sample size 
random <- rbind(data.frame(longitude=random[,1],latitude=random[,2]), coordinates(affinis)) ## add occurrences to background points
ext.ppt <- raster::extract(current, random) ## extract bioclim variables for those points
vifOut <- vifcor(ext.ppt) 
climateaff<- subset(current, as.character(vifOut@results$Variables)) #climate dataset without collinear var

# trying ENMeaval package

affinis<-as.data.frame(affinis)#needs occurrence points as a data frame 
affinis2<-affinis[,c(4,3)] #can only have two columns, longitude and latitude in that order

set.seed(0)
affBkg <- xyFromCell(bias, sample(ncell(bias), 10000, prob=values(bias))) #background points
regs<-c(0.5,1,2,3,4,5) #regularization multipliers
feats<-c("L", "Q", "P", "LQ", "H", "HQ", "T", "HQP", "HQPT", "HPLQ") #feature classes

aff<-ENMevaluate(affinis2, climateaff, bg.coords=affBkg, RMvalues= regs, fc=feats, 
                 method= "randomkfold", kfold=5, algorithm="maxent.jar")
predictOut <- predict( climateaff, aff@models[[39]], type="logistic",progress="text") ## predictOut raster
#function will need to normalize rasters after accounting for biasfile
#may need to assign CRS at some point

## find threshold, reviewer mentioned this but I don't think I'm going to do it
evalMaxent <- evaluate(affinis2, affBkg, aff@models[[1]], climateaff) ## threshold @ max TPR+TNR
threshold(evalMaxent)

writeRaster(predictOutCorrected, "affinisOut.grd", overwrite=T)

####################
###Climate Data Prepping

##Current climate
currentclim<-list.files("current", full.names = TRUE, pattern=".bil")
currentclim.all<-stack(currentclim)
currentclim.all <- crop(currentclim.all, poly) 
currentclim.all <- mask(currentclim.all, poly)
writeRaster(currentclim.all, "Current.tif")

##Prepping future climate data from worldclim

#RCP 2.6
bioclim <- list.files("RCP2.6", full.names = TRUE, pattern=".tif")
bioclim.all <- stack(bioclim)
bioclim.all <- crop(bioclim.all, poly) 
bioclim.all <- mask(bioclim.all, poly) 

## add letter to the end of the bioclim vars to separate out the 1s from the 10s
rasterNames <- paste0(names(bioclim.all),"e")

## Iterate through all 19 bioclim variables
future26 <- stack()
future26 <- lapply(1:19, function(i){
biovar <- paste0("bi50",i,"e") ## select bioclimate variable per iteration
bioclimTemp <- bioclim.all[[grep(biovar,rasterNames)]] ## select that climate variable
meanTemp <- mean(bioclimTemp) ## average across GCMs
future26 <- stack(future26, meanTemp)
})
## convert list of rasters to raster stack
all26 <- do.call(stack, future26)
names(all26) <- paste0(rep("bio",19),1:19)
writeRaster(all26, "Future2.6.grd")

## NOT RUN - prepped during initial manuscript submission and does not need to be revised. 

##Prepping bumble bee occurrence data
#allbumble<-read.csv("C:\\Users\\Amanda\\Documents\\PhD\\Chapter 3\\North America Datafiles\\all.shareable.bbna.03.10.2020.csv") 
#filter(species != "bombus") %>% #allbumble has 625485 observations
  #data <-allbumble%>%
  #filter(year > 2007 & year < 2019) 
#nrow(data) #227735 observations after subsetting to 2008-2018
#length(unique(data$species)) #47 species 
#data <- data %>% 
 # filter(species != "") %>%
  #filter(latitude != "") %>%
  #filter(longitude != "")
#nrow(data)#225130 observations after removing blank data cells
#length(unique(data$species)) #47 species 
#data <- data %>% distinct(species,latitude,longitude) %>%
#  droplevels()
#nrow(data) #70184 observations after removing duplicate records for the same species at the same location
#length(unique(data$species)) #47 species
#data %>% group_by(species) %>% summarize(n=length(species)) %>% arrange(n) %>% data.frame()
#species with fewer than ten records will be dropped
#cockerelli 5, distinguendus 6, variabilis  9
#data<-data %>%
  #filter(!species%in%c("cockerelli","distinguendus","variabilis"))
#write.csv(data, "NABumblebeeRecordsClean.csv")

#Prepping bioclim current climate data 
#bioclim <- list.files("C:\\Users\\Amanda\\Documents\\PhD\\Chapter 3\\North America Datafiles\\CurrentClimate2.5min", full.names = TRUE)
#bioclim.all <- stack(bioclim)
#names(bioclim.all) <- paste0(rep("bio",19),1:19) ## rename to logical names
#bioclim.all <- crop(bioclim.all, polyNA) 
#bioclim.all <- mask(bioclim.all, polyNA) 
#writeRaster(bioclim.all, "currentclimateNA2.5min.grd")

#creating a bias file
#allbumble<-read.csv("C:\\Users\\Amanda\\Documents\\PhD\\Chapter 3\\North America Datafiles\\all.shareable.bbna.03.10.2020.csv") #all bumble bee occurrences for Canada for all years
#allbumble<-allbumble %>%
 # filter(latitude != "") %>%
  #filter(longitude != "") %>%
  #filter(longitude > extent(polyNA)[1] & longitude < extent(polyNA)[2]) %>% filter(latitude > extent(polyNA)[3]  & latitude < extent(polyNA)[4]) ## crop out to polyNA
#coordinates(allbumble)<-~longitude+latitude #assigning coordinates columns 
#proj4string(allbumble)<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" #defining crs
#biasocc<- allbumble
#pres.locs <- data.frame(coordinates(biasocc))
#xband <- MASS::bandwidth.nrd(pres.locs[,1]) ## Find bandwidths for x of raster
#yband <- MASS::bandwidth.nrd(pres.locs[,2]) ## Find bandwidths for y of raster
#dens <- KernSmooth::bkde2D(cbind(pres.locs[,1], pres.locs[,2]), bandwidth=c(xband,yband), gridsize=c(nrow(current),ncol(current))) ## density function based on raster cells to create bias raster
#dens.ras <- raster(list(x=dens$x1,y=dens$x2,z=dens$fhat))
#crs(dens.ras) <- crs(current)
#dens.ras <- resample(dens.ras, current, method="bilinear") #make bias raster in the same resolution of landcover file
#writeRaster(dens.ras, "bias.raster.tif")

###

## prepping future climate data

