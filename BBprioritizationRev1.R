# code for revising bumble bee conservation areas for Canada under current and future climate scenarios 

#First section is making SDMs for each NA bumble bee species under current and future climate
#Prioritization analyses will be in the second half.

#SDM packages
library(raster) #for accessing raster data
library(tidyverse) #for combining outputs, data manipulation and cleaning
library(usdm) #vif function for checking collinearity
library(ENMeval)

#working directory for SDMS
setwd("C:\\Users\\Amanda\\Documents\\PhD\\Chapter 3\\Manuscript\\Revision 1\\Analyses\\Current") #when working on current climate
setwd("C:\\Users\\Amanda\\Documents\\PhD\\Chapter 3\\Manuscript\\Revision 1\\Analyses\\\Future2.6SDMs") #when working on future climate rcp 2.6
setwd("C:\\Users\\Amanda\\Documents\\PhD\\Chapter 3\\Manuscript\\Revision 1\\Analyses\\Future8.5SDMs") #when working on future climate rcp 8.5

#Data for SDMS
bombus<-read.csv("C:\\Users\\Amanda\\Documents\\PhD\\Chapter 3\\Manuscript\\Revision 1\\Analyses\\Data\\NABumblebeeRecordsClean.csv")#bumble bee occurrences. Code for prep below
coordinates(bombus)<-~longitude+latitude # assigns the occurrences to spatial points data frame
proj4string(bombus)<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" #sets the CRS for occurrence points
polyNA<-readOGR(dsn= "C:\\Users\\Amanda\\Documents\\PhD\\Chapter 3\\North America Datafiles\\NApolygon.shp", layer="NApolygon") #North America shapefile from rgdal
biasfile<-raster("C:\\Users\\Amanda\\Documents\\PhD\\Chapter 3\\North America Datafiles\\CurrentSDMs\\bias.raster.tif") #biasfile made from all bumble bees and all years. Code for prep below
current<-stack("C:\\Users\\Amanda\\Documents\\PhD\\Chapter 3\\North America Datafiles\\CurrentSDMs\\currentclimateNA2.5min.grd") # current climate raster from worldlclim at 2.5 km resolution. Code for prep below

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
allbumble<-read.csv("C:\\Users\\Amanda\\Documents\\PhD\\Chapter 3\\North America Datafiles\\all.shareable.bbna.03.10.2020.csv") #all bumble bee occurrences for Canada for all years
allbumble<-allbumble %>%
  filter(latitude != "") %>%
  filter(longitude != "") %>%
  filter(longitude > extent(polyNA)[1] & longitude < extent(polyNA)[2]) %>% filter(latitude > extent(polyNA)[3]  & latitude < extent(polyNA)[4]) ## crop out to polyNA
coordinates(allbumble)<-~longitude+latitude #assigning coordinates columns 
proj4string(allbumble)<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" #defining crs
biasocc<- allbumble
pres.locs <- data.frame(coordinates(biasocc))
xband <- MASS::bandwidth.nrd(pres.locs[,1]) ## Find bandwidths for x of raster
yband <- MASS::bandwidth.nrd(pres.locs[,2]) ## Find bandwidths for y of raster
dens <- KernSmooth::bkde2D(cbind(pres.locs[,1], pres.locs[,2]), bandwidth=c(xband,yband), gridsize=c(nrow(current),ncol(current))) ## density function based on raster cells to create bias raster
dens.ras <- raster(list(x=dens$x1,y=dens$x2,z=dens$fhat))
crs(dens.ras) <- crs(current)
dens.ras <- resample(dens.ras, current, method="bilinear") #make bias raster in the same resolution of landcover file
writeRaster(dens.ras, "bias.raster.tif")