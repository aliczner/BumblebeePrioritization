# code for revising bumble bee conservation areas for Canada under current and future climate scenarios 

#First section is making SDMs for each NA bumble bee species under current and future climate
#Prioritization analyses will be in the second half.

#SDM packages
library(raster) #for accessing raster data
library(tidyverse) #for combining outputs, data manipulation and cleaning
library(usdm) #vif function for checking collinearity
library(ENMeval) #function for MaxEnt to check feature class settings and regularization settings
library(dismo) #used to add background points to the SDM
library(rgdal) #needed for the buffer cropping

## specify temporary directory for rasters
#rasterOptions(tmpdir= "D://Documents//Temp")

#Data for SDMS
bombus<-read.csv("NABumblebeeRecordsClean.csv")#bumble bee occurrences. Code for prep below
coordinates(bombus)<-~longitude+latitude # assigns the occurrences to spatial points data frame
proj4string(bombus)<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" #sets the CRS for occurrence points
bias<-raster("bias.raster.tif")
poly<-shapefile("NApolygon.shp")
current<-stack("Currentclim.tif")
future26<-stack("Future2.6.grd")
future85<-stack("Future8.5.grd")

### Jan 2023 ####

#occidentalis has been split to subspecies. B.o.occidentalis is south of 55 deg N and 
#B.o.mckayi is north of 55 deg N. Need to separate bombus dataset

bombus2 <- bombus %>% 
  mutate(species=ifelse(species=="occidentalis" & latitude > 55, "mckayi", species))
# 230 mckay and 788 occidentalis

## Rename raster stacks to be the same names
names(current) <- paste0(rep("bio",19),1:19)
names(future26) <- paste0(rep("bio",19),1:19)
names(future85) <- paste0(rep("bio",19),1:19)

#########################
### SDM Climate Analyses for current and future
########################

speciesList <- sort(as.character(unique(bombus$species)))

###making species distribution buffers 
#all bumble bee records
library(dplyr)

allbumble<-read.csv("D:\\Documents\\PhD\\Chapter 3\\North America Datafiles\\all.shareable.bbna.03.10.2020.csv")
allbumble %>%
  filter(allbumble$species != "bombus")
data <-allbumble %>%
  dplyr::filter(species != "") %>%
  dplyr::filter(latitude != "") %>%
  dplyr::filter(longitude != "")
data<-data %>%
  filter(!species%in%c("cockerelli","distinguendus","variabilis", "impatiens/bimaculatus", "franklini"))
historical<-data
historical.sp<-historical[,c("species", "latitude", "longitude")]
coordinates(historical.sp)<-~longitude+latitude # assigns the occurrences to spatial points data frame
proj4string(historical.sp)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

library(adehabitatHR) #minimum convex polygon package

histmcp<-mcp(historical.sp, percent=100) #function to calculate minimum convex polygon

for (i in 1:44){
  tempBuffer<-raster::buffer(subset(histmcp, id==speciesList[i]), width=2.70, dissolve=TRUE)
  tempBuffer<-SpatialPolygonsDataFrame(tempBuffer, data=data.frame(id=speciesList[i]), match.ID = F)
  rgdal::writeOGR(tempBuffer, dsn="./buffers", layer=paste0(speciesList[i], "buffer"), driver="ESRI Shapefile", overwrite_layer = T)
  print(i)
}

#need to add buffer for mckayi and a new occidentalis buffer
mckayi<- subset(bombus2, species=="mckayi")

historical2 <- historical %>% 
  mutate(species=ifelse(species=="occidentalis" & latitude > 55, "mckayi", species))
         
histmckayi<-subset(historical2, species=="mckayi")
histocc <- subset(historical2, species=="occidentalis")

histmckayi.sp<-histmckayi[,c("species", "latitude", "longitude")]
coordinates(histmckayi.sp)<-~longitude+latitude # assigns the occurrences to spatial points data frame
proj4string(histmckayi.sp)<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

histocc.sp <- histocc[,c ("species", "latitude", "longitude")]
coordinates(histocc.sp) <-~longitude+latitude # assigns the occurrences to spatial points data frame
proj4string(histocc.sp)<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

library(adehabitatHR)
histmckayi.mcp <- mcp(histmckayi.sp, percent=100)
histocc.mcp <- mcp(histocc.sp, percent=100)

mckayibuffer<-raster::buffer(subset(histmckayi.mcp), width=2.70, dissolve=TRUE)
occidentalisbuffer <- raster::buffer(subset(histocc.mcp, width=2.70, dissolve=TRUE))

sf::st_write(sf::st_as_sf(mckayibuffer), "mckayibuffer.shp")
sf::st_write(sf::st_as_sf(occidentalisbuffer), "occidentalisbuffer.shp")

rgdal::writeOGR(mckayibuffer, dsn=".", layer="mckayibuffer", driver ="ESRI Shapefile")
rgdal::writeOGR(occidentalisbuffer, dsn=".", layer="occidentalisbuffer", driver ="ESRI Shapefile")

#make pdfs of the outputs
speciesfiles<-list.files("Outputs//current",full.names=T)
for(i in 1:44){
  sppname<-gsub(".tif", "",speciesfiles[i])
  sppname<-gsub("Outputs//current/Current_", "",sppname)
  pdf(paste0("maps/",sppname,"Map.pdf"), useDingbats = F)
  plot(raster(speciesfiles[i]),main=sppname)
  points<-subset(bombus,species==sppname)
  plot(points,add=TRUE)
  plot(polyCAN, add=TRUE)
  dev.off()
}

#####MaxEnt using ENMeval package function

BBsdm<-function(spp, spDF, biasfile)  {
  ## Load current climate
  climate=current
  
  ## Data set up
  tempSpp <- subset(spDF, species == spp)
  background <- nrow(tempSpp)*10
  
  ## load buffer - least polygon extent for restricting SDM predictions
  sppBuffer <- readOGR(dsn="./buffers", layer=paste0(spp,"buffer"))
  
  ## Mask climate to extent
  climate <- mask(climate, sppBuffer)
  
  ## Checking for collinearity among the climate variables using vifcor
  random <- sampleRandom(climate, background, xy=T) #adding random background points for comparison. 10 x sample size 
  random <- rbind(data.frame(longitude=random[,1],latitude=random[,2]), coordinates(tempSpp)) ## add occurrences to background points
  ext.ppt <- raster::extract(climate, random) ## extract bioclim variables for those points
  vifOut <- vifcor(ext.ppt) 
  climateSp<- subset(climate, as.character(vifOut@results$Variables)) #climate dataset without collinear variables
  varsExcluded <- vifOut@excluded
  print(paste0("Collinearity checked, variables removed", varsExcluded))

  ## Maxent using ENMeaval package
  
  tempSpp<-as.data.frame(tempSpp)
  tempSpp2<-tempSpp[,c(4,3)]
  
  set.seed(1)
  SppBkg <- xyFromCell(bias, sample(ncell(bias), 10000, prob=values(bias))) #background points
  regs<-c(0.5,1,2,3,4,5) #regularization multipliers
  feats<-c("L", "Q", "P", "LQ", "H", "HQ", "T", "HQP", "HQPT", "HPLQ") #feature classes

  SppMax<-ENMevaluate(tempSpp2, climateSp, bg.coords=SppBkg, RMvalues= regs, fc=feats, 
                   method= "randomkfold", kfold=5, algorithm="maxent.jar")
  ## Determine best MaxEnt Model
  bestMax <- SppMax@results[which (SppMax@results$delta.AICc==0)[1],]
  bestModel <- SppMax@models[[which (SppMax@results$delta.AICc==0)[1]]]
  varImp <- var.importance(bestModel)
  ## Mask future climate rasters to extent of species occurrences
  future26 <- mask(future26, sppBuffer)
  future85 <- mask(future85, sppBuffer)
  
  ## save the rasters
  writeRaster(predictOut, paste0("Outputs//Current_",spp,".tif", overwrite=T))
  writeRaster(predict26, paste0("Outputs//Future26_",spp,".tif", overwrite=T))
  writeRaster(predict85, paste0("Outputs//Future85_",spp,".tif", overwrite=T))

  
  #species area
  #range of species estimate
  calcsppRange <- function(x){
    spprange<- x
    spprange[spprange <= 0] <- NA
    spprange<-spprange[!is.na(spprange)]
    spprange<- ncell(spprange)/ncell(current)
    return(spprange)
  }
    
  ## Best MaxEnt Model
  bestMax <- SppMax@results[which (SppMax@results$delta.AICc==0),]
  bestModel <- SppMax@models[[which (SppMax@results$delta.AICc==0)]]
  varImp <- var.importance(bestModel)
  
  ## Output data
  outData <- data.frame(species= spp, ## label the species
                        auc=bestMax[,"train.AUC"], ## model statistic
                        features=bestMax[,"features"], regularization=bestMax[,"rm"], ## maxent parameters
                        npresence=nrow(tempSpp), nabsence=background,  ## number of presence or absences
                        VIFexclude=paste0(varsExcluded, collapse="", sep=";"), ## colinear variables removed
                        importantVars=paste0(varImp$importance, collapse=";"), importantValue=paste0(varImp$permutation.importance, collapse=";"), percentContribution=paste0(varImp$percent.contribution, collapse=";"), 
                        currentRange= calcsppRange(predictOut), ## Range characteristics
                        future26Range= calcsppRange(predict26), ## Range characteristics
                        future85Range= calcsppRange(predict85)) ## Range characteristics
                        
  print("Predict distribution complete - raster written - datefile saved")
  
  write.csv(outData, paste0("Outputs//",spp,"SDMoutput.csv"), row.names=FALSE)
}

for(i in 1:44){
  BBsdm( spp=speciesList[i], spDF = bombus,  biasfile=bias)
  print(i)
}

library(tidyverse)
MaxEntcsv <- list.files(path ="Outputs/speciesoutput", pattern = ".csv", full.names = TRUE)
MaxEntlist<-MaxEntcsv %>% 
  setNames(nm= .) %>%
  map_dfr(~read_csv(.x, col_types=cols(), col_names=FALSE), .id="file_name")
write.csv(MaxEntlist, "Outputs/speciesoutput/AllMaxEnts.csv")

########################################################################
#Prioritizr analyses
########################################################################

#libraries for prioritizr
library(prioritizr)
library(raster)
library(rgdal)
library(gurobi)
library(dplyr)

#datafile prep

polyCAN <- getData("GADM", country = "CAN", level = 1)


#need to create rasterstack of the masked outputs, crop to Canada, and remove low predicted values
canada1 <- raster("CanadaPolyRaster5km.grd")

currentlist<-list.files("Outputs//current",full.names=T)
currentstack<-stack(currentlist)
currentstackmin<-currentstack
currentstackmin <- crop(currentstackmin, canada1)
currentstackmin[canada1==1 & is.na(currentstackmin)] <- 0 ## assign non-predicted areas of Canada zero values
currentstackmin<-mask(currentstackmin, polyCAN)  ## clip predicted area to Canada boundaries
currentstackmin[currentstackmin<0.1]<-0
writeRaster(currentstackmin, "Outputs//CurrentStack.grd", overwrite=T)

future26list<-list.files("Outputs//future26",full.names=T)
future26stack<-stack(future26list)
future26stackmin<-future26stack
future26stackmin<-crop(future26stackmin, canada1)
future26stackmin[canada1==1 & is.na(future26stackmin)] <- 0
future26stackmin<-mask(future26stackmin, polyCAN) 
future26stackmin[future26stackmin<0.1]<-0
writeRaster(future26stackmin, "Outputs//Future26Stack.grd", overwrite=T)

future85list<-list.files("Outputs//future85",full.names=T)
future85stack<-stack(future85list)
future85stackmin<-future85stack
future85stackmin<-crop(future85stackmin, canada1)
future85stackmin[canada1==1 & is.na(future85stackmin)] <- 0
future85stackmin<-mask(future85stackmin, polyCAN) 
future85stackmin[future85stackmin<0.1]<-0
writeRaster(future85stackmin, "Outputs//Future85Stack.grd", overwrite=T)

####
###datafiles made from prep above
###
currentstack<-stack("Outputs//CurrentStack.gri")
future26stack<-stack("Outputs//Future26Stack.grd")
crs(future26stack)<-"+proj=longlat +datum=WGS84 +no_defs" 
future85stack<-stack("Outputs//Future85Stack.grd")
crs(future85stack)<-"+proj=longlat +datum=WGS84 +no_defs" 
canada <- raster("CanadaPolyRaster5km.grd") # cost file
canada[canada > 0] <- 1 # set cost of each pixel to be 1
terrprotectedareas <- raster("CanadaTerrestrialProtectedAreas.tif")
newprotect <- terrprotectedareas
newprotect[is.na(newprotect)] <- 0
newprotect[newprotect != 0 ] <- 1 # assigning areas that are not protected areas to 0 and areas that are to 1
newprotect2 <- resample(newprotect, currentstack, method="ngb")
newprotect2 <- crop(newprotect2, canada)
landcover <- stack("newlandstack.tif") # MODIS landcover raster


## add min set objective function
library(Rsymphony)

bb_priormin <- function(pu = canada, features = currentstack, rel_tar = 0.17, 
                        raster_name = "", np = newprotect2,
                        lc = landcover) {
  
  p0 <- problem(pu, features) %>%
    add_min_set_objective() %>%
    add_relative_targets(rel_tar) %>%
    add_binary_decisions() %>%
    add_rsymphony_solver(0.15)
  s0 <- solve(p0)
  if(raster_name != ""){
    ## Purge raster attributes and convert format
    rasterToWrite <- terra::rast(s0)
    terra::writeRaster(rasterToWrite, raster_name)
  }
  
  f0 <- eval_feature_representation_summary(p0, s0) # checking species representation
  
  # determining the amount of solution is within terrestrial protected areas
  zn_p <- zonal(s0, np, fun = "sum") # summing the amount of pixels within and outside of protected areas
  #     zone   sum
  # [1,]    0 86537
  # [2,]    1 16456   ###16,456 pixels from s2b are within protected areas
  fr_s <- freq(s0) ## 102993 pixels are in the prioritizr solution for s2b (values = 1)
  ## 15.98% of the solution is in protected areas
  
  # determining which landcover classes are within the solution
  s0landcover <- data.frame()
  for (i in 1:nlayers(lc)) {
    tempstats <- zonal(s0, lc[[i]], fun = "sum") ## create temporary dataset of zonal statistics
    temppercent <- tempstats[2, ] / freq(s0)[2, ]
    tempdata <- data.frame(lc = i, percent = temppercent[2] * 100)
    s0landcover <- rbind(s0landcover, tempdata)
    rownames(s0landcover) <- 1:nrow(s0landcover)
  }
  ## convert data classes to dataframes
  speciesData <- data.frame(f0)
  zn_p <- data.frame(zn_p)
  fr_s <- data.frame(fr_s)
  ## calculate percent protected area, write csv for output files, produce list
  speciesData[,"percentProtected"] <- zn_p[zn_p$zone==1, "sum"]/fr_s[2, "count"]*100
  write.csv(speciesData, paste0(gsub(".tif","", raster_name),"species.csv"), row.names=F)
  write.csv(s0landcover, paste0(gsub(".tif","", raster_name),"landcover.csv"), row.names=F)
  return(list(p0 = p0, s0 = s0, f0 = f0, zn_p = zn_p, fr_s = fr_s, s0landcover = s0landcover))
}


##########
##Min set function with higher targets for at-risk 
##########

#risk1current<-c(0.21, 0.17, 0.17, 0.17, 0.17, 0.21, 0.17, 0.21, 0.17, 0.17, 0.21, 0.21, 0.21, 0.17,
                #0.17, 0.17, 0.17, 0.17, 0.17, 0.21, 0.21, 0.17, 0.17, 0.21, 0.21, 0.17, 0.21, 0.21,
                #0.17, 0.21, 0.17, 0.17, 0.17, 0.17, 0.17, 0.21, 0.17, 0.21, 0.21)
#risk1future<-c(0.21, 0.17, 0.17, 0.17, 0.17, 0.21, 0.17, 0.21, 0.17, 0.17, 0.21, 0.21, 0.21, 0.21, 
               #0.17, 0.21, 0.17, 0.17, 0.17, 0.17, 0.17, 0.21, 0.21, 0.17, 0.17, 0.21, 0.21, 0.17,
               #0.21, 0.21, 0.17, 0.21, 0.17, 0.17, 0.17, 0.17, 0.17, 0.21, 0.17, 0.17, 0.17)
#risk2current<-c(0.37, 0.3, 0.3, 0.3, 0.3, 0.37, 0.3, 0.37, 0.3, 0.3, 0.37, 0.37, 0.37, 0.3,
                #0.3, 0.3, 0.3, 0.3, 0.3, 0.37, 0.37, 0.3, 0.3, 0.37, 0.37, 0.3, 0.37, 0.37,
                #0.3, 0.37, 0.3, 0.3, 0.3, 0.3, 0.3, 0.37, 0.3, 0.37, 0.37)
#risk2future<-c(0.37, 0.30, 0.30, 0.30, 0.30, 0.37, 0.30, 0.37, 0.30, 0.30, 0.37, 0.37, 0.37, 0.37, 
               #0.30, 0.37, 0.30, 0.30, 0.30, 0.30, 0.30, 0.37, 0.37, 0.30, 0.30, 0.37, 0.37, 0.30,
               #0.37, 0.37, 0.30, 0.37, 0.30, 0.30, 0.30, 0.30, 0.30, 0.37, 0.30, 0.30, 0.30)
#risk3current<-c(0.62, 0.5, 0.5, 0.5, 0.5, 0.62, 0.5, 0.62, 0.5, 0.5, 0.62, 0.62, 0.62, 0.5,
                #0.5, 0.5, 0.5, 0.5, 0.5, 0.62, 0.62, 0.5, 0.5, 0.62, 0.62, 0.5, 0.62, 0.62,
                #0.5, 0.62, 0.5, 0.5, 0.5, 0.5, 0.5, 0.62, 0.5, 0.62, 0.62)
#risk3future<-c(0.62, 0.50, 0.50, 0.50, 0.50, 0.62, 0.50, 0.62, 0.50, 0.50, 0.62, 0.62, 0.62, 0.62, 
               #0.50, 0.62, 0.50, 0.50, 0.50, 0.50, 0.50, 0.62, 0.62, 0.50, 0.50, 0.62, 0.62, 0.50,
               #0.62, 0.62, 0.50, 0.62, 0.50, 0.50, 0.50, 0.50, 0.50, 0.62, 0.50, 0.50, 0.50)

######
####Running min set prioritizr analyses
######

### Problem 1a add_min_set objective, 17%, current climate
p1a <- bb_priormin(raster_name =  "prioritizrResults//s1aMinSet17.tif") 

### Problem 2a add_min_set objective, 30% (Canada 2050 biodiversity target), current climate
p2a <- bb_priormin( rel_tar = 0.3, raster_name =  "prioritizrResults//s2aMinSet30.tif")

### Problem 3a add_min_set objective, 50% (nature needs half), current climate
p3a<- bb_priormin( rel_tar = 0.5, raster_name =  "prioritizrResults//s3aMinSet50.tif")

### Problem 1b add_min_set_objective, 17% target, future climate rcp 2.6
p1b <- bb_priormin(features = future26stack, raster_name =  "prioritizrResults//s1bMinSet17RCP26.tif")

### Problem 2b add_min_set_objective, 30% target, future climate rcp 2.6
p2b<- bb_priormin(features = future26stack, rel_tar = 0.3, raster_name =  "prioritizrResults//s2bMinSet30RCP26.tif")

### Problem 3b add_min_set objective, 50% (nature needs half), future climate rcp 2.6
p3b <- bb_priormin(features= future26stack, rel_tar = 0.5, raster_name =  "prioritizrResults//s3bMinSet50RCP26.tif")

### Problem 1c add_min_set_objective, 17% target, future climate rcp 8.5
p1c <- bb_priormin(features = future85stack, raster_name =  "prioritizrResults//s1cMinSet17RCP85.tif")

### Problem 2c add_min_set_objective, 30% target, future climate rcp 8.5
p2c<- bb_priormin(features = future85stack, rel_tar = 0.3, raster_name =  "prioritizrResults//s2cMinSet30RCP85.tif")

### Problem 3c add_min_set objective, 50% (nature needs half), future climate rcp 8.5
p3c <- bb_priormin(features= future85stack, rel_tar = 0.5, raster_name =  "prioritizrResults//s3cMinSet50RCP85.tif")

####

###function adding eval_replacement_importance to add min set
#instead of re-running above, writing a new function to test for evaluating solution importance
#changing  to proportional decisions to make it run faster with negligible difference

bb_replace <- function(pu = canada, features = currentstack2, rel_tar = 0.17) {
  
  p0 <- problem(pu, features) %>%
    add_min_set_objective() %>%
    add_relative_targets(rel_tar) %>%
    add_proportion_decisions() %>%
    add_gurobi_solver(threads = 6)
  s0 <- solve(p0)
  
r0<-eval_replacement_importance(p0, s0, threads=3)

writeRaster(r0, paste0("prioritizrResults/replacement/", expression(features),rel_tar,".tif"))
}

###
##re-running the prioritizr analysis with replacement importance

#current climate, 17% target
r1 <- bb_replace()
#current climate, 30% target
r2 <- bb_replace (rel_tar = 0.3)
#current climate, 50% target
r3<-bb_replace(rel_tar=0.5)

################################
###Prioritizr results figures
################################

library(raster)
library(ggplot2)
library(virdis)
library(rgdal)

prov<-readOGR(dsn=".", layer="provinces")
prov2<-spTransform(prov, CRS("+proj=longlat +datum=WGS84 +no_defs"))

#### prioritizr results <0.5 removed
s1a<-raster("prioritizrResults/s1aMinSet17.tif")
s2a<-raster("prioritizrResults/s2aMinSet30.tif")
s3a<-raster("prioritizrResults/s3aMinSet50.tif")
s1b<-raster("prioritizrResults/s1bMinSet17RCP26.tif")
s2b<-raster("prioritizrResults/s2bMinSet30RCP26.tif")
s3b<-raster("prioritizrResults/s3bMinset50RCP26.tif")
s1c<-raster("prioritizrResults/s1cMinSet17RCP85.tif")
s2c<-raster("prioritizrResults/s2cMinSet30RCP85.tif")
s3c<-raster("prioritizrResults/s3cMinSet50RCP85.tif")
ecoreg <- readOGR("NA_CEC_Eco_Level1.shp")
library(sf)

canreg <-readOGR("ecozones.shp")
canecoreg <- mask (ecoreg, s1a)

canecoplot <- ggplot(st_as_sf(canreg)) +
  geom_sf(aes(fill = ZONE_NAME)) +
  scale_fill_manual(
    values = c("#9E0142", "#F93943", "#E4572E", "#FDAE61", "#F5CB5C",
               "#DDE5B6", "#A7C957",  "#386641", "#66A182", "#66C2A5", 
               "#A9D6E5", "#3288BD", "#445E93","#6A4C93", "#6930C3",
               "#2E4057") , 
    name = "Ecoregions Level 1"
  )
canecoplot

pdf("LeftMinSet.pdf", width = 10, height = 8, useDingbats = F)
plot.new()
par(mfrow = c(2, 2))
par(mar=c(0, 0, 0, 0))
plot(s1a, col=c("#CFC9CF", "#232ED1"),  legend=F, xaxt = "n", yaxt = "n", box(col="white"))
plot(prov2, add=T)
plot(s1b, col=c("#D8CC3400", adjustcolor ("#E90023", alpha.f=0.7)), 
                 legend=F,  xaxt = "n", yaxt = "n", add=T, box(col="white"))
plot(s1c, col=c("#E3E3E300", adjustcolor ("#FFEE00", alpha.f=0.65)), 
      legend=F,  xaxt = "n", yaxt = "n", add=T, box(col="white"))
# dev.off()
# 
# pdf("MidMinSet.pdf", width = 16, height = 16, useDingbats = F)
plot(s2a, col=c("#CFC9CF", "#232ED1"), legend=F, xaxt = "n", yaxt = "n", box(col="white"))
plot(prov2, add=T)
plot(s2b, col=c("#D8CC3400", adjustcolor ("#E90023", alpha.f=0.7)), 
      legend=F,  xaxt = "n", yaxt = "n", add=T, box(col="white"))
plot(s2c, col=c("#E3E3E300", adjustcolor ("#FFEE00", alpha.f=0.65)), 
      legend=F,  xaxt = "n", yaxt = "n", add=T, box(col="white"))

# dev.off()
# 
# pdf("RightMinSet.pdf", width = 16, height = 16, useDingbats = F)
plot(s3a, col=c("#CFC9CF", "#232ED1"), legend=F, xaxt = "n", yaxt = "n", box(col="white"))
plot(prov2, add=T)
plot(s3b, col=c("#D8CC3400", adjustcolor ("#E90023", alpha.f=0.7)), 
      legend=F,  xaxt = "n", yaxt = "n", add=T, box(col="white"))
plot(s3c, col=c("#E3E3E300", adjustcolor ("#FFEE00", alpha.f=0.65)), 
      legend=F,  xaxt = "n", yaxt = "n", add=T, box(col="white"))
dev.off()


####landcover bar plot figure
library(ggplot2)
library(tidyverse)
library(purrr)

#first need to combine output files
landcovercsv <- list.files(path = "prioritizrResults/landcover", pattern = ".csv", full.names = TRUE)
landlist<-landcovercsv %>% 
  setNames(nm= .) %>%
  map_dfr(~read_csv(.x, col_types=cols(), col_names=FALSE), .id="file_name")
write.csv(landlist, "prioritizrResults/landcover/AllLandcover.csv")


## get proportionate land cover areas
propLandcover <- sapply(1:nlayers(landcover), function(i){
  tempRaster <- landcover[[i]]
  tempRaster[tempRaster > 0]  <- 1
  freq(tempRaster,  value = 1, useNA = "no")
})
landdat.new <- landdat %>% filter(X1 %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14",
                           "15", "16", "17"))
landdat.new2 <-  landdat.new %>%  rename (landcoverClass = X1, 
                     percentcover = X2)

landdat2 <- landdat.new2 %>% 
  mutate(modelCode = substr(file_name, 29, 31))

landdat2[,"propLand"] <- rep(unique(propLandcover), times = 9)

landdat3 <- landdat2 %>% 
  mutate(weighted = as.numeric(percentcover)/propLand)

propsum<- keeps <-landdat3 %>% group_by(landcoverClass) %>% summarize(n=min(weighted))
write.csv(propsum, "landpropweighted.csv")
keeps <-landdat3 %>% group_by(landcoverClass) %>% summarize(n=min(weighted)) %>% filter(n>0.00003)
landdat4 <- landdat3 %>% filter(landcoverClass %in% keeps$landcoverClass)

current <- c("s1a", "s2a", "s3a")
future2.6 <- c("s1b", "s2b", "s3b")
future8.5 <- c("s1c", "s2c", "s3c")

landdat5 <- landdat4 %>% 
  mutate(climate = ifelse(modelCode %in% current, "current", 
                          ifelse(modelCode %in% future2.6, "future2.6", "future8.5")))

#calculating standard error for barplot
error <-landdat5 %>% group_by(landcoverClass, climate) %>% summarise(mean = mean(weighted), 
                                                                     se=sd(weighted)/sqrt(3)) 

#making plot
ggplot(data=error, aes(x=reorder(landcoverClass, -mean), y=mean, fill=climate))+
  geom_bar(stat="summary",  color="black", position=position_dodge())  +
  geom_errorbar(data=error, aes(x=reorder(landcoverClass, -mean), ymin=mean-se, ymax=mean+se), 
                width=0, position=position_dodge(width=0.9)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_fill_manual(values=c("#1D29AD", "#E00031", "#FFD500")) +
  scale_x_discrete(breaks = c("12", "10", "7", "5", "15", "8"),
                   labels = c("Croplands", "Grasslands", "Open Shrublands", "Mixed Forest", "Snow or Ice",
                              "Woody Savannah")) + 
  xlab("Landcover Class")+ ylab("Weighted Percent Cover")+scale_y_continuous(expand=c(0,0))

###Proportion of common and at-risk species within protected areas

##combinbing csv files for species representation outputs
library(tidyverse)
MinSetcsv <- list.files(path ="prioritizrResults/speciesOutput", pattern = ".csv", full.names = TRUE)
MinSetlist<-MinSetcsv %>% 
  setNames(nm= .) %>%
  map_dfr(~read_csv(.x, col_types=cols(), col_names=FALSE), .id="file_name")
write.csv(MinSetlist, "prioritizrResults/speciesOutput/AllMinSet.csv")

##proportion of at-risk vs. common species within the solution
#summarized excel sheet in excel, added common/at-risk status
object<-read.csv("prioritizrResults/speciesOutput/AllMinSet.csv")

object.new <-  object %>%  rename (species = X2, 
                                   totalAmount = X3,
                                   absoluteHeld = X4,
                                   relativeHeld = X5,
                                   percentProtect = X6)

object.new2 <- object.new [!is.na(as.numeric(object.new$relativeHeld)),]

object2 <- object.new2 %>% 
  mutate(modelCode = substr(file_name, 33, 35))

current <- c("s1a", "s2a", "s3a")
future2.6 <- c("s1b", "s2b", "s3b")
future8.5 <- c("s1c", "s2c", "s3c")

target17 <-c("s1a", "s1b", "s1c")
target30 <-c("s2a", "s2b", "s2c")
target50 <-c("s3a", "s3b", "s3c")

object3 <- object2 %>% 
  mutate(climate = ifelse(modelCode %in% current, "current", 
                          ifelse(modelCode %in% future2.6, "future2.6", "future8.5")))

object3 <- object3 %>% 
  mutate(target = ifelse(modelCode %in% target17, "17", 
                         ifelse(modelCode %in% target30, "30", "50")))
write.csv(object3, "prioritizrResults/speciesOutput/AllMinSetCleaner.csv")

clean <- read.csv("prioritizrResults/speciesOutput/AllMinSetCleaner.csv")

atRisk <-c("affinis", "caliginosus", "cryptarum", "fervidus", "flavidus", "fraternus", "jonellus",
           "kirbiellus", "morrisoni", "neoboreus", "occidentalis", "pensylvanicus", "polaris", 
           "suckleyi", "terricola")

clean <- clean %>% 
  mutate(status = ifelse(species %in% atRisk, "atRisk", "stable"))


error <- clean %>% group_by(status, climate, target) %>% 
  summarise(mean = mean(relativeHeld, na.rm=TRUE), 
            se=sd(relativeHeld, na.rm=TRUE)/sqrt(length(relativeHeld[!is.na(relativeHeld)])))

ggplot(data=error, aes(x=climate, y=mean, fill=status))+
  geom_bar(stat="summary", color="black", position=position_dodge()) + facet_grid(cols=vars(target)) +
  geom_errorbar(data=error, aes (x=climate, ymin=mean-se, ymax=mean+se), width=0, position=position_dodge(width=0.9)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_fill_manual(values=c("#EAC435", "#345995"))+
  xlab("")+ ylab("Proportion of distribution")+scale_y_continuous(expand=c(0,0), breaks=seq(0, .70, .10))

#### percent area with >50% suitability figure
library(ggplot2)
library(ggbreak)
library(gridExtra)
library(forcats)
library(raster)

#the species range did not calculate properly due to a typo in the function so recalculating here
# IT DID IN FEB 2023 SO IGNORE NOW 
canada1 <- raster("CanadaPolyRaster5km.grd")

currentlist<-list.files("Outputs//current",full.names=T)
currentstack<-stack(currentlist)
currentstackmin<-currentstack
currentstackmin <- crop(currentstackmin, canada1)
currentstackmin[canada1==1 & is.na(currentstackmin)] <- 0 ## assign non-predicted areas of Canada zero values
currentstackmin<-mask(currentstackmin, polyCAN) 
currentstack2<-currentstackmin

future26list<-list.files("Outputs//future26",full.names=T)
future26stack<-stack(future26list)
future26stackmin<-future26stack
future26stackmin<-crop(future26stackmin, canada1)
future26stackmin[canada1==1 & is.na(future26stackmin)] <- 0
future26stackmin<-mask(future26stackmin, polyCAN) 
future26stack2<-future26stackmin

future85list<-list.files("Outputs//future85",full.names=T)
future85stack<-stack(future85list)
future85stackmin<-future85stack
future85stackmin<-crop(future85stackmin, canada1)
future85stackmin[canada1==1 & is.na(future85stackmin)] <- 0
future85stackmin<-mask(future85stackmin, polyCAN) 
future85stack2<-future85stackmin

bumblebeenames <- c("affinis", "appositus", "auricomus", "bifarius", "bimaculatus", "bohemicus", 
                    "borealis", "caliginosus", "centralis", "citrinus", "crotchii", "cryptarum", 
                    "fervidus", "flavidus", "flavifrons", "fraternus", "frigidus", "griseocollis",
                    "huntii", "impatiens", "insularis", "jonellus", "kirbiellus", 
                    "melanopygus", "mixtus", "morrisoni", "neoboreus", "nevadensis", 
                    "occidentalis", "pensylvanicus", "perplexus", "polaris", "rufocinctus", 
                    "sandersoni", "sitkensis", "suckleyi", "sylvicola", "ternarius", "terricola",
                    "vagans", "vandykei", "vosnesenskii")
names(currentrangestack) <- bumblebeenames

calcsppRange <- function(x){
  spprange<- x
  spprange[spprange <= 0.01] <- NA
  spprange<-spprange[!is.na(spprange)]
  spprange<- ncell(spprange)/ncell(x)
  spprange<-(spprange)*1717662
  return(spprange)
}
outRange <- lapply(1:nlayers(currentstack2), function(i){
  calcsppRange(currentstack2[[i]])
})
currentRanges <- data.frame(Species = bumblebeenames, CurrentRange= do.call(c, outRange))

outRange <- lapply(1:nlayers(future26stack2), function(i){
  calcsppRange(future26stack2[[i]])
})
future26Ranges <- data.frame(Species = bumblebeenames, future26Range= do.call(c, outRange))

outRange <- lapply(1:nlayers(future85stack2), function(i){
  calcsppRange(future85stack2[[i]])
})
future85Ranges <- data.frame(Species = bumblebeenames, future85Range= do.call(c, outRange))


prop50 <- currentRanges %>% left_join(future26Ranges) %>% left_join(future85Ranges)
prop50[,c("CurrentRange","future26Range","future85Range")] <- prop50[,c("CurrentRange","future26Range","future85Range")]  +1


Sppcsv <- list.files(path ="Outputs/speciesoutput", pattern = ".csv", full.names = TRUE)
Sppcsv<-Sppcsv %>% 
  setNames(nm= .) %>%
  map_dfr(~read_csv(.x, col_types=cols(), col_names=FALSE), .id="file_name")
write.csv(Sppcsv, "Outputs/speciesoutput/AllMaxents.csv")

alldem <- read.csv("Outputs/speciesoutput/AllMaxents.csv")

alldem.new <-  alldem %>%  rename (species = X1, 
                                   AUC = X2,
                                   CBI = X3,
                                   features = X4,
                                   regularization = X5,
                                   npresence = X6,
                                   nabsence = X7,
                                   VIFexclude = X8,
                                   importantVar = X9,
                                   importantValue = X10,
                                   percentContribution = X11,
                                   currentRange = X12,
                                   future2.6Range = X13,
                                   future8.5Range = X14)

alldem.new2 <- alldem.new [!is.na(as.numeric(alldem.new$currentRange)),]
write.csv(alldem.new2, "Outputs/speciesoutput/AllMaxentsCleaner.csv")

allrange <-subset(alldem.new2, select = c("species", "currentRange", "future2.6Range", "future8.5Range"))
allrange[,2:4]<- apply(allrange[,2:4], 2, as.numeric)

library(tidyverse)
prop <- allrange %>%  
  dplyr::select(species = species, current= currentRange,future2.6= future2.6Range, future8.5=future8.5Range)%>% 
  mutate(diff26 = (future2.6-current)/current, diff85 = (future8.5-current)/current) %>% 
  mutate(species = factor(species, levels=species[order(diff26)] )) %>% ## re-order by change in area
  mutate(group = ifelse(diff26 <= -0.001, "decreasing","increasing")) %>% 
  gather(climate, proportion, diff26:diff85) 


plot1 <- ggplot(prop %>% filter(group=="decreasing"), 
                aes(x= fct_rev(species), y=proportion*100, fill=fct_rev(climate))) + 
  geom_bar(stat="identity", position="dodge") + coord_flip() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none",
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_fill_manual(values=c("#FFD500", "#E00031", "#FFD500")) +
  theme(axis.text.y = element_text(face="italic", size=12))+
  theme(axis.text.x = element_text(size = 12)) +
  xlab("")+ ylab("Percent Change in Range from Current")+ ylim(-30,10)
  # scale_y_continuous(expand=c(0,0)) 

plot2 <- ggplot(prop %>% filter(group=="increasing"), 
                aes(x= fct_rev(species), y=proportion*100, fill=fct_rev(climate))) + 
  geom_bar(stat="identity", position="dodge") + coord_flip() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_fill_manual(values=c("#FFD500", "#E00031", "#FFD500"), 
                    labels = c("future RCP 8.5", "future RCP 2.6")) +
  theme(axis.text.y = element_text(face="italic", size=12))+
  theme(axis.text.x = element_text(size = 12)) +
  theme(legend.title=element_blank())+
  theme(legend.position = c(0.5, 0.95)) +
  xlab("")+ ylab("Percent Change in Range from Current")+ ylim(-1,30) 
  # scale_y_continuous(expand=c(0,0))

gridExtra::grid.arrange(plot1, plot2, ncol=2)

### change in latitude and elevation table
#libraries
library(raster)
library(dplyr)

#datafiles
elevation<-raster("CAN_msk_alt.grd")
s1a<-raster("prioritizrResults/s1aMinSet17.tif")
s2a<-raster("prioritizrResults/s2aMinSet30.tif")
s3a<-raster("prioritizrResults/s3aMinSet50.tif")
s1b<-raster("prioritizrResults/s1bMinSet17RCP26.tif")
s2b<-raster("prioritizrResults/s2bMinSet30RCP26.tif")
s3b<-raster("prioritizrResults/s3bMinset50RCP26.tif")
s1c<-raster("prioritizrResults/s1cMinSet17RCP85.tif")
s2c<-raster("prioritizrResults/s2cMinSet30RCP85.tif")
s3c<-raster("prioritizrResults/s3cMinSet50RCP85.tif")

east<-extent(-90, -54.27083, 42.85417, 83.0625) #separating into east and west to account for the Rockies
west<-extent(-141, -90, 41.66667, 83.125)

latitudefun<-function(future, current){
  currentp<-rasterToPoints(current, function(x) x>0, spatial=TRUE)
  futurep<-rasterToPoints(future, function(x) x>0, spatial=TRUE)
  ## Current
  currentE<-crop(current, east, keepres=F, snap="in")
  currentEp<-rasterToPoints(currentE, function(x) x>0, spatial=TRUE)
  currentW<-crop(current, west, keepres=T, snap="in")
  currentWp<-rasterToPoints(currentW, function(x) x>0, spatial=TRUE)
  ## Future
  futureE<-crop(future, east, keepres=T, snap="in")
  futureW<-crop(future, west, keepres=T, snap="in")
  differenceEp<-rasterToPoints(futureE, function(x) x>0, spatial=TRUE)
  differenceWp<-rasterToPoints(futureW, function(x) x>0, spatial=TRUE)
  dataOut <- data.frame( extent=c("current Canada","future Canada","current east", "future east","current west","future west"),
                         latitude=c(mean(coordinates(currentp)[,2]),mean(coordinates(futurep)[,2]),
                                    mean(coordinates(currentEp)[,2]),mean(coordinates(differenceEp)[,2]),
                                    mean(coordinates(currentWp)[,2]),mean(coordinates(differenceWp)[,2])))
  return(dataOut)
}

elevationfun<-function(future, current){
  currentp<-rasterToPoints(current, function(x) x>0, spatial=TRUE)
  futurep<-rasterToPoints(future, function(x) x>0, spatial=TRUE)
  ## Current
  currentE<-crop(current, east, keepres=F, snap="in")
  currentEp<-rasterToPoints(currentE, function(x) x>0, spatial=TRUE)
  currentW<-crop(current, west, keepres=T, snap="in")
  currentWp<-rasterToPoints(currentW, function(x) x>0, spatial=TRUE)
  ## Future
  futureE<-crop(future, east, keepres=T, snap="in")
  futureW<-crop(future, west, keepres=T, snap="in")
  futureEp<-rasterToPoints(futureE, function(x) x>0, spatial=TRUE)
  futureWp<-rasterToPoints(futureW, function(x) x>0, spatial=TRUE)
  
  elevfut<-raster::extract(elevation, futurep)
  elevcur<-raster::extract(elevation, currentp)
  elevcurE<-raster::extract(elevation, currentEp)
  elevcurW<-raster::extract(elevation, currentWp)
  elevfutE<-raster::extract(elevation, futureEp)
  elevfutW<-raster::extract(elevation, futureWp)
  dataOut <- data.frame( extent=c("future Canada","current Canada","current east", "current west","future east","future west"),
                         elevate=c(mean(elevfut, na.rm=T),mean(elevcur, na.rm=T),
                                   mean(elevcurE, na.rm=T),mean(elevcurW,na.rm=T),
                                   mean(elevfutE,na.rm=T),mean(elevfutW,na.rm=T)))
  return(dataOut)
}

currentList <- list(s1a, s1a, s2a, s2a, s3a, s3a)
futureList<-list(s1b, s1c, s2b, s2c, s3b, s3c)

#elevation
allChange <- data.frame()
for(i in 1:6){
  temp <- elevationfun(futureList[[i]], currentList[[i]])
  temp[,"iteration"] <- names(futureList[[i]])
  allChange <- rbind(allChange,temp)
}
allChange[,"scenario"]<- gsub(".*RCP","",allChange$iteration)


allDone<-allChange %>% group_by (extent, scenario) %>% summarise(mean=mean(elevate), se=sd(elevate)/sqrt(length(elevate))) %>% data.frame()
allDone[1:6,"scenario"] <- "current"
allDone <- allDone[!duplicated(allDone),]

#latitude
allChange <- data.frame()
for(i in 1:6){
  temp <- latitudefun(futureList[[i]], currentList[[i]])
  temp[,"iteration"] <- names(futureList[[i]])
  allChange <- rbind(allChange,temp)
}
allChange[,"scenario"]<- gsub(".*RCP","",allChange$iteration)

allDone2<-allChange %>% group_by (extent, scenario) %>% summarise(mean=mean(latitude), se=sd(latitude)/sqrt(length(latitude))) %>% data.frame()
allDone2[1:6,"scenario"] <- "current"
allDone2 <- allDone2[!duplicated(allDone2),]



###########################
###Climate Data Prepping
###########################

##Current climate
currentclim<-list.files("current", full.names = TRUE, pattern=".bil")
currentclim.all<-stack(currentclim)
currentclim.all <- crop(currentclim.all, poly) 
currentclim.all <- mask(currentclim.all, poly)
writeRaster(currentclim.all, "Current.tif")

###Prepping future climate data from worldclim

##RCP 2.6
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

##RCP 8.5
bioclim <- list.files("RCP8.5", full.names = TRUE, pattern=".tif")
bioclim.all <- stack(bioclim)
bioclim.all <- crop(bioclim.all, poly) 
bioclim.all <- mask(bioclim.all, poly) 

## add letter to the end of the bioclim vars to separate out the 1s from the 10s
rasterNames <- paste0(names(bioclim.all),"e")

## Iterate through all 19 bioclim variables
future85 <- stack()
future85 <- lapply(1:19, function(i){
  biovar <- paste0("bi50",i,"e") ## select bioclimate variable per iteration
  bioclimTemp <- bioclim.all[[grep(biovar,rasterNames)]] ## select that climate variable
  meanTemp <- mean(bioclimTemp) ## average across GCMs
  future85 <- stack(future85, meanTemp)
})
## convert list of rasters to raster stack
all85 <- do.call(stack, future85)
names(all85) <- paste0(rep("bio",19),1:19)
writeRaster(all85, "Future8.5.grd")

##################
#Prioritizr data prepping
##################

polyCAN <- getData("GADM", country = "CAN", level = 1)

#Current Climate
currents<-stack("Outputs/CurrentStackAll.tif")
cropstack<-crop(currents, polyCAN)
maskstack<-mask(cropstack, polyCAN) 
stackspmin<-maskstack
stackspmin <- reclassify(stackspmin, cbind(-Inf, 0.2, 0), right=FALSE)
writeRaster(stackspmin, "CurrentStack2.5min.grd", overwrite=T)

#Future climate RCP 2.6
futures26<-stack("Outputs/Future26StackAll.tif")
cropstack<-crop(futures26, polyCAN)
maskstack<-mask(cropstack, polyCAN) 
stackspmin<-maskstack
stackspmin <- reclassify(stackspmin, cbind(-Inf, 0.2, 0), right=FALSE)
writeRaster(stackspmin, "Future26stack2.5min.grd", overwrite=T)

#Future climate RCP 8.5
futures85<-stack("Outputs/Future85StackAll.tif")
cropstack<-crop(futures85, polyCAN)
maskstack<-mask(cropstack, polyCAN) 
stackspmin<-maskstack
stackspmin <- reclassify(stackspmin, cbind(-Inf, 0.2, 0), right=FALSE)
writeRaster(stackspmin, "Future85stack2.5min.grd", overwrite=T)

## what land cover classes are within protected areas?

prozone <- zonal(landcover, newprotect2, fun= "sum")
prozone
#ncells of landcover = 2110395
prozone2 <- prozone/2110395



zn_p <- zonal(s0, np, fun = "sum") # summing the amount of pixels within and outside of protected areas
#     zone   sum
# [1,]    0 86537
# [2,]    1 16456   ###16,456 pixels from s2b are within protected areas
fr_s <- freq(s0) #



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

