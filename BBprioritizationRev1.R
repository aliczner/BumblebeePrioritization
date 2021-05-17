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
current<-stack("Currentclim.tif")
future26<-stack("Future2.6.grd")
future85<-stack("Future8.5.grd")

## Rename raster stacks to be the same names
names(current) <- paste0(rep("bio",19),1:19)
names(future26) <- paste0(rep("bio",19),1:19)
names(future85) <- paste0(rep("bio",19),1:19)

########################
### SDM Climate Analyses for current and future
#######################

speciesList <- sort(as.character(unique(bombus$species)))

currents<-stack("Outputs/CurrentStack.tif")

BBsdm<-function(spp, spDF, biasfile)  {
  ## Load current climate
  climate=current
  
  ## Data set up
  tempSpp <- subset(spDF, species == spp)
  background <- nrow(tempSpp)*10
  
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
  predictOut <- predict(climateSp, SppMax@models[[which (SppMax@results$delta.AICc==0)]], type="logistic") ## predictOut raster
  predict26 <- predict(future26, SppMax@models[[which (SppMax@results$delta.AICc==0)]], type="logistic") ## predictOut raster
  predict85 <- predict(future85, SppMax@models[[which (SppMax@results$delta.AICc==0)]], type="logistic") ## predictOut raster
  
  ## save the rasters
  writeRaster(predictOut, paste0("Outputs//Current_",spp,".tif", overwrite=T))
  writeRaster(predict26, paste0("Outputs//Future26_",spp,".tif", overwrite=T))
  writeRaster(predict85, paste0("Outputs//Future85_",spp,".tif", overwrite=T))
  
  #may need to assign CRS at some point
  
  #species area for current climate
  #range of species estimate
  calcsppRange <- function(x){
    spprange<- x
    spprange[spprange>0] <- NA
    spprange<-spprange[!is.na(spprange)]
    spprange<- ncell(spprange)/ncell(current)
    spprange<-(spprange)*1717662
    return(spprange)
  }
    
  #>50%
  calcspp50 <- function(x){
    spp50<- x
    spp50[spp50<0.5] <-  NA 
    spp50<-spp50[!is.na(spp50)]
    spp50<- ncell(spp50)/ncell(current)
    spp50<-(spp50)*1717662  #area of Canada multiplied by the percent pixels occupied by the species
    return(spp50)
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
                        currentRange= calcsppRange(predictOut), currentArea50km = calcspp50(predictOut), ## Range characteristics
                        future26Range= calcsppRange(predict26), future26Area50km = calcspp50(predict26), ## Range characteristics
                        future85Range= calcsppRange(predict85), future85Area50km = calcspp50(predict85)) ## Range characteristics
                        
  print("Predict distribution complete - raster written - datefile saved")
  
  write.csv(outData, paste0("Outputs//",spp,"SDMoutput.csv"), row.names=FALSE)
}

for(i in 1:44){
  BBsdm( spp=speciesList[i], spDF = bombus,  biasfile=bias)
  print(i)
}
#######################################################################
#Prioritizr analyses
#######################################################################

#libraries
library(prioritizr)
library(raster)
library(rgdal)
library(gurobi)

#datafiles
currentstackmin<-stack("CurrentStack2.5min.grd")
currentstackmin[[16]][currentstackmin[[16]]>0] <-0 #very small value predicted in current for fraternus which does not occur in Canada so removing
future26stackmin<-stack("Future26Stack2.5min.grd")
crs(future26stackmin)<-"+proj=longlat +datum=WGS84 +no_defs" 
future85stackmin<-stack("Future85Stack2.5min.grd")
crs(future85stackmin)<-"+proj=longlat +datum=WGS84 +no_defs" 
canada <- raster("CanadaPolyRaster5km.grd") # cost file
canada[canada > 0] <- 1 # set cost of each pixed to be 1
terrprotectedareas <- raster("CanadaTerrestrialProtectedAreas.tif")
newprotect <- terrprotectedareas
newprotect[is.na(newprotect)] <- 0
newprotect[newprotect != 0 ] <- 1 # assigning areas that are not protected areas to 0 and areas that are to 1
newprotect2 <- resample(newprotect, currentstackmin, method="ngb")

landcover <- stack("newlandstack.tif") # MODIS landcover raster

## add min set objective function

bb_priormin <- function(pu = canada, features = currentstackmin, rel_tar = 0.17, 
                        raster_name = "", np = newprotect2,
                        lc = landcover) {
  
  p0 <- problem(pu, features) %>%
    add_min_set_objective() %>%
    add_relative_targets(rel_tar) %>%
    add_binary_decisions() %>%
    add_gurobi_solver(threads = 6)
  s0 <- solve(p0)
  if(raster_name != ""){
    writeRaster(s0, raster_name)
  }
  
  f0 <- eval_feature_representation_summary(p0, s0) # checking species representation
  
  # determining the amount of solution is within terrestrial protected areas
  zn_p <- zonal(s0, np, fun = "sum") # summing the amount of pixes within and outside of protected areas
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


####Min set all species equal

### Problem 1a add_min_set objective, 17%, current climate
p1a <- bb_priormin(raster_name =  "prioritizrResults//s1aMinSet17.tif")

### Problem 2a add_min_set objective, 30% (Canada 2050 biodiversity target), current climate
p2a <- bb_priormin( rel_tar = 0.3, raster_name =  "prioritizrResults//s2aMinSet30.tif")

### Problem 3a add_min_set objective, 50% (nature needs half), current climate
p3a<- bb_priormin( rel_tar = 0.5, raster_name =  "prioritizrResults//s3aMinSet50.tif")

### Problem 1b add_min_set_objective, 17% target, future climate rcp 2.6
p1b <- bb_priormin(features = future26stackmin, raster_name =  "prioritizrResults//s1bMinSet17RCP26.tif")

### Problem 2b add_min_set_objective, 30% target, future climate rcp 2.6
p2b<- bb_priormin(features = future26stackmin, rel_tar = 0.3, raster_name =  "prioritizrResults//s2bMinSet30RCP26.tif")

### Problem 3b add_min_set objective, 50% (nature needs half), future climate rcp 2.6
p3b <- bb_priormin(features= future26stackmin, rel_tar = 0.5, raster_name =  "prioritizrResults//s3bMinSet50RCP26.tif")

### Problem 1c add_min_set_objective, 17% target, future climate rcp 8.5
p1c <- bb_priormin(features = future85stackmin, raster_name =  "prioritizrResults//s1cMinSet17RCP85.tif")

### Problem 2c add_min_set_objective, 30% target, future climate rcp 8.5
p2c<- bb_priormin(features = future85stackmin, rel_tar = 0.3, raster_name =  "prioritizrResults//s2cMinSet30RCP85.tif")

### Problem 3c add_min_set objective, 50% (nature needs half), future climate rcp 8.5
p3c <- bb_priormin(features= future85stackmin, rel_tar = 0.5, raster_name =  "prioritizrResults//s3cMinSet50RCP85.tif")

################
#Prioritizr results figures
###############

## prioritizr results figure 1 add min set objective
s1a<-raster("prioritizrResults/s1aMinSet17.tif")
s2a<-raster("prioritizrResults/s2aMinSet30.tif")
s3a<-raster("prioritizrResults/s3aMinSet50.tif")
s1b<-raster("prioritizrResults/s1bMinSet17RCP26.tif")
s2b<-raster("prioritizrResults/s2bMinSet30RCP26.tif")
s3b<-raster("prioritizrResults/s3bMinset50RCP26.tif")
s1c<-raster("prioritizrResults/s1cMinSet17RCP85.tif")
s2c<-raster("prioritizrResults/s2cMinSet30RCP85.tif")
s3c<-raster("prioritizrResults/s3cMinSet50RCP85.tif")

topleft<-(s1b*1.1)-s1a
topmid<-(s2b*1.1)-s2a
topright<-(s3b*1.1)-s3a
midleft<-(s1c*1.1)-s1a
midmid<-(s2c*1.1)-s2a
midright<-(s3c*1.1)-s3a
botleft<-s1a*s1b*s1c
botmid<-s2a*s2b*s2c
botright<-s3a*s3b*s3c

pdf("FigureMinSet.pdf", width = 16, height = 16, useDingbats = F)
par(mfrow = c(3, 3))
par(mar = c(4.5, 4.5, 0, 0))
plot(topleft, col = c("#D4B483","#E3E3E3", "#416788","#416788"), legend = FALSE, xaxt = "n", cex.axis = 2, las = 1, xlim=c(-145,-52))
plot(topmid, col = c("#D4B483","#E3E3E3", "#416788","#416788"), legend = FALSE, xaxt = "n", yaxt = "n", xlim=c(-145,-52))
plot(topright, col = c("#D4B483","#E3E3E3", "#416788","#416788"), legend = FALSE, xaxt = "n", yaxt = "n", xlim=c(-145,-52))
plot(midleft, col = c("#D4B483","#E3E3E3", "#416788","#416788"), legend = FALSE, xaxt = "n", cex.axis = 2, las = 1, xlim=c(-145,-52))
plot(midmid, col = c("#D4B483","#E3E3E3", "#416788","#416788"), legend = FALSE, xaxt = "n", yaxt = "n", xlim=c(-145,-52))
plot(midright, col = c("#D4B483","#E3E3E3", "#416788","#416788"), legend = FALSE, xaxt = "n", yaxt = "n", xlim=c(-145,-52))
plot(botleft, col = c("#E3E3E3", "#51A3A3"), legend = FALSE, cex.axis = 2, las = 1, xlim=c(-145,-52))
plot(botmid, col = c("#E3E3E3", "#51A3A3"), legend = FALSE, yaxt = "n", cex.axis = 2, xlim=c(-145,-52))
plot(botright, col = c("#E3E3E3", "#51A3A3"), legend = FALSE, yaxt = "n", cex.axis = 2, xlim=c(-145,-52))
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

#reorganized AllLandcover csv and reclassified landcover classes to names
#making the figure
#calculating percent cover of each landcover class
landdat<-read.csv("prioritizrResults/landcover/AllLandcover.csv")
keeps <-landdat %>% group_by(landcover) %>% summarize(n=min(percent)) %>% filter(n>1)
landdat2 <- landdat %>% filter(landcover %in% keeps$landcover)

#calculating standard error for barplot
error <-landdat2 %>% group_by(landcover, Climate) %>% summarise(mean = mean(percent), se=sd(percent)/sqrt(3)) 

#making plot
ggplot(data=error, aes(x=reorder(landcover, -mean), y=mean, fill=Climate))+
  geom_bar(stat="summary",  color="black", position=position_dodge())  +
  geom_errorbar(data=error, aes(x=reorder(landcover, -mean), ymin=mean-se, ymax=mean+se), width=0, position=position_dodge(width=0.9)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_fill_manual(values=c("#FAFF81", "#FC7753", "#403D58"))+
  xlab("Landcover Class")+ ylab("Percent Cover")+scale_y_continuous(expand=c(0,0))


#####################
###Climate Data Prepping
#####################

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

###############
#Prioritizr data prepping
###############

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

