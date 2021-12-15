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
  
  ## load buffer - least polygon extent
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

  
  #species area for current climate
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


########################################################################
#Prioritizr analyses
########################################################################

#libraries
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
currentstackmin[currentstackmin<0.2]<-0
writeRaster(currentstackmin, "Outputs//CurrentStack.grd", overwrite=T)

future26list<-list.files("Outputs//future26",full.names=T)
future26stack<-stack(future26list)
future26stackmin<-future26stack
future26stackmin<-crop(future26stackmin, canada1)
future26stackmin[canada1==1 & is.na(future26stackmin)] <- 0
future26stackmin<-mask(future26stackmin, polyCAN) 
future26stackmin[future26stackmin<0.2]<-0
writeRaster(future26stackmin, "Outputs//Future26Stack.grd", overwrite=T)

future85list<-list.files("Outputs//future85",full.names=T)
future85stack<-stack(future85list)
future85stackmin<-future85stack
future85stackmin<-crop(future85stackmin, canada1)
future85stackmin[canada1==1 & is.na(future85stackmin)] <- 0
future85stackmin<-mask(future85stackmin, polyCAN) 
future85stackmin[future85stackmin<0.2]<-0
writeRaster(future85stackmin, "Outputs//Future85Stack.grd", overwrite=T)

###datafiles
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

currentstackless<-currentstack[[-c(11, 16, 17, 32)]] #dropping crotchii, fraternus, frigidus, polaris

## add min set objective function

bb_priormin <- function(pu = canada, features = currentstack, rel_tar = 0.17, 
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

####Min set

### Problem 1a add_min_set objective, 17%, current climate
p1a <- bb_priormin(raster_name =  "prioritizrResults//s1aMinSet17.tif")

p1a2 <- bb_priormin(features = currentstackless, raster_name =  "prioritizrResults//s1a2MinSet17.tif")
### Problem 2a add_min_set objective, 30% (Canada 2050 biodiversity target), current climate
p2a <- bb_priormin( rel_tar = 0.3, raster_name =  "prioritizrResults//s2aMinSet30.tif")

p2a2<-bb_priormin( features= currentstackless, rel_tar = 0.3, raster_name =  "prioritizrResults//s2a2MinSet30.tif")

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

###################
#Prioritizr results figures
###################

## prioritizr results figure 1 add min set objective with problematic speceis removed
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
#calculating weighted  percent cover of each landcover class
#weighted by proprotion of each landcover type throughout Canada
landdat<-read.csv("prioritizrResults/landcover/AllLandcover.csv")
keeps <-landdat %>% group_by(landcover) %>% summarize(n=min(weightedPercent)) %>% filter(n>1)
landdat2 <- landdat %>% filter(landcover %in% keeps$landcover)

#calculating standard error for barplot
error <-landdat2 %>% group_by(landcover, climate) %>% summarise(mean = mean(weightedPercent), se=sd(weightedPercent)/sqrt(3)) 

#making plot
ggplot(data=error, aes(x=reorder(landcover, -mean), y=mean, fill=climate))+
  geom_bar(stat="summary",  color="black", position=position_dodge())  +
  geom_errorbar(data=error, aes(x=reorder(landcover, -mean), ymin=mean-se, ymax=mean+se), width=0, position=position_dodge(width=0.9)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_fill_manual(values=c("#FAFF81", "#FC7753", "#403D58"))+
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
error <-object %>% group_by(status, climate, target) %>% summarise(mean = mean(relative_held, na.rm=TRUE), se=sd(relative_held, na.rm=TRUE)/sqrt(length(relative_held[!is.na(relative_held)])))

ggplot(data=error, aes(x=climate, y=mean, fill=status))+
  geom_bar(stat="summary", color="black", position=position_dodge()) + facet_grid(cols=vars(target)) +
  geom_errorbar(data=error, aes (x=climate, ymin=mean-se, ymax=mean+se), width=0, position=position_dodge(width=0.9)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_fill_manual(values=c("#EAC435", "#345995"))+
  xlab("")+ ylab("Proportion of distribution")+scale_y_continuous(expand=c(0,0), breaks=seq(0, .70, .10))

#### percent area with >50% suitability figure
library(ggplot2)
library(gridExtra)
library(forcats)
library(raster)

#the species range did not calculate properly due to a typo in the function so recalculating here
currentrangestack <- stack("Outputs/CurrentUpdated.grd")
future26rangestack<-stack("Outputs/Future26Updated.grd")
future85rangestack<-stack("Outputs/Future85Updated.grd")
bumblebeenames <- c("affinis", "appositus", "auricomus", "bifarius", "bimaculatus", "bohemicus", 
                    "borealis", "caliginosus", "centralis", "citrinus", "crotchii", "cryptarum", 
                    "fervidus", "flavidus", "flavifrons", "fraternus", "frigidus", "griseocollis",
                    "huntii", "impatiens", "insularis", "jonellus", "kirbiellus", "kluanensis", 
                    "melanopygus", "mixtus", "morrisoni", "natvigi", "neoboreus", "nevadensis", 
                    "occidentalis", "pensylvanicus", "perplexus", "polaris", "rufocinctus", 
                    "sandersoni", "sitkensis", "suckleyi", "sylvicola", "ternarius", "terricola",
                    "vagans", "vandykei", "vosnesenskii")
names(currentrangestack) <- bumblebeenames

calcsppRange <- function(x){
  spprange<- x
  spprange[spprange <= 0.001] <- NA
  spprange<-spprange[!is.na(spprange)]
  spprange<- ncell(spprange)/ncell(x)
  spprange<-(spprange)*1717662
  return(spprange)
}
outRange <- lapply(1:nlayers(currentrangestack), function(i){
  calcsppRange(currentrangestack[[i]])
})
currentRanges <- data.frame(Species = bumblebeenames, CurrentRange= do.call(c, outRange))

outRange <- lapply(1:nlayers(future26rangestack), function(i){
  calcsppRange(future26rangestack[[i]])
})
future26Ranges <- data.frame(Species = bumblebeenames, future26Range= do.call(c, outRange))

outRange <- lapply(1:nlayers(future85rangestack), function(i){
  calcsppRange(future85rangestack[[i]])
})
future85Ranges <- data.frame(Species = bumblebeenames, future85Range= do.call(c, outRange))


prop50 <- read.csv("Outputs/allModels.csv", stringsAsFactors=FALSE)
library(tidyverse)
prop50 <- prop50 %>% filter(!(species %in% c("bohemicus", "natvigi","neoboreus", "polaris", "suckleyi"))) %>% 
  dplyr::select(species, Current= currentRange,Future26= future26Range, Future85=future85Range) %>% 
  mutate(diff26 = (Future26-Current)/Current, diff85 = (Future85-Current)/Current) %>% 
  mutate(species = factor(species, levels=species[order(diff26)] )) %>% ## re-order by change in area
  mutate(group = ifelse(diff26 < 0, "decreasing","increasing")) %>% 
  gather(climate, proportion, diff26:diff85) 


plot1 <- ggplot(prop50 %>% filter(group=="decreasing"), aes(x= fct_rev(species), y=proportion*100, fill=fct_rev(climate))) + geom_bar(stat="identity", position="dodge") + coord_flip() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_fill_manual(values=c("#403D58", "#FC7753", "#FAFF81")) +
  xlab("")+ ylab("Percent Change in Range from Current")+ ylim(-100,5)
  # scale_y_continuous(expand=c(0,0)) 

plot2 <- ggplot(prop50 %>% filter(group=="increasing"), aes(x= fct_rev(species), y=proportion*100, fill=fct_rev(climate))) + geom_bar(stat="identity", position="dodge") + coord_flip() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_fill_manual(values=c("#403D58", "#FC7753", "#FAFF81")) +
  xlab("")+ ylab("Percent Change in Range from Current")+ ylim(0,120)
  # scale_y_continuous(expand=c(0,0))

gridExtra::grid.arrange(plot1, plot2, ncol=2)

### change in latitude and elevation table
#libraries
library(raster)
library(dplyr)

#datafiles
elevation<-raster("CAN_msk_alt.grd")
s1a<-raster("prioritizrResults/s1a3MinSet17.tif")
s2a<-raster("prioritizrResults/s2a3MinSet30.tif")
s3a<-raster("prioritizrResults/s3a3MinSet50.tif")
s1b<-raster("prioritizrResults/s1b3MinSet17RCP26.tif")
s2b<-raster("prioritizrResults/s2b3MinSet17RCP26.tif")
s3b<-raster("prioritizrResults/s3b3Minset50RCP26.tif")
s1c<-raster("prioritizrResults/s1c3MinSet17RCP85.tif")
s2c<-raster("prioritizrResults/s2c3MinSet30RCP85.tif")
s3c<-raster("prioritizrResults/s3c3MinSet50RCP85.tif")

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




########################
###Climate Data Prepping
########################

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

