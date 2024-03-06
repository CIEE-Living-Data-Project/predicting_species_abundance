# add in elevation from worldclim for each lat long coordinate, MAT, MAP, Tseas
library(terra)
library(dplyr)


long.lat <- readRDS("Revision 1 ecography/output/prep_data/within.study.updated.interactions.020724ENB.RDS") 
long.lat <- distinct(select(long.lat,LATITUDE, LONGITUDE))

#### read in in elevation data #####
# saved locally
elev.tif <-'/Users/Gavia/Documents/14 U of T/Mahler lab/Ugrads/Michelle/Analysis/bioclim wrangling/wc2.1_30s_elev.tif' 
elev_SpatRaster=rast(elev.tif) # could also use raster but maybe best to stick with terra

##### read in Bioclim data ####
# bioclim data downloaded from https://www.worldclim.org/data/worldclim21.html, Jan 2022
# saved locally
rasts <- rast(list.files("/Users/Gavia/Documents/14 U of T/Mahler lab/Ugrads/Michelle/Analysis/bioclim wrangling/wc2.1_30s_bio", full.names=T)) # in terra
name1 <- as.character(names(rasts)) #name of bioclim rasters

##### Bioclim values #####

worldclim <- matrix(ncol=5, nrow=nrow(long.lat))
colnames(worldclim) <- c("Elevation", "MAT", "MAP", "Tseas", "Pseas")

# 54 and 57 are NAs bc long lats are in ocean, both tropical. fill with close by values
# lat long taken from point near ocean in Viaques National wildlife refuge on PR
# both NA studies in ocean off of wildlife refuge in PR
data.NA <- data.frame(LATIDUE=18.157528, LONGITUDE= -65.413876)

for(i in 1:nrow(long.lat)) {
  print(i)
  # match elevation to each long.lat occurrences
  worldclim[i,1] <- terra::extract(elev_SpatRaster, long.lat[i, c(2,1)])[[2]] # pulls out elevations at occurrence in a vector
  if(is.na(worldclim[i,1])==TRUE){
    worldclim[i,1] <- terra::extract(elev_SpatRaster, data.NA[, c(2,1)])[[2]]
  }
  # summarize bioclim variables
  worldclim[i,2] <- terra::extract(rasts[[1]], long.lat[i, c(2,1)])[[2]] # bio1
  if(is.nan(worldclim[i,2])==TRUE){
    worldclim[i,2] <- terra::extract(rasts[[1]], data.NA[, c(2,1)])[[2]]
  }
  worldclim[i,3]  <- terra::extract(rasts[[4]], long.lat[i, c(2,1)])[[2]] # bio12
  if(is.nan(worldclim[i,3])==TRUE){
    worldclim[i,3] <- terra::extract(rasts[[4]], data.NA[, c(2,1)])[[2]]
  }
  worldclim[i,4]  <- terra::extract(rasts[[14]], long.lat[i, c(2,1)])[[2]] # bio4
  if(is.nan(worldclim[i,4])==TRUE){
    worldclim[i,4] <- terra::extract(rasts[[14]], data.NA[, c(2,1)])[[2]]
  }
  worldclim[i,5]  <- terra::extract(rasts[[7]], long.lat[i, c(2,1)])[[2]] # bio15
  if(is.nan(worldclim[i,5])==TRUE){
    worldclim[i,5] <- terra::extract(rasts[[7]], data.NA[, c(2,1)])[[2]]
  }
}

worldclim <- data.frame(cbind(long.lat, worldclim))

# calculate correlations
# check correlations in data between MAP, MAT, elev, abs lat
# all low, oddly elev and abs.lat highest cor, but not such that can't include in same model
cor(worldclim$MAP, worldclim$MAT) #0.11254
cor(worldclim$MAP, worldclim$Elevation) #-0.2734882
cor(worldclim$MAT, worldclim$Elevation) #-0.09762336
cor(worldclim$MAT, worldclim$Tseas) #-0.5267236
cor(worldclim$MAT, worldclim$Pseas) #0.1461913

cor(worldclim$Tseas, worldclim$Pseas) #0.2121774
cor(worldclim$Tseas, worldclim$Elevation) # 0.1057463
cor(worldclim$Pseas, worldclim$Elevation) # 0.3247045

cor(abs(worldclim$LATITUDE), worldclim$MAT) #-0.8647094
cor(abs(worldclim$LATITUDE), worldclim$MAP) #-0.1770841
cor(abs(worldclim$LATITUDE), worldclim$Elevation) #-0.2307737
cor(abs(worldclim$LATITUDE), worldclim$Tseas) #0.2882644
cor(abs(worldclim$LATITUDE), worldclim$Pseas) #-0.2792891

  
write.csv(data.frame(cbind(long.lat, worldclim)), file="Revision 1 ecography/output/prep_data/worldclim.csv")
