############################################################
#create event data from the isoseismal map of Al-Tarazi 2000
############################################################
library(rgdal)
library(raster)
library(sp)

setwd("/home/mhaas/PhD/oq-deserve/seismicity/")
crsys <- crs("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

rp475y <- read.csv("smoothed_seismicity.csv",sep=',')

coordinates(rp475y)=~lon+lat
coordinates(rp2475y)=~lon+lat

gridded(rp475y) = TRUE
gridded(rp2475y) = TRUE

int475a <- raster(rp475y)
int2475a <- raster(rp2475y)
projection(int475a) = crsys
projection(int2475a) = crsys
writeRaster(int475a,"int475a.tif")
writeRaster(int2475a,"int2475a.tif")