library(meteR)
library(plyr)
library(sp)
library(rgdal)

setwd('~/Dropbox/Research/hawaiiMETE_plant')

## read raw data
dat <- read.csv('~/Dropbox/Research/data/hawaiiPlant/Hawaii_PlantPlots_22022017.csv')

## organize plot info
datInfo <- ddply(dat, 'PlotIDn', function(x) 
    unique(x[, c('geo_entity', 'Study', 'Long_Dec', 'Lat_Dec', 'max_PlotArea', 'Elev', 
                 'MAP', 'MAT', 'PET', 'HII', 'HFP', 'PrecipSeasonality', 'TempSeasonality', 
                 'Age_FromRange', 'Age_FromCode', 'Prop_Introduced')]))

## turn into spatial object
datInfo <- SpatialPointsDataFrame(coords = datInfo[, c('Long_Dec', 'Lat_Dec')], data = datInfo, 
                                  proj4string = 
                                      CRS('+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs'))

## add missing geologic ages
oldwd <- setwd('~/Dropbox/hawaiiDimensions/geodata/env_data/geol')
geoPoly <- readOGR('hawaii_state_geol_ageClean.shp', 'hawaii_state_geol_ageClean')
setwd(oldwd)

## note that 'Age_FromRange' and 'age_mid' are the same except 'age_mid' has more non-NA values
datInfo <- spTransform(datInfo, CRS(proj4string(geoPoly)))
datInfo@data <- cbind(datInfo@data, 
                      over(datInfo, geoPoly)[, c('age_min', 'age_max', 'age_mid')])

## write out the plot info
save(datInfo, file = 'datInfo.RData')
