library(meteR)
library(plyr)
library(sp)
library(rgdal)
library(raster)
library(xlsx)

setwd('~/Dropbox/Research/hawaiiMETE_plant')

dataWD <- '~/Dropbox/Research/data/hawaiiPlant'
files <- c('big island.xlsx', 'maui.xlsx', 'molokai.xlsx', 'oahu_final.xlsx', 'kauai.xlsx')

plots <- lapply(files, function(f) {
    wb <- loadWorkbook(file.path(dataWD, f))
    sheets <- names(getSheets(wb))
    sheets <- sheets[!grepl('list', sheets)]
    
    out <- lapply(sheets, function(s) {
        print(paste(f, s, sep = ': '))
        
        x <- read.xlsx2(file.path(dataWD, f), sheetName = s, header = TRUE)
        
        ## insure there are no factors
        x <- as.data.frame(apply(x, 2, as.character), stringsAsFactors = FALSE)
        
        if(!('Easting' %in% names(x))) {
            ## need to convert lat lon to utm
            x$Longitude <- as.numeric(gsub('[^0-9\\.-]', '', x$Longitude))
            x$Latitude <- as.numeric(gsub('[^0-9\\.-]', '', x$Latitude))
            noCoord <- is.na(x$Longitude) | is.na(x$Latitude)
            
            xsp <- SpatialPoints(x[!noCoord, 1:2], proj4string = CRS('+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs'))
            xsp <- spTransform(xsp, CRS('+proj=utm +zone=5 +datum=NAD83 +units=m +no_defs +towgs84=0,0,0'))
            
            x[!noCoord, 1:2] <- coordinates(xsp)
            names(x)[1:2] <- c('Easting', 'Northing')
        }
        
        ## get correct column classes
        x$Easting <- as.numeric(x$Easting)
        x$Northing <- as.numeric(x$Northing)
        
        ## add island and site
        x <- cbind(island = gsub('.xlsx', '', f), site = s, x)
        
        ## could be missing common name column
        if(ncol(x) < 11) x$common.name <- NA
        
        return(x)
    })
    
    ## correct variability in names
    out <- lapply(out, function(x) {
        names(x) <- names(out[[1]])
        return(x)
    })
    
    out <- do.call(rbind, out)
    
    return(out)
})

## correct variability in names
plots <- lapply(plots, function(x) {
    names(x) <- names(plots[[1]])
    return(x)
})
plots <- do.call(rbind, plots)

## final column class correction
plots$dbh <- as.numeric(plots$dbh)
plots$dbh[plots$dbh == 999] <- NA
plots$Point_ID <- as.numeric(plots$Point_ID)


## get plot centroids

plotsInfo <- ddply(plots, 'site', function(x) data.frame(island = x$island[1], 
                                                         lon = mean(x$Easting, na.rm = TRUE), 
                                                         lat = mean(x$Northing, na.rm = TRUE)))
plotsInfo$island <- gsub('_.*', '', plotsInfo$island)
plotsInfo <- SpatialPointsDataFrame(plotsInfo[, c('lon', 'lat')], 
                                    plotsInfo[, c('island', 'site')], 
                                    proj4string = CRS('+proj=utm +zone=5 +datum=NAD83 +units=m +no_defs +towgs84=0,0,0'))

## match to geol data

envDataWD <- '~/Dropbox/hawaiiDimensions/geodata/env_data'

oldwd <- setwd(file.path(envDataWD, 'geol'))
geol <- readOGR('hawaii_state_geol_ageClean.shp', 'hawaii_state_geol_ageClean')
plotsInfo@data <- cbind(plotsInfo@data, over(spTransform(plotsInfo, CRS(proj4string(geol))), 
                                             geol)[, c('age_min', 'age_max', 'age_mid')])
setwd(oldwd)


## match to precip
precip <- raster(file.path(envDataWD, 'precip/StateRFGrids_mm2/staterf_mmann/w001001.adf'))
plotsInfo@data$MAP <- extract(precip, spTransform(plotsInfo, CRS(proj4string(precip))))

## match to elevation

## helper function
elev <- function(island) {
    origIsland <- island
    island <- gsub(' ', '', island)
    
    info <- readLines(sprintf(file.path(envDataWD, '%sDEM/README.%s'), island, island), n = 1)
    info <- strsplit(info, ', ')[[1]]
    
    thisCRS <- sprintf('+proj=utm +zone=%s +datum=%s', gsub('UTM zone ', '', info[2]), info[3])
    
    r <- raster(sprintf(file.path(envDataWD, '%sDEM/%s.bil'), island, island), crs = CRS(thisCRS))
    
    e <- extract(r, spTransform(plotsInfo[plotsInfo$island == origIsland, ], CRS(proj4string(r))))
    
    if(gsub('.* ', '', info[1]) == 'feet') e <- round(e*0.3048)
    
    plotsInfo@data$elev[plotsInfo$island == origIsland] <<- e
}

for(i in unique(plotsInfo@data$island)) elev(i)

## combine with plot data in flat format
plots$age_mid <- plotsInfo@data$age_mid[match(plots$site, plotsInfo@data$site)]
plots$MAP <- plotsInfo@data$MAP[match(plots$site, plotsInfo@data$site)]
plots$elev <- plotsInfo@data$elev[match(plots$site, plotsInfo@data$site)]
plots$Lon <- plots$Lat <- NA
plots[, c('Lon', 'Lat')] <- coordinates(spTransform(plotsInfo, 
                                                    CRS('+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs')))[match(plots$site, plotsInfo@data$site),]

## write it all out
save(plotsInfo, file = 'plotsInfo.RData')
write.csv(plots[, 1:(ncol(plots)-5)], file = file.path(dataWD, 'KnightBarton_clean.csv'), 
          row.names = FALSE)
write.csv(plots, file = file.path(dataWD, 'KnightBarton_flat.csv'), 
          row.names = FALSE)
