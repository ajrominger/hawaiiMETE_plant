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


## remove special characters from names
plots$species <- gsub('[^[:alnum: ]]', '', plots$species)

## remove unwanted spaced
substring(plots$species, nchar(plots$species))[
    grep('[^[:alnum:]]', substring(plots$species, nchar(plots$species)))] <- ''

## capitalize first letter of all names
substring(plots$species, 1, 1) <- toupper(substring(plots$species, 1, 1))

## get native/non-native
sppList <- read.csv(file.path(dataWD, 'Hawaii_SPPLIST.csv'), as.is = TRUE)

plots$native <- sppList$Native_Status_Flora_HawaiianIslands_Synth[
    match(tolower(gsub(' ', '', plots$species)), 
          tolower(gsub(' ', '', sppList$Accepted_name_species)))]

## by hand native and spelling correction
plots$native[grep('Ageratina adenophora', plots$species)] <- 'Introduced'
plots$species[grep('Alectryon', plots$species)] <- 'Alectryon macrococcus'
plots$native[grep('Alectryon macrococcus', plots$species)] <- 'Native_Hawaii'
plots$species[grep('Alectryon macrococcus', plots$species)] <- 'Alectryon macrococcus'
plots$native[grep('Alyxia', plots$species)] <- 'Native_Hawaii'
plots$native[grep('Antidesma', plots$species)] <- 'Native_Hawaii'
plots$species[grep('Antidesma', plots$species)] <- 'Antidesma platyphyllum'
plots$native[grep('Cheirodendron', plots$species)] <- 'Native_Hawaii'
plots$species[grep('Claoxylon', plots$species)] <- 'Claoxylon sandwicense'
plots$native[grep('Claoxylon sandwicense', plots$species)] <- 'Native_Hawaii'
plots$native[grep('Clermontia', plots$species)] <- 'Native_Hawaii'
plots$native[grep('Coprosma', plots$species)] <- 'Native_Hawaii'
plots$species[grep('ordyline fruticosa', plots$species)] <- 'Cordyline fruticosa'
plots$native[grep('Cyanea', plots$species)] <- 'Native_Hawaii'
plots$native[grep('Cyrtandra', plots$species)] <- 'Native_Hawaii'
plots$species[grepl('Cyrtandra', plots$species) & grepl('longifolia', plots$species)] <- 
    'Cyrtandra longifolia'
plots$native[grep('Dubautia', plots$species)] <- 'Native_Hawaii'
plots$native[grep('Fig unknown', plots$species)] <- 'Introduced'
plots$native[grep('Flindersia brayleyana', plots$species)] <- 'Introduced'
plots$native[grep('Kadua', plots$species)] <- 'Native_Hawaii'
plots$native[grep('Labordia tinifolia', plots$species)] <- 'Native_Hawaii'
plots$species[grep('Leptocophy', plots$species)] <- 'Leptecophylla tameiameiae'
plots$native[grep('Leptecophylla tameiameiae', plots$species)] <- 'Native_Hawaii'
plots$native[grep('Lysiloma watsonii', plots$species)] <- 'Introduced'
plots$native[grep('Melastoma', plots$species)] <- 'Introduced'
plots$native[grep('Melestome', plots$species)] <- 'Introduced'
plots$native[grep('melicope', plots$species, ignore.case = TRUE)] <- 'Native_Hawaii'
plots$species[grep('clusiifolia', plots$species)] <- 'Melicope clusiifolia'
plots$native[grep('Myrsine', plots$species)] <- 'Native_Hawaii'
plots$species[grep('lessertiana', plots$species)] <- 'Myrsine lessertiana'
plots$species[grep('Nestegis', plots$species)] <- 'Nestegis sandwicensis'
plots$native[grep('Nestegis sandwicensis', plots$species)] <- 'Native_Hawaii'
plots$native[grep('Pinus', plots$species)] <- 'Introduced'
plots$native[grep('Platydesma spathulata', plots$species)] <- 'Native_Hawaii'
plots$native[grep('Psychotria', plots$species)] <- 'Native_Hawaii'
plots$native[grep('Rubus niveus', plots$species)] <- 'Introduced'
plots$native[grep('Rubus rosifolius', plots$species)] <- 'Introduced'
plots$native[grep('Scaevola gaudichaidiana', plots$species)] <- 'Native_Hawaii'
plots$native[grep('Smilax melastomifolia', plots$species)] <- 'Native_Hawaii'
plots$species[grep('Sophora', plots$species)] <- 'Sophora chrysophylla'
plots$native[grep('Sophora chrysophylla', plots$species)] <- 'Native_Hawaii'
plots$native[grep('Stachytarpheta dichotoma', plots$species)] <- 'Introduced'
plots$native[grep('Tibouchina herbacea', plots$species)] <- 'Introduced'
plots$species[grep('vaccinium reticulatum', plots$species)] <- 'Vaccinium reticulatum'
plots$species[grep('Vaccium reticulatum', plots$species)] <- 'Vaccinium reticulatum'
plots$species[grep('Vaccinum dentatum', plots$species)] <- 'Vaccinium dentatum'
plots$native[grep('Vaccinium', plots$species)] <- 'Native_Hawaii'
plots$native[grep('Xylosma hawaiiense', plots$species)] <- 'Native_Hawaii'
plots$species[grepl('Kadua', plots$species) & grepl('affinis', plots$species)] <- 
    'Kadua affinis'
plots$species[grepl('Coprosma', plots$species) & grepl('ochracea', plots$species)] <- 
    'Coprosma ochracea'

## make native binary
foo <- numeric(nrow(plots))
foo[plots$native == 'Native_Hawaii'] <- 1
foo[plots$native == 'Introduced'] <- 0
foo[is.na(plots$native)] <- NA
plots$native <- foo


## write it all out
save(plotsInfo, file = 'plotsInfo.RData')
write.csv(plots[, 1:(ncol(plots)-5)], file = file.path(dataWD, 'KnightBarton_clean.csv'), 
          row.names = FALSE)
write.csv(plots, file = file.path(dataWD, 'KnightBarton_flat.csv'), 
          row.names = FALSE)
