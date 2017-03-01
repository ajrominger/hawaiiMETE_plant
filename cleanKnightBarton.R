library(meteR)
library(plyr)
library(sp)
library(rgdal)
library(xlsx)

setwd('~/Dropbox/Research/hawaiiMETE_plant')

dataWD <- '~/Dropbox/Research/data/hawaiiPlant'
files <- c('big island.xlsx', 'maui.xlsx', 'molokai.xlsx', 'oahu_final.xlsx', 'kauai.xlsx')

plotCenters <- sapply(files[files != 'molokai.xlsx'], function(f) {
    wb <- loadWorkbook(file.path(dataWD, f))
    sheets <- names(getSheets(wb))
    sheets <- sheets[!grepl('list', sheets)]
    
    t(sapply(sheets, function(s) {
        x <- read.xlsx2(file.path(dataWD, f), sheetName = s, header = TRUE)
        
        if('Easting' %in% names(x)) {
            x$Easting <- as.numeric(as.character(x$Easting))
            x$Northing <- as.numeric(as.character(x$Northing))
            x$dbh <- as.numeric(as.character(x$dbh))
            x <- x[!is.na(x$Easting) & !is.na(x$Northing), ]
            if(class(x) == 'try-error') browser()
            
            xsp <- SpatialPointsDataFrame(coords = x[, 1:2], data = x,
                                          proj4string = CRS('+proj=utm +zone=5 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0'))
            
            # xsp <- spTransform(xsp, CRS(proj4string(islands)))
            return(colMeans(coordinates(xsp)))
        } else {
            return(rep(NA, 2))
        }
    }))
})

plotCenters <- do.call(rbind, plotCenters)
plotCenters <- data.frame(plotCenters, name = rownames(plotCenters))
rownames(plotCenters) <- NULL
plotCenters <- plotCenters[!is.na(plotCenters$Easting), ]
plotCenters <- SpatialPointsDataFrame(plotCenters[, 1:2], plotCenters[, -(1:2), drop = FALSE], 
                                      proj4string = CRS('+proj=utm +zone=5 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0'))

writeOGR(plotCenters, 'plotCenters', 'plotCenters', driver = 'ESRI Shapefile', 
         overwrite_layer = TRUE)
