library(meteR)
library(parallel)
library(socorro)

setwd('~/Dropbox/Research/hawaiiMETE_plant')

## read data
dat <- read.csv('~/Dropbox/Research/data/hawaiiPlant/KnightBarton_clean.csv')
load('plotsInfo.RData')

meteKB <- lapply(split(dat, dat$site), 
                   # mc.cores = 6, 
                   FUN = function(x) {
                       esf <- try(meteESF(x$species, rep(1, nrow(x), x$dbh^2)))
                       
                       if(class(esf) == 'try-error') {
                           sadZ <- NA
                           ipdMSE <- NA
                       } else {
                           sadZ <- try(logLikZ(sad(esf), nrep = 2)$z)
                           if(class(sadZ) == 'try-error') sadZ <- NA
                           ipdMSE <- try(mse(ipd(esf)))
                           if(class(ipdMSE) == 'try-error') ipdMSE <- NA
                       }
                       
                       return(c(sadZ = sadZ, ipdMSE = ipdMSE))
                   })

meteKB <- do.call(rbind, meteKB)
meteKB <- data.frame(site = rownames(meteKB), meteKB)
rownames(meteKB) <- NULL

meteKB <- cbind(meteKB, plotsInfo@data[match(meteKB$site, plotsInfo@data$site), 
                                       c('age_mid', 'MAP', 'elev')])
                
plot(meteKB$age_mid, meteKB$sadZ)
plot(meteKB$elev, meteKB$sadZ)
plot(meteKB$MAP, meteKB$sadZ)
plot(meteKB[!is.na(meteKB$sadZ), c('age_mid', 'elev')])


fitCol <- function(x) {
    m <- 1/(max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
    
    newx <- m * (x - min(x, na.rm = TRUE))
    
    colval <- colorRamp(hsv(c(0.17, 0.6), c(0.5, 1), c(1, 0.7)))(newx)
    colval[is.na(colval)] <- 0
    
    col <- rgb(colval, maxColorValue = 255)
    col[col == '#000000'] <- 'transparent'
    
    return(col)
}

plot(meteKB[, c('age_mid', 'elev')], col = fitCol(meteKB$sadZ), pch = 16, cex = 2)
