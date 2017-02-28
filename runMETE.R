library(meteR)
library(sp)
library(rgdal)
library(parallel)

setwd('~/Dropbox/Research/hawaiiMETE_plant')

## read data
dat <- read.csv('~/Dropbox/Research/data/hawaiiPlant/Hawaii_PlantPlots_22022017.csv')
load('datInfo.RData')

## add metabolic rate
dat$metab <- dat$DBH^2

## get METE predictions for SAD and IPD
meteSADIPD <- mclapply(split(dat[, c('PlotIDn', 'ID', 'SPP_CODE2', 'Abundance', 'metab')], 
                             dat$PlotIDn), 
                       mc.cores = 6,
                       FUN = function(x) {
                           if(length(unique(x$SPP_CODE2)) < 5) {
                               state.var <- c(S0 = length(unique(x$SPP_CODE2)), 
                                              N0 = sum(x$Abundance), 
                                              E0 = sum(x$metab))
                               return(list(sad = list(state.var = state.var),
                                           ipd = list(state.var = state.var)))
                           } else {
                               esf <- try(meteESF(spp = x$SPP_CODE2, abund = x$Abundance, 
                                                  power = x$metab))
                               if(class(esf) == 'try-error') {
                                   state.var <- c(S0 = length(unique(x$SPP_CODE2)), 
                                                  N0 = sum(x$Abundance), 
                                                  E0 = sum(x$metab))
                                   return(list(sad = list(state.var = state.var), 
                                               ipd = list(state.var = state.var)))
                               } else {
                                   thisSAD <- try(sad(esf))
                                   thisIPD <- try(ipd(esf))
                                   if(class(thisSAD)[1] == 'try-error') {
                                       thisSAD <- list(state.var = esf$state.var)
                                   }
                                   if(class(thisIPD)[1] == 'try-error') {
                                       thisIPD <- list(state.var = esf$state.var)
                                   }
                                   return(list(sad = thisSAD, ipd = thisIPD))
                               }
                           }
                       })

## proportion failures
sum(sapply(meteSADIPD, function(x) class(x[[1]])[1] == 'list')) / length(unique(dat$PlotIDn))

## extract mse
mseSADIPD <- t(sapply(meteSADIPD, function(x) {
    if(class(x$sad)[1] == 'list') {
        mseSAD <- NA
    } else {
        mseSAD <- mse(x$sad, relative = FALSE)
    }
    
    if(class(x$ipd)[1] == 'list') {
        mseIPD <- NA
    } else {
        mseIPD <- mse(x$ipd, relative = FALSE)
    }
    
    return(c(mseSAD = mseSAD, mseIPD = mseIPD))
}))

## extract R2

R2 <- function(x) {
    pred <- meteDist2Rank(x)
    obs <- x$data
    
    return(1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2))
}

r2SADIPD <- t(sapply(meteSADIPD, function(x) {
    if(class(x$sad)[1] == 'list') {
        r2SAD <- NA
    } else {
        r2SAD <- R2(x$sad)
    }
    
    if(class(x$ipd)[1] == 'list') {
        r2IPD <- NA
    } else {
        r2IPD <- R2(x$ipd)
    }
    
    return(c(r2SAD = r2SAD, r2IPD = r2IPD))
}))


## log likelihood z^2 values
z2SADIP <- t(sapply(meteSADIPD, function(x) {
    if(class(x$sad)[1] == 'list') {
        z2SAD <- NA
    } else {
        z2SAD <- logLikZ(x$sad)$z
    }
    
    if(class(x$ipd)[1] == 'list') {
        z2IPD <- NA
    } else {
        z2IPD <- NA
    }
    
    return(c(z2SAD = z2SAD, z2IPD = z2IPD))
}))


## extract state vars
stateVars <- t(sapply(meteSADIPD, function(x) x[[1]]$state.var))

## combine it all
meteSumm <- datInfo
meteSumm@data <- cbind(meteSumm@data, mseSADIPD, r2SADIPD, stateVars, z2SADIP)

save(meteSumm, file = 'meteSumm.RData')

## plotting
plot(meteSumm@data[, c('age_mid', 'z2SAD')], log = 'x')
plot(meteSumm@data[, c('MAP', 'z2SAD')], log = 'x')
plot(meteSumm@data[, c('MAT', 'z2SAD')], log = 'x')

fitCol <- function(x) {
    m <- 1/(max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
    
    newx <- m * (x - min(x, na.rm = TRUE))
    
    colval <- colorRamp(hsv(c(0.17, 0.6), c(0.5, 1), c(1, 0.7)))(newx)
    colval[is.na(colval)] <- 0
    
    col <- rgb(colval, maxColorValue = 255)
    col[col == '#000000'] <- 'transparent'
    
    return(col)
}

layout(matrix(2:1, nrow = 2), heights = c(1, 4))
par(mar = c(3, 3, 0, 0) + 0.5, mgp = c(2, 0.75, 0))
plot(meteSumm@data[, c('age_mid', 'MAP')], log = 'x', pch = 16, cex = 0.5, 
     col = fitCol(sqrt(meteSumm@data$z2SAD)))
par(mar = c(0, 3, 0, 0) + 0.5, mgp = c(2, 0.75, 0))
plot(sort(meteSumm@data$z2SAD), col = fitCol(sort(sqrt(meteSumm@data$z2SAD))), 
     pch = 16, cex = 0.5, ylab = 'z^2 value', xaxt = 'n')



plot(meteSumm, col = fitCol(sqrt(meteSumm@data$z2SAD)), pch = 16, cex = 0.5)


plot(meteSumm@data[, c('S0', 'z2SAD')])
plot(meteSumm@data[, c('N0', 'z2SAD')], log = 'x')
plot(meteSumm@data[, c('E0', 'z2SAD')], log = 'x')
plot(meteSumm@data[, c('Prop_Introduced', 'z2SAD')])

pairs(meteSumm@data[, c('S0', 'N0', 'E0')], log = 'xy')

plot(meteSumm@data[, c('max_PlotArea', 'z2SAD')], log = 'x')
plot(meteSumm@data[, c('max_PlotArea', 'S0')], log = 'xy')
plot(meteSumm@data[, c('max_PlotArea', 'N0')], log = 'xy')
plot(meteSumm@data[, c('max_PlotArea', 'E0')], log = 'xy')

plot(meteSumm@data[, c('age_mid', 'z2SAD')], log = 'x')
plot(meteSumm@data[, c('age_mid', 'S0')], log = 'xy')
plot(meteSumm@data[, c('age_mid', 'N0')], log = 'xy')
plot(meteSumm@data[, c('age_mid', 'E0')], log = 'xy')


layout(matrix(2:1, nrow = 2), heights = c(1, 4))
par(mar = c(3, 3, 0, 0) + 0.5, mgp = c(2, 0.75, 0))
plot(meteSumm@data[, c('age_mid', 'max_PlotArea')], log = 'xy', pch = 16, cex = 1, 
     col = fitCol(sqrt(meteSumm@data$z2SAD)))
par(mar = c(0, 3, 0, 0) + 0.5, mgp = c(2, 0.75, 0))
plot(sort(meteSumm@data$z2SAD), col = fitCol(sort(sqrt(meteSumm@data$z2SAD))), 
     pch = 16, cex = 0.5, ylab = 'z^2 value', xaxt = 'n')

layout(matrix(2:1, nrow = 2), heights = c(1, 4))
par(mar = c(3, 3, 0, 0) + 0.5, mgp = c(2, 0.75, 0))
plot(meteSumm@data[, c('age_mid', 'max_PlotArea')], log = 'xy', pch = 16, cex = 1, 
     col = fitCol(log(meteSumm@data$S0)))
par(mar = c(0, 3, 0, 0) + 0.5, mgp = c(2, 0.75, 0))
plot(sort(meteSumm@data$S0), col = fitCol(sort(log(meteSumm@data$S0))), 
     pch = 16, cex = 0.5, ylab = 'S0', xaxt = 'n', log = 'y')


par(mfrow = c(4, 4), mar = c(2, 2, 0, 0), oma = c(2, 2, 0.5, 0.5), 
    mgp = c(1.5, 0.3, 0), tcl = -0.25)
for(i in sample(which(!is.na(meteSumm@data$z2SAD)), prod(par('mfrow')))) {
    plot(meteSADIPD[[i]]$sad, ptype = 'rad', add.legend = FALSE, 
         xlab = '', ylab = '', log = 'y', ylim = c(1, max(meteSADIPD[[i]]$sad$data)))
}
mtext('Rank', side = 1, outer = TRUE)
mtext('Abundance', side = 2, outer = TRUE)
