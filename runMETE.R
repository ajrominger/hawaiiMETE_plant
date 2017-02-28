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
z2SADIP <- mclapply(meteSADIPD, mc.cores = 6, FUN = function(x) {
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
})

z2SADIP <- do.call(rbind, z2SADIP)

## extract state vars
stateVars <- t(sapply(meteSADIPD, function(x) x[[1]]$state.var))

## combine it all
meteSumm <- datInfo
meteSumm@data <- cbind(meteSumm@data, mseSADIPD, r2SADIPD, stateVars, z2SADIP)

## save output
save(meteSumm, file = 'meteSumm.RData')

## save a few example SAD and IPDs
meteEg <- meteSADIPD[sample(which(!is.na(meteSumm@data$z2SAD)), 16)]
save(meteEg, file = 'meteEg.RData')
