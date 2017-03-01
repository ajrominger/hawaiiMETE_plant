library(meteR)
library(parallel)
library(socorro)

setwd('~/Dropbox/Research/hawaiiMETE_plant')

## read data
dataAbund <- read.csv('~/Dropbox/Research/data/hawaiiPlant/KnightBarton_summAbund.csv')
# dataMetab <- read.csv('~/Dropbox/Research/data/hawaiiPlant/KnightBarton_summBasalArea.csv')

meteKB <- mclapply(1:nrow(dataAbund), 
                   mc.cores = 6, 
                   FUN = function(i) {
                       abund <- as.numeric(dataAbund[i, -(1:4)])
                       spp <- names(dataAbund)[-(1:4)]
                       
                       spp <- spp[abund > 0]
                       abund <- abund[abund > 0]
                       
                       out <- try(sad(meteESF(spp, abund)))
                       if(class(out) == 'try-error') {
                           return(list(state.var = c(S0 = length(spp),
                                                     N0 = sum(abund),
                                                     E0 = NA)))
                       } else {
                           return(out)
                       }
                   })


meteKBsumm <- mclapply(meteKB, mc.cores = 6, FUN = function(x) {
    if(class(x)[1] == 'list') {
        z2 <- NA
    } else {
        z2 <- logLikZ(x)$z
    }
    
    return(c(x$state.var, z2 = z2))
})

meteKBsumm <- do.call(rbind, meteKBsumm)
meteKBsumm <- cbind(dataAbund[, 1:3], meteKBsumm)

## plotting

pdf('fig_sadKnightBarton.pdf', width = 5, height = 6)

nplot <- tapply(!is.na(meteKBsumm$z2), meteKBsumm$Island, sum)
par(mfcol = c(max(nplot), length(nplot)), mar = c(1, 1, 0, 0), oma = c(1, 1, 3, 0) + 0.5)

for(i in c('Kauai', 'Oahu', 'Molokai', 'Maui', 'Big')) {
    for(j in 1:nplot[i]) {
        if(i == 'Kauai' & j == 1) {
            par(mfg = c(6, 1))
            plot(1, axes = FALSE, type = 'n')
        }
        par(mfg = c(max(nplot) + 1 - j, switch(i, 'Kauai' = 1, 
                              'Oahu' = 2,
                              'Molokai' = 3,
                              'Maui' = 4,
                              'Big' = 5)))
        plot(meteKB[[which(meteKBsumm$Island == i & !is.na(meteKBsumm$z2))[j]]], 
             ptype = 'rad', add.legend = FALSE, log = 'y', axes = FALSE,
             xlim = c(0.8, max(meteKBsumm$S0)), ylim = c(0.8, max(meteKBsumm$N0)))
        logAxis(2, col.axis = 'transparent')
        axis(1, col.axis = 'transparent')
        box()
        
        legend('topright', legend = round(meteKBsumm$z2[which(meteKBsumm$Island == i &
                                                                  !is.na(meteKBsumm$z2))[j]],
                                          3),
               bty = 'n')
    }
    
    par(xpd = NA, 
        mfg = c(1, switch(i, 'Kauai' = 1, 
                                           'Oahu' = 2,
                                           'Molokai' = 3,
                                           'Maui' = 4,
                                           'Big' = 5)))
    plot(1, type = 'n', xlab = '', ylab = '', axes = FALSE)
    mtext(i, side = 3, line = 1)
    par(xpd = FALSE)
}

mtext('Abundance', side = 2, outer = TRUE, line = 0)
mtext('Rank', side = 1, outer = TRUE, line = 0)

dev.off()


plot(meteKBsumm$S0, !is.na(meteKBsumm$z2))
