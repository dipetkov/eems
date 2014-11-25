

source('default.eems.plots.R')


mcmcpath <- paste('../data/sat-barrier-nIndiv150-nSites16-EEMS-nDemes200-simno',1:3,sep='')
plotpath <- 'sat-barrier-nIndiv150-nSites16-EEMS-nDemes200-simno1_3'
 

## mcmcpath is actually a list of 3 `mcmcpath`s and the surfaces are averaged
eemsplots(mcmcpath,plotpath,add.map=FALSE)
