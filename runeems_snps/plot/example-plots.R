

source('default.eems.plots.R')


mcmcpath <- paste('../data/barrier-schemeX-nIndiv300-nSites3000-EEMS-nDemes153-simno',1:3,sep='')
plotpath <- 'barrier-schemeX-nIndiv300-nSites3000-EEMS-nDemes153-simno1_3'


## mcmcpath is a list of three output directories; the results are averaged
eemsplots(mcmcpath,plotpath,add.map=FALSE)
