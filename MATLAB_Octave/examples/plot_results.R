

source('../rscripts/default.eems.plots.R')


mcmcpath <-
    paste('./data/barrier-schemeZ-nIndiv300-nSites3000-g13x7-simno',1:3,sep='')
plotpath <- './data/barrier-schemeZ-nIndiv300-nSites3000-g13x7-simno1_3'


## mcmcpath is a list of three output directories; the results are averaged
eemsplots(mcmcpath,plotpath,add.map=FALSE)

