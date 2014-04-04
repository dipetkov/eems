

rm(list=ls( ))
source('src/myfunctions.R')


## Path to the OCTAVE/MATLAB output 
datapath <- '../data/hap-tribars-nIndiv300-nSites3000-gridSize12x8'
mcmcpath <- '..//data/hap-tribars-nIndiv300-nSites3000-gridSize12x8-g12x8-simno1'
plotpath <- './hap-tribars-nIndiv300-nSites3000-gridSize12x8-g12x8-simno1'


dimns <- read.dimns(datapath)

pdf(height=5,width=5*(dimns$xspan/dimns$yspan),
    file=paste(plotpath,'-draws.pdf',sep=''))
mcmc.voronoi(mcmcpath,dimns)

dev.off( )
