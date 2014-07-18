

rm(list=ls( ))
source('src/myeems.plots.R')


## Path to the OCTAVE/MATLAB output 
datapath <- '../data/dip-uniform-nIndiv100-nSites1000-gridSize8x7'
mcmcpath <- '../data/dip-uniform-nIndiv100-nSites1000-gridSize8x7-g12x8'
plotpath <- './dip-uniform-nIndiv100-nSites1000-gridSize8x7-g12x8'

simnos <- 1:3

dimns <- read.dimns(datapath)
plot.height <- 5
plot.width <- 5*(dimns$xspan/dimns$yspan)
plot.filename <- paste(plotpath,'-simno',min(simnos),'_',max(simnos),'-avesoln.png',sep='')

png(file=plot.filename,height=plot.height,width=plot.width,units="in",res=150)
mcmc.mrates.simnos(mcmcpath,simnos,dimns)
dev.off( )
