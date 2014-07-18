

rm(list=ls( ))
source('src/myeems.plots.R')


## Path to the OCTAVE/MATLAB output 
datapath <- '../data/hap-barrier-nIndiv300-nSites3000-gridSize12x8'
mcmcpath <- '../data/hap-barrier-nIndiv300-nSites3000-gridSize12x8-g12x8-simno1'
plotpath <- './hap-barrier-nIndiv300-nSites3000-gridSize12x8-g12x8-simno1'


dimns <- read.dimns(datapath)
plot.height <- 5
plot.width <- 6
plot.filename <- paste(plotpath,'-rdist.png',sep='')

png(file=plot.filename,height=plot.height,width=plot.width,units="in",res=150)
dist.scatterplot(mcmcpath)
dev.off( )
