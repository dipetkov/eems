

rm(list=ls( ))
source('src/myeems.plots.R')


## Path to the OCTAVE/MATLAB output 
datapath <- '../examples/data/uniform-schemeZ-nIndiv300-s12x8-u4Nm1-L3000'
mcmcpath <- '../examples/data/uniform-schemeZ-nIndiv300-s12x8-u4Nm1-L3000-g13x7-simno1'
plotpath <- '../examples/data/uniform-schemeZ-nIndiv300-s12x8-u4Nm1-L3000-g13x7-simno1'


dimns <- read.dimns(datapath)
plot.height <- 5
plot.width <- 5*(dimns$xspan/dimns$yspan)
plot.filename <- 

png(file=paste(plotpath,'-mVoronoi%03d.png',sep=''),
    height=plot.height,width=plot.width,units="in",res=150)
mlegend <- mcmc.mrates.voronoi(mcmcpath,dimns)
dev.off( )
png(file=paste(plotpath,'-qVoronoi%03d.png',sep=''),
    height=plot.height,width=plot.width,units="in",res=150)
mlegend <- mcmc.qrates.voronoi(mcmcpath,dimns)
dev.off( )
