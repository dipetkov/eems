

rm(list=ls( ))
source('src/myeems.plots.R')


## Path to the OCTAVE/MATLAB output 
datapath <- '../examples/data/barrier-schemeZ-nIndiv300-s12x8-u4Nm1-L3000'
mcmcpath <- '../examples/data/barrier-schemeZ-nIndiv300-s12x8-u4Nm1-L3000-g13x7-simno1'
plotpath <- '../examples/data/barrier-schemeZ-nIndiv300-s12x8-u4Nm1-L3000-g13x7-simno1'


dimns <- read.dimns(datapath)
plot.height <- 5
plot.width <- 5*(dimns$xspan/dimns$yspan)

bitmap(paste(plotpath,'-mVoronoi%03d.png',sep=''),type='png16m',res=300,
       height=plot.height,width=plot.width,units='in')
mlegend <- mcmc.mrates.voronoi(mcmcpath,dimns)
dev.off( )
bitmap(paste(plotpath,'-qVoronoi%03d.png',sep=''),type='png16m',res=300,
       height=plot.height,width=plot.width,units='in')
mlegend <- mcmc.qrates.voronoi(mcmcpath,dimns)
dev.off( )
