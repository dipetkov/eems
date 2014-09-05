

rm(list=ls( ))
source('src/myeems.plots.R')


## Path to the OCTAVE/MATLAB output 
datapath <- '../examples/data/uniform-schemeZ-nIndiv300-s12x8-u4Nm1-L3000'
mcmcpath <- '../examples/data/uniform-schemeZ-nIndiv300-s12x8-u4Nm1-L3000-g13x7'
plotpath <- '../examples/data/uniform-schemeZ-nIndiv300-s12x8-u4Nm1-L3000-g13x7'

simnos <- 1:5

dimns <- read.dimns(datapath)
plot.height <- 5
plot.width <- 5*(dimns$xspan/dimns$yspan)

png(file=paste(plotpath,'-simno',min(simnos),'_',max(simnos),'-avesoln%02d.png',sep=''),
    height=plot.height,width=plot.width,units="in",res=150)
mlegend <- mcmc.mrates.simnos(mcmcpath,simnos,dimns)
qlegend <- mcmc.qrates.simnos(mcmcpath,simnos,dimns)
dev.off( )

png(file=paste(plotpath,'-simno',min(simnos),'_',max(simnos),'-legend%02d.png',sep=''),
    height=plot.height,width=0.3*plot.width,units="in",res=150)
mcmc.mrates.legend(datapath,mcmcpath,dimns,mmrks,mlegend)
mcmc.qrates.legend(datapath,mcmcpath,dimns,qmrks,qlegend)
dev.off( )
