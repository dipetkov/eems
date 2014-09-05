

rm(list=ls( ))
source('src/myeems.plots.R')


## Path to the OCTAVE/MATLAB output 
datapath <- '../examples/data/uniform-schemeZ-nIndiv300-s12x8-u4Nm1-L3000'
mcmcpath <- '../examples/data/uniform-schemeZ-nIndiv300-s12x8-u4Nm1-L3000-g13x7'
plotpath <- '../examples/data/uniform-schemeZ-nIndiv300-s12x8-u4Nm1-L3000-g13x7'

## Ideally, we have several realizations of the Markov chain
simnos <- 1

dimns <- read.dimns(datapath)
plot.height <- 5
plot.width <- 7

png(file=paste(plotpath,'-simno',min(simnos),'_',max(simnos),'-pilogl%02d.png',sep=''),
    height=plot.height,width=plot.width,units="in",res=150)
plot.logposterior.simnos(mcmcpath,simnos)
dev.off( )
