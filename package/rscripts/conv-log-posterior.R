

rm(list=ls( ))
source('src/myfunctions.R')


## Path to the OCTAVE/MATLAB output 
mcmcpath <- '../data/hap-tribars-nIndiv300-nSites3000-gridSize12x8-g12x8'
plotpath <- './hap-tribars-nIndiv300-nSites3000-gridSize12x8-g12x8'
simnos <- 1:3


pdf(height=7,width=14,file=paste(plotpath,'-lpost.pdf',sep=''))
plot.log.posterior(mcmcpath,simnos)
dev.off( )
