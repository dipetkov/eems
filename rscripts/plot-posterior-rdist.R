

rm(list=ls( ))
source('src/myeems.plots.R')


datapath <- '../examples/data/barrier-schemeZ-nIndiv300-s12x8-u4Nm1-L3000'
mcmcpath <- '../examples/data/barrier-schemeZ-nIndiv300-s12x8-u4Nm1-L3000-g13x7-simno1'
plotpath <- '../examples/data/barrier-schemeZ-nIndiv300-s12x8-u4Nm1-L3000-g13x7-simno1'


plot.height <- 5
plot.width <- 6

png(file=paste(plotpath,'-rdist%02d.png',sep=''),
    height=plot.height,width=plot.width,units="in",res=150)
dist.scatterplot(mcmcpath)
dev.off( )
