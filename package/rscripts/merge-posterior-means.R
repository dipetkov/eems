

rm(list=ls( ))
source('src/myfunctions.R')


## Path to the OCTAVE/MATLAB output 
datapath <- '../data/hap-tribars-nIndiv300-nSites3000-gridSize12x8'
mcmcpath <- '../data/hap-tribars-nIndiv300-nSites3000-gridSize12x8-g12x8'
plotpath <- './hap-tribars-nIndiv300-nSites3000-gridSize12x8-g12x8-simno1'

simnos <- 1:3

nxmrks <- 50
nymrks <- 38
aveZmrks <- matrix(0,nxmrks,nymrks)
dimns <- read.dimns(datapath,nxmrks=nxmrks,nymrks=nymrks)
count <- 0


for (simno in simnos) {
    Zmrks <- mcmc.solution.pt1(paste(mcmcpath,'-simno',simno,sep=''),dimns)
    aveZmrks <- aveZmrks + Zmrks
    count <- count + 1
}

aveZmrks <- aveZmrks/count

x11(type="cairo",height=5,width=5*(dimns$xspan/dimns$yspan))

figure <- paste(plotpath,'-simno',min(simnos),'_',max(simnos),'-avesoln',sep='')
mcmc.solution.pt2(paste(mcmcpath,'-simno',max(simnos),sep=''),dimns,aveZmrks)
savePlot(filename=paste(figure,".png",sep=""),type="png")

dev.off( )
