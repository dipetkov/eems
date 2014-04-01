

rm(list=ls( ))

datapath <- '~/package/data'
plotpath <- '.'

source('src/mydichromat.R')
source('src/myfilled.contour.R')
library(fields)


mcmc.solution.pt1 <- function(datapath,dataset,config,dimns) {

    filepath <- paste(datapath,'/',dataset,'-',config,sep='')
    print(filepath)
    
    mrates <- scan(paste(filepath,'.mcmcmrates',sep=''))
    ntiles <- scan(paste(filepath,'.mcmcntiles',sep=''))
    xcoord <- scan(paste(filepath,'.mcmcxcoord',sep=''))
    ycoord <- scan(paste(filepath,'.mcmcycoord',sep=''))
    
    mrates <- log10(mrates)

    xcoord <- as.numeric(xcoord)
    ycoord <- as.numeric(ycoord)
    mrates <- as.numeric(mrates)
    
    xmrks <- dimns$xmrks
    ymrks <- dimns$ymrks
    nxmrks <- length(xmrks)
    nymrks <- length(ymrks)
    
    marks <- cbind(rep(xmrks,nymrks),rep(ymrks,each=nxmrks))
    aveZmrks <- matrix(0,nxmrks,nymrks)
    totmrks <- nxmrks*nymrks
    
    niter <- length(ntiles)
    count <- 0
    
    for (i in 1:niter) {
        curr.ntiles <- ntiles[i]
        curr.mrates <- mrates[(count+1):(count+curr.ntiles)]
        curr.xcoord <- xcoord[(count+1):(count+curr.ntiles)]
        curr.ycoord <- ycoord[(count+1):(count+curr.ntiles)]
        centers <- cbind(curr.xcoord,curr.ycoord)
        dist <- rdist(centers,marks)
        nuclei <- apply(dist,2,which.min)

        weights <- numeric(curr.ntiles)
        for (t in 1:curr.ntiles) { weights[t] <- sum(nuclei==t) }
        curr.effcts <- curr.mrates - weighted.mean(curr.mrates,weights)

        zmrks <- matrix(curr.effcts[nuclei],nxmrks,nymrks,byrow=FALSE)
        aveZmrks <- aveZmrks + zmrks
        count <- count + curr.ntiles
    }
    
    aveZmrks <- aveZmrks/niter
    return (aveZmrks)
}

mcmc.solution.pt2 <- function(datapath,dataset,config,dimns,aveZmrks) {
    
    filepath <- paste(datapath,'/',dataset,'-',config,sep='')
    print(filepath)
    
    coord <- read.table(paste(filepath,'.coord',sep=''),header=FALSE)
    edges <- read.table(paste(filepath,'.edges',sep=""),header=FALSE)
    
    nv <- nrow(edges)
    nn <- ncol(edges)
    nodes <- 1:nv
    
    xmrks <- dimns$xmrks
    ymrks <- dimns$ymrks
    nxmrks <- length(xmrks)
    nymrks <- length(ymrks)
    xrange <- dimns$xrange
    yrange <- dimns$yrange
    
    plot.new( )
    par(new = "TRUE",las = 1,cex.axis = 1, plt = c(0.05,0.84,0.05,0.95),tck = -0.02)
    
    myfilled.contour(xmrks,ymrks,aveZmrks,asp=1,
                     xlim=xrange,ylim=yrange,
                     col=mycolors,levels=mylevels,frame.plot=FALSE,
                     plot.axes = { xlab = ''; ylab = ''; 
                                   for (a in 1:nv) {
                                   for (i in 1:nn) {
                                       b <- edges[a,i]
                                       if (sum(nodes==b)) {
                                           lines(coord[c(a,b),1],
                                                 coord[c(a,b),2],col="gray80")
                                       }
                                   } }
                               })
        
    par(new = "TRUE", plt = c(0.87,0.91,0.1,0.9), las = 1, cex.axis = 1.2)
    
    myfilled.legend(xmrks,ymrks,aveZmrks,asp=1,
                    xlim=xrange,ylim=yrange,
                    colors=mycolors,levels=mylevels,
                    key.axes = axis(4,at=the.pts,labels=the.lab,
                        tick=FALSE,line=-0.5),
                    key.title = mtext(expression(paste("e"["c"],sep="")),
                        side=3,line=1))
}


dataset <- 'hap-tribars-nIndiv300-nSites3000-gridSize12x8'
xPop <- 12
yPop <- 8
simnos <- 1:6

nxmrks <- 50
nymrks <- 38
aveZmrks <- matrix(0,nxmrks,nymrks)
dimns <- read.dimns(datapath,dataset,xPop,yPop,nxmrks=nxmrks,nymrks=nymrks)
count <- 0

for (simno in simnos) {

    config <- paste('simno',simno,'-g',xPop,'x',yPop,sep='')    
    aveZmrks <- aveZmrks + mcmc.solution.pt1(datapath,dataset,config,dimns)
    count <- count + 1
    
}

x11(type="cairo",height=5,width=5*(dimns$xspan/dimns$yspan))
figure <- paste(plotpath,'/',dataset,'-g',xPop,'x',yPop,
        '-simno',min(simnos),'_',max(simnos),'-avesoln',sep='')
aveZmrks <- aveZmrks/count
mcmc.solution.pt2(datapath,dataset,config,dimns,aveZmrks)
savePlot(filename=paste(figure,".png",sep=""),type="png")

dev.off( )
