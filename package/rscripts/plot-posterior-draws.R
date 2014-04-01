

rm(list=ls( ))

datapath <- '~/package/data'
plotpath <- '.'

source('src/myfilled.contour.R')
source('src/mydichromat.R')
library(fields) ## Computes Euclidean distances ##
library(deldir) ## Calculates Voronoi tessellation ##


mcmc.voronoi <- function(datapath,dataset,config,dimns) {
    
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

    xrange <- dimns$xrange
    yrange <- dimns$yrange
    
    niter <- length(ntiles)
    count <- 0
  
    for (i in 1:niter) {
    
        curr.ntiles <- ntiles[i]
        curr.mrates <- mrates[(count+1):(count+curr.ntiles)]
        curr.xcoord <- xcoord[(count+1):(count+curr.ntiles)]
        curr.ycoord <- ycoord[(count+1):(count+curr.ntiles)]

        curr.effcts <- curr.mrates - mean(curr.mrates)

        L <- length(mylevels)
        indices <- which(curr.effcts<mylevels[1])
        curr.effcts[indices] <- 0.999*mylevels[1]
        indices <- which(curr.effcts>mylevels[L])
        curr.effcts[indices] <- 0.999*mylevels[L]
        
        plot(0,0,xlim=xrange,ylim=yrange,asp=1,xlab='',ylab='',type='n')
        
        centers <- cbind(curr.xcoord,curr.ycoord)
        indlvls <- 1:(numlevels+1)

        if (curr.ntiles==1) {

            ## There is only one tile
            k <- max(indlvls[mylevels<curr.effcts])
            polygon(xrange,yrange,col=mycolors[k],border=FALSE)
            
        } else {
    
            D <- deldir(curr.xcoord,curr.ycoord,rw=c(1.1*xrange,1.1*yrange))
            tilelist <- tile.list(D)
            centroids <- tile.centroids(tilelist)
    
            ## centroids != centers (See deldir.pdf)
            for (c in 1:nrow(centroids)) {
                centroid <- matrix(tilelist[[c]]$pt,nrow=1,ncol=2)
                euDist <- rdist(centers,centroid)
                cc <- apply(euDist,2,which.min)
                k <- max(indlvls[mylevels<curr.effcts[cc]])
                polygon(tilelist[[c]]$x,tilelist[[c]]$y,col=mycolors[k],border=FALSE)
            }
        }

        ## Plot centers
        points(centers,pch=19)
        count <- count + curr.ntiles
    }
}


dataset <- 'hap-tribars-nIndiv300-nSites3000-gridSize12x8'
simno <- 1
xPop <- 12
yPop <- 8


dimns <- read.dimns(datapath,dataset,xPop,yPop)

png(height=5,width=5*(dimns$xspan/dimns$yspan),units="in",res=150,
    file=paste(plotpath,'/',dataset,'-g',xPop,'x',yPop,'-simno',simno,'-draws%02d.png',sep=''))

config <- paste('simno',simno,'-g',xPop,'x',yPop,sep='')
mcmc.voronoi(datapath,dataset,config,dimns)

dev.off( )
