

source('src/myfilled.contour.R')
source('src/mydichromat.R')
library(fields) ## the rdist function computes Euclidean distances ##
library(deldir) ## the deldir function calculates Voronoi tessellation ##


round2 <- function(x) { trunc(x+0.5) }

read.dimns <- function(datapath,nxmrks=NULL,nymrks=NULL) {
    dimns <- read.table(paste(datapath,'.dimns',sep=''))
    xrange <- as.numeric(dimns[1,])
    yrange <- as.numeric(dimns[2,])
    xspan <- xrange[2]-xrange[1]
    yspan <- yrange[2]-yrange[1]
    if (is.null(nxmrks)&&is.null(nymrks)) {
        nxmrks <- 50
        nymrks <- round2(nxmrks*(yspan/xspan))
    }
    xmrks <- seq(xrange[1],xrange[2],length=nxmrks)
    ymrks <- seq(yrange[1],yrange[2],length=nymrks)
    return(list(xmrks=xmrks,xrange=xrange,xspan=xspan,
                ymrks=ymrks,yrange=yrange,yspan=yspan))
}    

mcmc.solution.pt1 <- function(mcmcpath,dimns) {

    print(mcmcpath)
    
    mrates <- scan(paste(mcmcpath,'.mcmcmrates',sep=''))
    ntiles <- scan(paste(mcmcpath,'.mcmcntiles',sep=''))
    xcoord <- scan(paste(mcmcpath,'.mcmcxcoord',sep=''))
    ycoord <- scan(paste(mcmcpath,'.mcmcycoord',sep=''))
    
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

mcmc.solution.pt2 <- function(mcmcpath,dimns,aveZmrks) {
    
    print(mcmcpath)
    
    ipmap <- read.table(paste(mcmcpath,'.ipmap',sep=''),header=FALSE)
    coord <- read.table(paste(mcmcpath,'.coord',sep=''),header=FALSE)
    edges <- read.table(paste(mcmcpath,'.edges',sep=""),header=FALSE)
    ipmap <- as.numeric(ipmap)
    
    nv <- nrow(edges)
    nn <- ncol(edges)
    nodes <- 1:nv
    
    xmrks <- dimns$xmrks
    ymrks <- dimns$ymrks
    nxmrks <- length(xmrks)
    nymrks <- length(ymrks)
    xrange <- dimns$xrange
    yrange <- dimns$yrange
    
    sizes <- table(ipmap)
    demes <- as.numeric(names(sizes))
    sizes <- as.numeric(sizes)

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
                                   points(coord[demes,1],pch=19,
                                          coord[demes,2],col="gray10",
                                          cex=1+sizes/max(sizes));
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

mcmc.solution <- function(mcmcpath,dimns) {
    aveZmrks <- mcmc.solution.pt1(mcmcpath,dimns)
    mcmc.solution.pt2(mcmcpath,dimns,aveZmrks)
}

mcmc.voronoi <- function(mcmcpath,dimns) {
    
    print(mcmcpath)
  
    mrates <- scan(paste(mcmcpath,'.mcmcmrates',sep=''))
    ntiles <- scan(paste(mcmcpath,'.mcmcntiles',sep=''))
    xcoord <- scan(paste(mcmcpath,'.mcmcxcoord',sep=''))
    ycoord <- scan(paste(mcmcpath,'.mcmcycoord',sep=''))
    
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
    
            ## centroids != centers
            for (c in 1:nrow(centroids)) {
                centroid <- matrix(tilelist[[c]]$pt,nrow=1,ncol=2)
                euDist <- rdist(centers,centroid)
                cc <- apply(euDist,2,which.min)
                k <- max(indlvls[mylevels<curr.effcts[cc]])
                polygon(tilelist[[c]]$x,tilelist[[c]]$y,col=mycolors[k],border=FALSE)
            }
        }

        ## Plot the centers
        points(centers,pch=19)
        count <- count + curr.ntiles
    }
}

plot.log.posterior <- function(mcmcpath,simnos) {
    
    print(mcmcpath)
    
    pilogl <- read.table(paste(mcmcpath,'-simno',simnos[1],'.mcmcpilogl',sep=''),header=FALSE)
    nsimno <- length(simnos)
    niters <- nrow(pilogl)
    posterior <- matrix(0,niters,nsimno)
    
    for (i in 1:nsimno) {
        pilogl <- read.table(paste(mcmcpath,'-simno',simnos[i],'.mcmcpilogl',sep=''),header=FALSE)
        posterior[,i] <- pilogl[,1] + pilogl[,2]
    }
    
    matplot(posterior,type="l",col=varcolors,
            xlab="iteration",ylab="log posterior",main=mcmcpath,
            lwd=3,lty=1,axes=FALSE,xlim=c(1,1.06*niters),bty="n")
    legend(1.01*niters,max(posterior),legend=simnos,col=varcolors,lwd=5,bty="n")
}
