

source('src/myfilled.contour.R')
source('src/mydichromat.R')
library(fields) ## the rdist function computes Euclidean distances ##
library(deldir) ## the deldir function calculates Voronoi tessellation ##


round2 <- function(x) { trunc(x+0.5) }

read.dimns <- function(datapath,nxmrks=NULL,nymrks=NULL) {
    dimns <- scan(paste(datapath,'.dimns',sep=''),what=numeric())
    xrange <- as.numeric(dimns[1:2])
    yrange <- as.numeric(dimns[3:4])
    xspan <- xrange[2]-xrange[1]
    yspan <- yrange[2]-yrange[1]
    if (is.null(nxmrks)&&is.null(nymrks)) {
        nxmrks <- 50
        nymrks <- round2(nxmrks*(yspan/xspan))
    }
    xmrks <- seq(xrange[1],xrange[2],length=nxmrks)
    ymrks <- seq(yrange[1],yrange[2],length=nymrks)
    return(list(nxmrks=nxmrks,xmrks=xmrks,xrange=xrange,xspan=xspan,
                nymrks=nymrks,ymrks=ymrks,yrange=yrange,yspan=yspan))
}    

mcmc.mrates.pt1 <- function(mcmcpath,dimns) {
    print('Computing posterior mean surface of mrates')
    print(mcmcpath)
    
    mrates <- scan(paste(mcmcpath,'.mcmcmrates',sep=''),what=numeric())
    mtiles <- scan(paste(mcmcpath,'.mcmcmtiles',sep=''),what=numeric())
    xcoord <- scan(paste(mcmcpath,'.mcmcxcoord',sep=''),what=numeric())
    ycoord <- scan(paste(mcmcpath,'.mcmcycoord',sep=''),what=numeric())
    
    mrates <- log10(mrates)

    xmrks <- dimns$xmrks
    ymrks <- dimns$ymrks
    nxmrks <- length(xmrks)
    nymrks <- length(ymrks)
    
    marks <- cbind(rep(xmrks,nymrks),rep(ymrks,each=nxmrks))
    Zmrks <- matrix(0,nxmrks,nymrks)
    totmrks <- nxmrks*nymrks
    
    niter <- length(mtiles)
    count <- 0
    
    for (i in 1:niter) {
        curr.mtiles <- mtiles[i]
        curr.mrates <- mrates[(count+1):(count+curr.mtiles)]
        curr.xcoord <- xcoord[(count+1):(count+curr.mtiles)]
        curr.ycoord <- ycoord[(count+1):(count+curr.mtiles)]
        centers <- cbind(curr.xcoord,curr.ycoord)
        distances <- rdist(centers,marks)
        closest <- apply(distances,2,which.min)
        ## Standardize the log-transformed migration rates
        ## by weighting each tile by its relative size
        ## Rather than compute the exact area of each polygon,
        ## approximate its area by the number of zmrks that
        ## fall within the tile  
        weights <- numeric(curr.mtiles)
        for (t in 1:curr.mtiles) { weights[t] <- sum(closest==t) }
        curr.effcts <- curr.mrates - weighted.mean(curr.mrates,weights)
        zmrks <- matrix(curr.effcts[closest],nxmrks,nymrks,byrow=FALSE)
        Zmrks <- Zmrks + zmrks/niter
        count <- count + curr.mtiles
    }
    return (Zmrks)
}

mcmc.mrates.pt2 <- function(mcmcpath,dimns,Zmrks) {
    print('Plotting posterior mean surface of mrates')
    print(mcmcpath)
    
    ipmap <- scan(paste(mcmcpath,'.ipmap',sep=''),what=numeric())
    coord <- scan(paste(mcmcpath,'.demes',sep=''),what=numeric())
    edges <- scan(paste(mcmcpath,'.edges',sep=''),what=numeric())
    coord <- matrix(coord,ncol=2,byrow=TRUE)
    edges <- matrix(edges,ncol=6,byrow=TRUE)
    
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
    par(new = "TRUE",las = 1,cex.axis = 1,plt = c(0.05,0.84,0.05,0.95),tck = -0.02,xpd = TRUE)
    
    myfilled.contour(xmrks,ymrks,Zmrks,asp=1,
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
        
    par(new = "TRUE",las = 1,cex.axis = 1.2,plt = c(0.87,0.91,0.1,0.9))
    
    myfilled.legend(xmrks,ymrks,Zmrks,asp=1,
                    xlim=xrange,ylim=yrange,
                    colors=mycolors,levels=mylevels,
                    key.axes = axis(4,at=the.pts,labels=the.lab,
                        tick=FALSE,line=-0.5),
                    key.title = mtext(expression(paste("e"["m"],sep="")),
                        side=3,line=1))
}

mcmc.mrates <- function(mcmcpath,dimns) {
    Zmrks <- mcmc.mrates.pt1(mcmcpath,dimns)
    mcmc.mrates.pt2(mcmcpath,dimns,Zmrks)
}

mcmc.mrates.simnos <- function(mcmcpath,simnos,dimns) {
    Zmrks <- matrix(0,dimns$nxmrks,dimns$nymrks)
    for (simno in simnos) {
        Zmrks <- Zmrks + mcmc.mrates.pt1(paste(mcmcpath,'-simno',simno,sep=''),dimns)
    }
    Zmrks <- Zmrks/length(simnos)
    mcmc.mrates.pt2(paste(mcmcpath,'-simno',simnos[1],sep=''),dimns,Zmrks)
}

mcmc.voronoi <- function(mcmcpath,dimns) {
    print('Plotting Voronoi tessellation')    
    print(mcmcpath)
  
    mrates <- scan(paste(mcmcpath,'.mcmcmrates',sep=''),what=numeric())
    mtiles <- scan(paste(mcmcpath,'.mcmcmtiles',sep=''),what=numeric())
    xcoord <- scan(paste(mcmcpath,'.mcmcxcoord',sep=''),what=numeric())
    ycoord <- scan(paste(mcmcpath,'.mcmcycoord',sep=''),what=numeric())
    
    mrates <- log10(mrates)

    xrange <- dimns$xrange
    yrange <- dimns$yrange
    
    niter <- length(mtiles)
    count <- 0
  
    for (i in 1:niter) {
    
        curr.mtiles <- mtiles[i]
        curr.mrates <- mrates[(count+1):(count+curr.mtiles)]
        curr.xcoord <- xcoord[(count+1):(count+curr.mtiles)]
        curr.ycoord <- ycoord[(count+1):(count+curr.mtiles)]
        ## Standardize the log-transformed migration rates
        ## without taking into account the relative size of the tiles
        ## (only because it is harder to do without the Zmrks grid) 
        curr.effcts <- curr.mrates - mean(curr.mrates)

        L <- length(mylevels)
        indices <- which(curr.effcts<mylevels[1])
        curr.effcts[indices] <- 0.999*mylevels[1]
        indices <- which(curr.effcts>mylevels[L])
        curr.effcts[indices] <- 0.999*mylevels[L]
        centers <- cbind(curr.xcoord,curr.ycoord)
        
        plot.new( )
        par(new = TRUE,plt = c(0,1,0,1),las = 1,cex.axis = 1)
        plot.window(xlim=xrange,ylim=yrange,asp=1,xlab='',ylab='')

        if (curr.mtiles==1) {

            ## There is only one tile
            which.tile <- max((1:L)[mylevels<curr.effcts])
            polygon(xrange,yrange,col=mycolors[which.tile])
            
        } else {
            
            Voronoi <- deldir(curr.xcoord,curr.ycoord,rw=c(xrange,yrange))
            tilelist <- tile.list(Voronoi)
            centroids <- tile.centroids(tilelist)
            
            ## centroids != centers
            for (c in 1:nrow(centroids)) {

                centroid <- matrix(tilelist[[c]]$pt,nrow=1,ncol=2)
                euDist <- rdist(centers,centroid)
                closest <- apply(euDist,2,which.min)
                which.tile <- max((1:L)[mylevels<curr.effcts[closest]])
                polygon(tilelist[[c]]$x,tilelist[[c]]$y,col=mycolors[which.tile])
            }
        }

        ## Plot the centers
        points(centers,pch=19)
        count <- count + curr.mtiles
    }
}

dist.scatterplot <- function(mcmcpath) {
    print('Plotting average distances between demes')
    print(mcmcpath)

    JtDobsJ <- read.table(paste(mcmcpath,'.rdistJtDobsJ',sep=''),header=FALSE)
    JtDhatJ <- read.table(paste(mcmcpath,'.rdistJtDhatJ',sep=''),header=FALSE)
    JtDobsJ <- as.matrix(JtDobsJ)
    JtDhatJ <- as.matrix(JtDhatJ)

    nPop <- nrow(JtDobsJ)    
    Bobs <- JtDobsJ
    Bhat <- JtDhatJ
    nPop <- nrow(JtDobsJ)
    Wobs <- matrix(diag(Bobs),nPop,nPop)
    What <- matrix(diag(Bhat),nPop,nPop)
    Wobs <- (Wobs+t(Wobs))/2
    What <- (What+t(What))/2

    Dobs.ab <- Bobs-Wobs
    Dist.ab <- Bhat-What
    Dobs.ab <- Dobs.ab[upper.tri(Dobs.ab,diag=TRUE)]
    Dist.ab <- Dist.ab[upper.tri(Dist.ab,diag=TRUE)]

    plot(Dist.ab,Dobs.ab,type="n",
         xlab=expression(paste("Fitted between distances  ",hat(D)[ab],"-(",hat(D)[aa],"+",hat(D)[bb],")/2",sep="")),
         ylab=expression(paste("Observed between distances  ",D[ab],"-(",D[aa],"+",D[bb],")/2",sep="")))
    abline(a=0,b=1,col="red",lwd=2)
    points(Dist.ab,Dobs.ab)
    print(summary(lm(Dobs.ab~Dist.ab))$r.squared)
}
