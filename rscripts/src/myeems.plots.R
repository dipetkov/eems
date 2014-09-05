
###########################################################################################
###########################################################################################
## This is the default Dark Orange to Blue color scheme, with "white"
## as the midpoint color.
## It combines two color schemes from the dichromat package, which
## itself is based on a collection of color schemes for scientific
## data graphics:
## http://geog.uoregon.edu/datagraphics/color_scales.htm
## and 
## Light A and Bartlein PJ (2004). The End of the Rainbow? Color Schemes
## for Improved Data Graphics.
## EOS Transactions of the American Geophysical Union, 85(40), 385.

## The default Dark Orange to Blue color scheme:
eems.colors <- c("#994000","#CC5800","#FF8F33","#FFAD66","#FFCA99","#FFE6CC",
                 "#FFFFFF",
                 "#CCFDFF","#99F8FF","#66F0FF","#33E4FF","#00AACC","#007A99")
numlevels <- length(eems.colors)

###########################################################################################
###########################################################################################
## modification by Ian Taylor of the filled.contour function
## to remove the key and facilitate overplotting with contour( )
## further modified by Carey McGilliard and Bridget Ferris
## to allow multiple plots on one page

## Downloaded filled.contour3.R and filled.legend.R from:
## http://wiki.cbr.washington.edu/qerm/index.php/R/Contour_Plots

myfilled.contour <-
    function (x = seq(0, 1, length.out = nrow(z)),
              y = seq(0, 1, length.out = ncol(z)), z, xlim = range(x, finite = TRUE), 
              ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE), 
              levels = pretty(zlim, nlevels), nlevels = 20, color.palette = cm.colors, 
              col = color.palette(length(levels) - 1), plot.title, plot.axes, 
              key.title, key.axes, asp = NA, xaxs = "i", yaxs = "i", las = 1, 
              axes = TRUE, frame.plot = axes,mar, ...) 
{
    
    if (missing(z)) {
        if (!missing(x)) {
            if (is.list(x)) {
                z <- x$z
                y <- x$y
                x <- x$x
            } else {
                z <- x
                x <- seq.int(0, 1, length.out = nrow(z))
            }
        } else stop("no 'z' matrix specified")
    } else if (is.list(x)) {
        y <- x$y
        x <- x$x
    }
    if (any(diff(x) <= 0) || any(diff(y) <= 0)) {
        stop("increasing 'x' and 'y' values expected")
    }
    plot.new( )

    plot.window(xlim=xlim, ylim=ylim, xaxs = xaxs, yaxs = yaxs, asp = asp)

    if (!is.matrix(z) || nrow(z) <= 1 || ncol(z) <= 1) 
        stop("no proper 'z' matrix specified")
    if (!is.double(z)) 
        storage.mode(z) <- "double"
    .filled.contour(as.double(x), as.double(y), z, as.double(levels),
                    col = col)
    if (frame.plot) { box( ) }
    if (missing(plot.title)) { title(...) }
    else { plot.title }
    if (!missing(plot.axes)) { plot.axes }
    invisible( )
}
## modification of filled.contour by Carey McGilliard and Bridget Ferris
## designed to just plot the legend
myfilled.legend <-
    function (x = seq(0, 1, length.out = nrow(z)),
              y = seq(0, 1, length.out = ncol(z)), z, xlim = range(x, finite = TRUE), 
              ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE), 
              levels = pretty(zlim, nlevels), nlevels = 20 , color.palette = cm.colors, 
              colors = color.palette(length(levels) - 1), plot.title, plot.axes, 
              key.title, key.axes, asp = NA, xaxs = "i", yaxs = "i", las = 1, 
              axes = TRUE, frame.plot = axes, ...) 
{

    if (missing(z)) {
        if (!missing(x)) {
            if (is.list(x)) {
                z <- x$z
                y <- x$y
                x <- x$x
            } else {
                z <- x
                x <- seq.int(0, 1, length.out = nrow(z))
            }
        }
        else stop("no 'z' matrix specified")
    } else if (is.list(x)) {
        y <- x$y
        x <- x$x
    }
    if (any(diff(x) <= 0) || any(diff(y) <= 0)) {
        stop("increasing 'x' and 'y' values expected")
    }
    plot.new( )
    nlevels = length(levels)

    plot.window(xlim = c(0,1), ylim = range(levels), xaxs = "i", yaxs = "i")
    
    rect(0, levels[-nlevels], 1, levels[-1], col = colors, border="white")
    if (!missing(key.axes)) {
        key.axes
    } else {
        axis(4,tick = FALSE)
    }
    if (!missing(key.title)) {
        key.title
    }
}

###########################################################################################
###########################################################################################
library(fields) ## the rdist function computes Euclidean distances ##
library(deldir) ## the deldir function calculates Voronoi tessellation ##

round2 <- function(x) { trunc(x+0.5) }
center.mrates <- function(Zmrks) {
    minZ <- min(Zmrks)
    maxZ <- max(Zmrks)
    numlevels <- length(eems.colors)
    if ((minZ > -2.5)&&(maxZ < +2.5)) {
        minZ <- -2.5
        maxZ <- +2.5
        step <- (maxZ-minZ)/numlevels
        eems.levels <- seq(from=minZ,to=maxZ,by=step)
    } else {
        ## The mrates should be centered at about zero
        maxZ2 <- max(maxZ,-minZ)
        maxZ2 <- ceiling(maxZ2*1000)/1000
        minZ <- -maxZ2
        maxZ <- maxZ2
        eems.levels <- seq(from=minZ,to=maxZ,length=numlevels)
    }
    return(eems.levels)
}
center.qrates <- function(Zmrks) {
    minZ <- min(Zmrks)
    maxZ <- max(Zmrks)
    numlevels <- length(eems.colors)
    ## The qrates should be centered at about zero
    maxZ2 <- max(maxZ,-minZ)
    maxZ2 <- ceiling(maxZ2*1000)/1000
    minZ <- -maxZ2
    maxZ <- maxZ2
    if ((maxZ - minZ)>0.001) {
        eems.levels <- seq(from=minZ,to=maxZ,length=numlevels)
    } else {
        eems.levels <- seq(from=minZ-0.001/2,to=maxZ+0.001/2,length=numlevels)
    }
    return(eems.levels)
}
read.dimns <- function(datapath,nxmrks=NULL,nymrks=NULL) {
    dimns <- scan(paste(datapath,'.dimns',sep=''),what=numeric())
    xrange <- as.numeric(dimns[1:2])
    yrange <- as.numeric(dimns[3:4])
    xspan <- xrange[2]-xrange[1]
    yspan <- yrange[2]-yrange[1]
    if (is.null(nxmrks)&&is.null(nymrks)) {
        nxmrks <- 100
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
mcmc.mrates.pt2 <- function(mcmcpath,dimns,Zmrks,add.grid=TRUE,add.samples=TRUE) {
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
    sizes <- table(ipmap)
    demes <- as.numeric(names(sizes))
    sizes <- as.numeric(sizes)
    eems.levels <- center.mrates(Zmrks)
    plot.new( )
    par(new = TRUE, plt = c(0,1,0,1), las = 1, cex.axis = 1)
    myfilled.contour(dimns$xmrks,dimns$ymrks,Zmrks,asp=1,
                     xlim=dimns$xrange,ylim=dimns$yrange,
                     col=eems.colors,levels=eems.levels,frame.plot=FALSE,
                     plot.axes = { xlab = ''; ylab = '';
                                   if (add.grid) {
                                       for (a in 1:nv) {
                                       for (i in 1:nn) {
                                           b <- edges[a,i]
                                           if (b %in% nodes) {
                                               lines(coord[c(a,b),1],
                                                     coord[c(a,b),2],col="gray80")
                                           }
                                       } }
                                   }
                                   if (add.samples) {
                                       points(coord[demes,1],pch=19,
                                              coord[demes,2],col="gray10",
                                              cex=1+sizes/max(sizes));
                                   }
                               })
    return(list(colors=eems.colors,levels=eems.levels))
}
mcmc.mrates <- function(mcmcpath,dimns,add.grid=TRUE,add.samples=TRUE) {
    Zmrks <- mcmc.mrates.pt1(mcmcpath,dimns)
    eems.legend <- NULL
    if (file.exists(paste(mcmcpath,'.mcmcmrates',sep=''))) {
        eems.legend <- mcmc.mrates.pt2(mcmcpath,dimns,Zmrks,add.grid,add.samples)
    }
    return(eems.legend)
}
mcmc.mrates.simnos <- function(mcmcpath,simnos,dimns,add.grid=TRUE,add.samples=TRUE) {
    Zmrks <- matrix(0,dimns$nxmrks,dimns$nymrks)
    nSimno <- 0
    existsSimno <- 0
    for (simno in simnos) {
        File <- paste(mcmcpath,'-simno',simno,sep='')
        if (file.exists(paste(File,'.mcmcmrates',sep=''))) {
            Zmrks <- Zmrks + mcmc.mrates.pt1(File,dimns)
            nSimno <- nSimno + 1
            existsSimno <- simno
        }
    }
    eems.legend <- NULL
    if (existsSimno) {
        Zmrks <- Zmrks/nSimno
        eems.legend <- mcmc.mrates.pt2(paste(mcmcpath,'-simno',existsSimno,sep=''),dimns,Zmrks,add.grid,add.samples)
    }
    return(eems.legend)
}
mcmc.qrates.pt1 <- function(mcmcpath,dimns) {
    print('Computing posterior mean surface of qrates')
    print(mcmcpath)
    mrates <- scan(paste(mcmcpath,'.mcmcqrates',sep=''),what=numeric())
    qtiles <- scan(paste(mcmcpath,'.mcmcqtiles',sep=''),what=numeric())
    xcoord <- scan(paste(mcmcpath,'.mcmcwcoord',sep=''),what=numeric())
    ycoord <- scan(paste(mcmcpath,'.mcmczcoord',sep=''),what=numeric())
    mrates <- log10(mrates)
    xmrks <- dimns$xmrks
    ymrks <- dimns$ymrks
    nxmrks <- length(xmrks)
    nymrks <- length(ymrks)
    marks <- cbind(rep(xmrks,nymrks),rep(ymrks,each=nxmrks))
    Zmrks <- matrix(0,nxmrks,nymrks)
    niter <- length(qtiles)
    count <- 0
    for (i in 1:niter) {
        curr.qtiles <- qtiles[i]
        curr.mrates <- mrates[(count+1):(count+curr.qtiles)]
        curr.xcoord <- xcoord[(count+1):(count+curr.qtiles)]
        curr.ycoord <- ycoord[(count+1):(count+curr.qtiles)]
        centers <- cbind(curr.xcoord,curr.ycoord)
        distances <- rdist(centers,marks)
        closest <- apply(distances,2,which.min)
        weights <- numeric(curr.qtiles)
        for (t in 1:curr.qtiles) { weights[t] <- sum(closest==t) }
        curr.effcts <- curr.mrates - weighted.mean(curr.mrates,weights)
        zmrks <- matrix(curr.effcts[closest],nxmrks,nymrks,byrow=FALSE)
        Zmrks <- Zmrks + zmrks/niter
        count <- count + curr.qtiles
    }
    return (Zmrks)
}
mcmc.qrates.pt2 <- function(mcmcpath,dimns,Zmrks,add.grid=TRUE,add.samples=TRUE) {
    print('Plotting posterior mean surface of qrates')
    print(mcmcpath)
    ipmap <- scan(paste(mcmcpath,'.ipmap',sep=''),what=numeric())
    coord <- scan(paste(mcmcpath,'.demes',sep=''),what=numeric())
    edges <- scan(paste(mcmcpath,'.edges',sep=''),what=numeric())
    coord <- matrix(coord,ncol=2,byrow=TRUE)
    edges <- matrix(edges,ncol=6,byrow=TRUE)
    nv <- nrow(edges)
    nn <- ncol(edges)
    nodes <- 1:nv
    sizes <- table(ipmap)
    demes <- as.numeric(names(sizes))
    sizes <- as.numeric(sizes)
    eems.levels <- center.qrates(Zmrks)
    plot.new( )
    eems.levels <- center.qrates(Zmrks)
    par(new = TRUE, plt = c(0,1,0,1), las = 1, cex.axis = 1)
    myfilled.contour(dimns$xmrks,dimns$ymrks,Zmrks,asp=1,
                     xlim=dimns$xrange,ylim=dimns$yrange,
                     col=eems.colors,levels=eems.levels,frame.plot=FALSE,
                     plot.axes = { xlab = ''; ylab = '';
                                   if (add.grid) {
                                       for (a in 1:nv) {
                                       for (i in 1:nn) {
                                           b <- edges[a,i]
                                           if (b %in% nodes) {
                                               lines(coord[c(a,b),1],
                                                     coord[c(a,b),2],col="gray80")
                                           }
                                       } }
                                   }
                                   if (add.samples) {
                                       points(coord[demes,1],pch=19,
                                              coord[demes,2],col="gray10",
                                              cex=1+sizes/max(sizes))
                                   }
                               })
    return(list(colors=eems.colors,levels=eems.levels))
}
mcmc.qrates <- function(mcmcpath,dimns,add.grid=TRUE,add.samples=TRUE) {
    Zmrks <- mcmc.qrates.pt1(mcmcpath,dimns)
    eems.legend <- NULL
    if (file.exists(paste(mcmcpath,'.mcmcmrates',sep=''))) {
        eems.legend <- mcmc.qrates.pt2(mcmcpath,dimns,Zmrks,add.grid,add.samples)
    }
    return(eems.legend)
}
mcmc.qrates.simnos <- function(mcmcpath,simnos,dimns,add.grid=TRUE,add.samples=TRUE) {
    Zmrks <- matrix(0,dimns$nxmrks,dimns$nymrks)
    nSimno <- 0
    existsSimno <- 0
    for (simno in simnos) {
        File <- paste(mcmcpath,'-simno',simno,sep='')
        if (file.exists(paste(File,'.mcmcqrates',sep=''))) {
            Zmrks <- Zmrks + mcmc.qrates.pt1(File,dimns)
            nSimno <- nSimno + 1
            existsSimno <- simno
        }
    }
    eems.legend <- NULL
    if (existsSimno) {
        Zmrks <- Zmrks/nSimno
        eems.legend <- mcmc.qrates.pt2(paste(mcmcpath,'-simno',simnos[1],sep=''),dimns,Zmrks)
    }
    return(eems.legend)
}
mcmc.mrates.legend <- function(datapath,mcmcpath,dimns,Zmrks,legend) {
    plot.new( )
    par(new = TRUE, plt = c(0,0.3,0,1), las = 1, cex.axis = 1.2)
    myfilled.legend(dimns$xmrks,dimns$ymrks,Zmrks,asp=1,
                    xlim=dimns$xrange,ylim=dimns$yrange,
                    colors=legend$colors,levels=legend$levels,
                    key.axes = axis(4,tick=FALSE,line=-0.5),
                    key.title = mtext(expression(paste("e"["m"],sep="")),side=3,line=1))
}
mcmc.qrates.legend <- function(datapath,mcmcpath,dimns,Zmrks,legend) {
    plot.new( )
    par(new = TRUE, plt = c(0,0.3,0,1), las = 1, cex.axis = 1.2)
    myfilled.legend(dimns$xmrks,dimns$ymrks,Zmrks,asp=1,
                    xlim=dimns$xrange,ylim=dimns$yrange,
                    colors=legend$colors,levels=legend$levels,
                    key.axes = axis(4,tick=FALSE,line=-0.5),
                    key.title = mtext(expression(paste("e"["q"],sep="")),side=3,line=1))
}
mcmc.mrates.voronoi <- function(mcmcpath,dimns) {
    print('Plotting mrates Voronoi tessellation')
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
    eems.levels <- center.mrates(mrates)
    for (i in 1:niter) {
        curr.mtiles <- mtiles[i]
        curr.mrates <- mrates[(count+1):(count+curr.mtiles)]
        curr.xcoord <- xcoord[(count+1):(count+curr.mtiles)]
        curr.ycoord <- ycoord[(count+1):(count+curr.mtiles)]
        ## Standardize the log-transformed migration rates
        ## without taking into account the relative size of the tiles
        ## (only because it is harder to do without the Zmrks grid) 
        curr.effcts <- curr.mrates - mean(curr.mrates)
        L <- length(eems.levels)
        indices <- which(curr.effcts<eems.levels[1])
        curr.effcts[indices] <- 0.999*eems.levels[1]
        indices <- which(curr.effcts>eems.levels[L])
        curr.effcts[indices] <- 0.999*eems.levels[L]
        centers <- cbind(curr.xcoord,curr.ycoord)
        plot.new( )
        par(new = TRUE,plt = c(0,1,0,1),las = 1,cex.axis = 1)
        plot.window(xlim=xrange,ylim=yrange,asp=1,xlab='',ylab='')
        if (curr.mtiles==1) {
            ## There is only one tile
            which.tile <- max((1:L)[eems.levels<curr.effcts])
            polygon(xrange,yrange,col=eems.colors[which.tile],border=FALSE)
        } else {
            Voronoi <- deldir(curr.xcoord,curr.ycoord,rw=c(xrange,yrange))
            tilelist <- tile.list(Voronoi)
            centroids <- tile.centroids(tilelist)
            ## centroids != centers
            for (c in 1:nrow(centroids)) {
                centroid <- matrix(tilelist[[c]]$pt,nrow=1,ncol=2)
                euDist <- rdist(centers,centroid)
                closest <- apply(euDist,2,which.min)
                which.tile <- max((1:L)[eems.levels<curr.effcts[closest]])
                polygon(tilelist[[c]]$x,tilelist[[c]]$y,col=eems.colors[which.tile],border=FALSE)
            }
        }
        ## Plot the centers
        points(centers,pch=19)
        count <- count + curr.mtiles
    }
    return(list(colors=eems.colors,levels=eems.levels))    
}
mcmc.qrates.voronoi <- function(mcmcpath,dimns) {
    print('Plotting qrates Voronoi tessellation')
    print(mcmcpath)
    qrates <- scan(paste(mcmcpath,'.mcmcqrates',sep=''),what=numeric())
    qtiles <- scan(paste(mcmcpath,'.mcmcqtiles',sep=''),what=numeric())
    xcoord <- scan(paste(mcmcpath,'.mcmcwcoord',sep=''),what=numeric())
    ycoord <- scan(paste(mcmcpath,'.mcmczcoord',sep=''),what=numeric())
    qrates <- log10(qrates)
    xrange <- dimns$xrange
    yrange <- dimns$yrange
    niter <- length(qtiles)
    count <- 0
    eems.levels <- center.qrates(qrates)
    for (i in 1:niter) {
        curr.qtiles <- qtiles[i]
        curr.qrates <- qrates[(count+1):(count+curr.qtiles)]
        curr.xcoord <- xcoord[(count+1):(count+curr.qtiles)]
        curr.ycoord <- ycoord[(count+1):(count+curr.qtiles)]
        ## Standardize the log-transformed migration rates
        ## without taking into account the relative size of the tiles
        ## (only because it is harder to do without the Zmrks grid) 
        curr.effcts <- curr.qrates - mean(curr.qrates)
        L <- length(eems.levels)
        indices <- which(curr.effcts<eems.levels[1])
        curr.effcts[indices] <- 0.999*eems.levels[1]
        indices <- which(curr.effcts>eems.levels[L])
        curr.effcts[indices] <- 0.999*eems.levels[L]
        centers <- cbind(curr.xcoord,curr.ycoord)
        plot.new( )
        par(new = TRUE,plt = c(0,1,0,1),las = 1,cex.axis = 1)
        plot.window(xlim=xrange,ylim=yrange,asp=1,xlab='',ylab='')
        if (curr.qtiles==1) {
            ## There is only one tile
            which.tile <- max((1:L)[eems.levels<curr.effcts])
            polygon(xrange,yrange,col=eems.colors[which.tile],border=FALSE)
        } else {
            Voronoi <- deldir(curr.xcoord,curr.ycoord,rw=c(xrange,yrange))
            tilelist <- tile.list(Voronoi)
            centroids <- tile.centroids(tilelist)
            ## centroids != centers
            for (c in 1:nrow(centroids)) {
                centroid <- matrix(tilelist[[c]]$pt,nrow=1,ncol=2)
                euDist <- rdist(centers,centroid)
                closest <- apply(euDist,2,which.min)
                which.tile <- max((1:L)[eems.levels<curr.effcts[closest]])
                polygon(tilelist[[c]]$x,tilelist[[c]]$y,col=eems.colors[which.tile],border=FALSE)
            }
        }
        ## Plot the centers
        points(centers,pch=19)
        count <- count + curr.qtiles
    }
    return(list(colors=eems.colors,levels=eems.levels))    
}
plot.logposterior.simnos <- function(mcmcpath,simnos) {
    print('Plotting (log) posterior')
    print(mcmcpath)
    nSimno <- 0
    existsSimno <- 0
    for (simno in simnos) {
        File <- paste(mcmcpath,'-simno',simno,'.mcmcpilogl',sep='')
        if (file.exists(File)) {
            nSimno <- nSimno + 1
            existsSimno <- simno
        }
    }
    if (existsSimno) {
        pilogl <- read.table(paste(mcmcpath,'-simno',existsSimno,'.mcmcpilogl',sep=''),colClasses=numeric())
        pilogl <- as.matrix(pilogl)
        nIters <- nrow(pilogl)
        posteriors <- matrix(0,nIters,0)
        plotSimnos <- numeric( )
        for (simno in simnos) {
            File <- paste(mcmcpath,'-simno',simno,'.mcmcpilogl',sep='')
            if (file.exists(File)) {
                pilogl <- read.table(File,colClasses=numeric())
                pilogl <- as.matrix(pilogl)
                if (nrow(pilogl)!=nIters) {
                    print('Warning: nrow(pilogl)!=nIters: ',File,sep='')
                } else {
                    posteriors <- cbind(posteriors,apply(pilogl,1,sum))
                    plotSimnos <- c(plotSimnos,simno)
                }
            }
        }
        par(xpd=TRUE)
        x = strsplit(mcmcpath,'/')[[1]]
        filename = x[length(x)]
        matplot(posteriors,type="l",col=1:8,
                xlab="iteration (after thinning)",ylab="log posterior",main=filename,
                lwd=3,lty=1,xlim=c(1,1.06*nIters),bty="n")
        legend(1.01*nIters,max(posteriors),legend=plotSimnos,col=1:8,lwd=5,bty="n")
    }
}
dist.scatterplot <- function(mcmcpath,singletons=TRUE) {
    print('Plotting average distances between demes')
    print(mcmcpath)
    Sizes <- scan(paste(mcmcpath,'.rdistSizes',sep=''))
    nPops <- length(Sizes)
    Size1 <- matrix(Sizes,nPops,nPops)
    Size2 <- t(Size1)
    JtDobsJ <- read.table(paste(mcmcpath,'.rdistJtDobsJ',sep=''),header=FALSE)
    JtDhatJ <- read.table(paste(mcmcpath,'.rdistJtDhatJ',sep=''),header=FALSE)
    JtDobsJ <- as.matrix(JtDobsJ)
    JtDhatJ <- as.matrix(JtDhatJ)
    minSize <- pmin(Size1,Size2)
    ypts <- JtDobsJ[upper.tri(JtDobsJ,diag=TRUE)]
    xpts <- JtDhatJ[upper.tri(JtDhatJ,diag=TRUE)]
    par(font.main=1)
    plot(xpts,ypts,type="n",
         ylab=expression(paste("Observed distances between and within demes  ",D[ab],sep="")),
         xlab=expression(paste("Fitted distances between and within demes  ",hat(D)[ab],sep="")),
         main="Distances between and within demes")
    abline(a=0,b=1,col="red",lwd=2)
    points(xpts,ypts)
    Wobs <- diag(JtDobsJ)
    What <- diag(JtDhatJ)
    ones <- matrix(1,nPops,1)
    Bobs <- JtDobsJ - (Wobs%*%t(ones) + ones%*%t(Wobs))/2
    Bhat <- JtDhatJ - (What%*%t(ones) + ones%*%t(What))/2
    if (!singletons) {
        ## Remove singleton locations
        remove <- (1:nPops)[Sizes<=1]
        if (length(remove)) {
            Bobs <- Bobs[-remove,-remove]
            Bhat <- Bhat[-remove,-remove]
            Wobs <- Wobs[-remove]
            What <- What[-remove]
            Sizes <- Sizes[-remove]
            minSize <- minSize[-remove,-remove]
        }
    }
    ypts <- Bobs[upper.tri(Bobs,diag=TRUE)]
    xpts <- Bhat[upper.tri(Bhat,diag=TRUE)]
    cnts <- minSize[upper.tri(minSize,diag=TRUE)]
    plot(xpts,ypts,type="n",
         xlab=expression(paste("Fitted distances between demes  ",hat(D)[ab],"-(",hat(D)[aa],"+",hat(D)[bb],")/2",sep="")),
         ylab=expression(paste("Observed distances between demes  ",D[ab],"-(",D[aa],"+",D[bb],")/2",sep="")))
    if (!singletons) {
        title(main="Distances between demes\nSingleton demes excluded from plot (not from EEMS)")
    } else {
        title(main="Distances between demes\nGray means a or b has a single individual sampled from")
    }
    abline(a=0,b=1,col="red",lwd=2)
    points(xpts,ypts,col=c("black","gray60")[1+1*(cnts==1)])
    ypts <- Wobs
    xpts <- What
    plot(xpts,ypts,type="n",
         xlab=expression(paste("Fitted distances within demes  ",hat(D)[aa],sep="")),
         ylab=expression(paste("Observed distances within demes  ",D[aa],sep="")))
    if (!singletons) {
        title(main="Distances within demes\nSingleton demes excluded from plot (not from EEMS)")
    } else {
        title(main="Distances within demes\nGray means a has a single individual sampled from")
    }
    abline(a=0,b=1,col="red",lwd=2)
    points(xpts,ypts,col=c("black","gray60")[1+1*(Sizes==1)])
}
