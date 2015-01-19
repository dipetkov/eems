
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
x <- 1-1:numlevels/(numlevels+1)
var.colors <- rgb(cbind(x,x,x))
var.levels <- seq(from=0,to=1,length.out=numlevels)
map.color <- "gray60" ## the color of the geographic map
grid.color <- "gray80" ## the color of the population grid
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
    plot.window(xlim = c(0,1), ylim = range(levels),xaxs = "i",yaxs = "i")
    rect(0, levels[-nlevels], 1, levels[-1], col = colors, border = NA)
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
        eems.levels <- seq(from=minZ,to=maxZ,length=numlevels+1)
    }
    return(eems.levels)
}
center.qrates <- function(Zmrks) {
    minZ <- min(Zmrks)
    maxZ <- max(Zmrks)
    ## The qrates should be centered at about zero
    maxZ2 <- max(maxZ,-minZ)
    maxZ2 <- ceiling(maxZ2*1000)/1000
    minZ <- -maxZ2
    maxZ <- maxZ2
    numlevels <- length(eems.colors)
    if ((maxZ - minZ)>0.001) {
        eems.levels <- seq(from=minZ,to=maxZ,length=numlevels+1)
    } else {
        eems.levels <- seq(from=minZ-0.001/2,to=maxZ+0.001/2,length=numlevels+1)
    }
    return(eems.levels)
}
read.dimns <- function(mcmcpath,longlat,nxmrks=NULL,nymrks=NULL) {
    outer <- scan(paste(mcmcpath,'/outer.txt',sep=''),what=numeric(),quiet=TRUE)
    outer <- matrix(outer,ncol=2,byrow=TRUE)
    if (!longlat) {
        outer <- outer[,c(2,1)]
    }
    xrange <- range(outer[,1])
    yrange <- range(outer[,2])
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
standardize.mmrks <- function(mcmcpath,dimns,longlat) {
    mrates <- scan(paste(mcmcpath,'/mcmcmrates.txt',sep=''),what=numeric(),quiet=TRUE)
    mtiles <- scan(paste(mcmcpath,'/mcmcmtiles.txt',sep=''),what=numeric(),quiet=TRUE)
    xcoord <- scan(paste(mcmcpath,'/mcmcxcoord.txt',sep=''),what=numeric(),quiet=TRUE)
    ycoord <- scan(paste(mcmcpath,'/mcmcycoord.txt',sep=''),what=numeric(),quiet=TRUE)
    if (!longlat) {
        tempic <- xcoord
        xcoord <- ycoord
        ycoord <- tempic
    }
    mrates <- log10(mrates)
    xmrks <- dimns$xmrks
    ymrks <- dimns$ymrks
    nxmrks <- length(xmrks)
    nymrks <- length(ymrks)
    marks <- cbind(rep(xmrks,nymrks),rep(ymrks,each=nxmrks))
    Zmrks <- matrix(0,nxmrks,nymrks)
    niters <- length(mtiles)
    count <- 0
    for (i in 1:niters) {
        curr.mtiles <- mtiles[i]
        curr.mrates <- mrates[(count+1):(count+curr.mtiles)]
        curr.xcoord <- xcoord[(count+1):(count+curr.mtiles)]
        curr.ycoord <- ycoord[(count+1):(count+curr.mtiles)]
        centers <- cbind(curr.xcoord,curr.ycoord)
        distances <- rdist(centers,marks)
        closest <- apply(distances,2,which.min)
        zmrks <- matrix(curr.mrates[closest],nxmrks,nymrks,byrow=FALSE)
        zmrks <- zmrks - mean(zmrks)
        Zmrks <- Zmrks + zmrks
        count <- count + curr.mtiles
    }
    return(list(Zmrks=Zmrks,niters=niters))
}
standardize.mmrks.var <- function(mcmcpath,dimns,Zmmu,longlat) {
    mrates <- scan(paste(mcmcpath,'/mcmcmrates.txt',sep=''),what=numeric(),quiet=TRUE)
    mtiles <- scan(paste(mcmcpath,'/mcmcmtiles.txt',sep=''),what=numeric(),quiet=TRUE)
    xcoord <- scan(paste(mcmcpath,'/mcmcxcoord.txt',sep=''),what=numeric(),quiet=TRUE)
    ycoord <- scan(paste(mcmcpath,'/mcmcycoord.txt',sep=''),what=numeric(),quiet=TRUE)
    if (!longlat) {
        tempic <- xcoord
        xcoord <- ycoord
        ycoord <- tempic
    }
    mrates <- log10(mrates)
    xmrks <- dimns$xmrks
    ymrks <- dimns$ymrks
    nxmrks <- length(xmrks)
    nymrks <- length(ymrks)
    marks <- cbind(rep(xmrks,nymrks),rep(ymrks,each=nxmrks))
    Zmvar <- matrix(0,nxmrks,nymrks)
    niters <- length(mtiles)
    count <- 0
    for (i in 1:niters) {
        curr.mtiles <- mtiles[i]
        curr.mrates <- mrates[(count+1):(count+curr.mtiles)]
        curr.xcoord <- xcoord[(count+1):(count+curr.mtiles)]
        curr.ycoord <- ycoord[(count+1):(count+curr.mtiles)]
        centers <- cbind(curr.xcoord,curr.ycoord)
        distances <- rdist(centers,marks)
        closest <- apply(distances,2,which.min)
        mmrks <- matrix(curr.mrates[closest],nxmrks,nymrks,byrow=FALSE)
        mmrks <- mmrks - mean(mmrks)
        Zmvar <- Zmvar + (mmrks - Zmmu)^2
        count <- count + curr.mtiles
    }
    return(list(Zmvar=Zmvar,niters=niters))
}
mcmc.mrates0 <- function(mcmcpath,dimns,Zmmu,Zmvar,longlat,
                         add.map=FALSE,col.map=map.color,lwd.map=1,
                         add.grid=TRUE,col.grid=grid.color,lwd.grid=1,
                         add.samples=TRUE,col.samples="black",pch.samples=19,cex.samples=1,max.cex.samples=2,
                         add.outline=TRUE,col.outline="red",lwd.outline=3) {
    ipmap <- scan(paste(mcmcpath,'/ipmap.txt',sep=''),what=numeric(),quiet=TRUE)
    coord <- scan(paste(mcmcpath,'/demes.txt',sep=''),what=numeric(),quiet=TRUE)
    edges <- scan(paste(mcmcpath,'/edges.txt',sep=''),what=numeric(),quiet=TRUE)
    outer <- scan(paste(mcmcpath,'/outer.txt',sep=''),what=numeric(),quiet=TRUE)
    coord <- matrix(coord,ncol=2,byrow=TRUE)
    edges <- matrix(edges,ncol=6,byrow=TRUE)
    outer <- matrix(outer,ncol=2,byrow=TRUE)

    if (!longlat) {
        coord <- coord[,c(2,1)]
        outer <- outer[,c(2,1)]
    }
    
    nv <- nrow(edges)
    nn <- ncol(edges)
    nodes <- 1:nv
    sizes <- table(ipmap)
    demes <- as.numeric(names(sizes))
    sizes <- as.numeric(sizes)
    eems.levels <- center.mrates(Zmmu)

    layout(matrix(c(1,2),1,2,byrow=TRUE),c(1,1/4*(dimns$yspan/dimns$xspan)),1,respect=FALSE)

    par(plt = c(0.03,0.97,0.03,0.92),las=1,font.main=1,cex.main=1.5,xpd=TRUE)
    myfilled.contour(dimns$xmrks,dimns$ymrks,Zmmu,asp=1,
                     xlim=dimns$xrange,ylim=dimns$yrange,
                     col=eems.colors,levels=eems.levels,frame.plot=FALSE,
                     plot.axes = { xlab = ''; ylab = '';
                                   if (add.map) {
                                       map(database="worldHires",fill=FALSE,add=TRUE,col=col.map,lwd=lwd.map)
                                   }
                                   if (add.outline) {
                                       for (v in 2:nrow(outer)) {
                                           lines(outer[c(v-1,v),1],outer[c(v-1,v),2],col=col.outline,lwd=lwd.outline)
                                       }
                                   }
                                   if (add.grid) {
                                       for (a in 1:nv) {
                                       for (i in 1:nn) {
                                           b <- edges[a,i]
                                           if (b %in% nodes) {
                                               lines(coord[c(a,b),1],coord[c(a,b),2],col=col.grid,lwd=lwd.grid)
                                           }
                                       } }
                                   }
                                   if (add.samples) {
                                       points(coord[demes,1],coord[demes,2],
                                              col=col.samples,pch=pch.samples,
                                              cex=cex.samples+max.cex.samples*sizes/max(sizes))
                                   }
                               })
    title(main='Effective migration rates m : posterior mean',line=0.5)
    par(plt = c(0.03,0.45,0.03,0.92),cex.axis=1.7,las=1)
    myfilled.legend(dimns$xmrks,dimns$ymrks,Zmrks,asp=1,
                    xlim=dimns$xrange,ylim=dimns$yrange,
                    colors=eems.colors,levels=eems.levels,
                    key.axes = axis(4,tick=FALSE,line=3,hadj=1),
                    key.title = mtext(expression(paste("e"["m"],sep="")),side=3,line=1,cex=1.5))

    min.Zmvar <- min(Zmvar)
    max.Zmvar <- max(Zmvar)
    Zmvar <- (Zmvar - min.Zmvar)/(max.Zmvar - min.Zmvar)
    par(plt = c(0.03,0.97,0.03,0.92),las=1,font.main=1,cex.main=1.5,xpd=TRUE)
    myfilled.contour(dimns$xmrks,dimns$ymrks,Zmvar,asp=1,
                     xlim=dimns$xrange,ylim=dimns$yrange,
                     col=var.colors,levels=var.levels,frame.plot=FALSE,
                     plot.axes = { xlab = ''; ylab = '';
                                   if (add.map) {
                                       map(database="worldHires",fill=FALSE,add=TRUE,col=col.map,lwd=lwd.map)
                                   }
                                   if (add.outline) {
                                       for (v in 2:nrow(outer)) {
                                           lines(outer[c(v-1,v),1],outer[c(v-1,v),2],col=col.outline,lwd=lwd.outline)
                                       }
                                   }                                       
                                   if (add.grid) {
                                       for (a in 1:nv) {
                                       for (i in 1:nn) {
                                           b <- edges[a,i]
                                           if (b %in% nodes) {
                                               lines(coord[c(a,b),1],coord[c(a,b),2],col=col.grid,lwd=lwd.grid)
                                           }
                                       } }
                                   }
                                   if (add.samples) {
                                       points(coord[demes,1],coord[demes,2],
                                              col=col.samples,pch=pch.samples,
                                              cex=cex.samples+max.cex.samples*sizes/max(sizes))
                                   }
                               })
    title(main='Effective migration rates m : posterior variance',line=0.5)
}
mcmc.mrates <- function(mcmcpath,dimns,longlat,
                        add.map=FALSE,col.map=map.color,lwd.map=1,
                        add.grid=TRUE,col.grid=grid.color,lwd.grid=1,
                        add.samples=TRUE,col.samples="black",pch.samples=19,cex.samples=1,max.cex.samples=2,
                        add.outline=TRUE,col.outline="red",lwd.outline=3) {
    print('Plotting effective migration rates m : posterior mean and variance')
    mcmcpath1 <- character()
    for (file in mcmcpath) {
        if (file.exists(paste(file,'/mcmcmtiles.txt',sep=''))&&
            file.exists(paste(file,'/mcmcxcoord.txt',sep=''))&&
            file.exists(paste(file,'/mcmcycoord.txt',sep=''))&&
            file.exists(paste(file,'/mcmcmrates.txt',sep=''))) {
            mcmcpath1 <- c(mcmcpath1,file)
        }
    }
    mcmcpath <- mcmcpath1
    nsimnos <- length(mcmcpath)
    niters <- 0
    if (nsimnos==0) { return(0) }
    Zmmu <- matrix(0,dimns$nxmrks,dimns$nymrks)
    for (file in mcmcpath) {
        rslt <- standardize.mmrks(file,dimns,longlat)
        niters <- niters + rslt$niters
        Zmmu <- Zmmu + rslt$Zmrks
    }
    Zmmu <- Zmmu/niters
    Zmvar <- matrix(0,dimns$nxmrks,dimns$nymrks)
    for (file in mcmcpath) {
        rslt <- standardize.mmrks.var(file,dimns,Zmmu,longlat)
        Zmvar <- Zmvar + rslt$Zmvar
    }
    Zmvar <- Zmvar/(niters-1)
    mcmc.mrates0(mcmcpath[1],dimns,Zmmu,Zmvar,longlat,
                 add.map=add.map,col.map=col.map,lwd.map=lwd.map,
                 add.grid=add.grid,col.grid=col.grid,lwd.grid=lwd.grid,
                 add.samples=add.samples,col.samples=col.samples,
                 pch.samples=pch.samples,cex.samples=cex.samples,max.cex.samples=max.cex.samples,
                 add.outline=add.outline,col.outline=col.outline,lwd.outline=lwd.outline)
}
standardize.qmrks <- function(mcmcpath,dimns,longlat) {
    qrates <- scan(paste(mcmcpath,'/mcmcqrates.txt',sep=''),what=numeric(),quiet=TRUE)
    qtiles <- scan(paste(mcmcpath,'/mcmcqtiles.txt',sep=''),what=numeric(),quiet=TRUE)
    xcoord <- scan(paste(mcmcpath,'/mcmcwcoord.txt',sep=''),what=numeric(),quiet=TRUE)
    ycoord <- scan(paste(mcmcpath,'/mcmczcoord.txt',sep=''),what=numeric(),quiet=TRUE)
    if (!longlat) {
        tempic <- xcoord
        xcoord <- ycoord
        ycoord <- tempic
    }
    qrates <- log10(qrates)
    xmrks <- dimns$xmrks
    ymrks <- dimns$ymrks
    nxmrks <- length(xmrks)
    nymrks <- length(ymrks)
    marks <- cbind(rep(xmrks,nymrks),rep(ymrks,each=nxmrks))
    Zmrks <- matrix(0,nxmrks,nymrks)
    niters <- length(qtiles)
    count <- 0
    for (i in 1:niters) {
        curr.qtiles <- qtiles[i]
        curr.qrates <- qrates[(count+1):(count+curr.qtiles)]
        curr.xcoord <- xcoord[(count+1):(count+curr.qtiles)]
        curr.ycoord <- ycoord[(count+1):(count+curr.qtiles)]
        centers <- cbind(curr.xcoord,curr.ycoord)
        distances <- rdist(centers,marks)
        closest <- apply(distances,2,which.min)
        zmrks <- matrix(curr.qrates[closest],nxmrks,nymrks,byrow=FALSE)
        zmrks <- zmrks - mean(zmrks)
        Zmrks <- Zmrks + zmrks
        count <- count + curr.qtiles
    }
    return(list(Zmrks=Zmrks,niters=niters))
}
standardize.qmrks.var <- function(mcmcpath,dimns,Zqmu,longlat) {
    qrates <- scan(paste(mcmcpath,'/mcmcqrates.txt',sep=''),what=numeric(),quiet=TRUE)
    qtiles <- scan(paste(mcmcpath,'/mcmcqtiles.txt',sep=''),what=numeric(),quiet=TRUE)
    xcoord <- scan(paste(mcmcpath,'/mcmcwcoord.txt',sep=''),what=numeric(),quiet=TRUE)
    ycoord <- scan(paste(mcmcpath,'/mcmczcoord.txt',sep=''),what=numeric(),quiet=TRUE)
    if (!longlat) {
        tempic <- xcoord
        xcoord <- ycoord
        ycoord <- tempic
    }
    qrates <- log10(qrates)
    xmrks <- dimns$xmrks
    ymrks <- dimns$ymrks
    nxmrks <- length(xmrks)
    nymrks <- length(ymrks)
    marks <- cbind(rep(xmrks,nymrks),rep(ymrks,each=nxmrks))
    Zqvar <- matrix(0,nxmrks,nymrks)
    niters <- length(qtiles)
    count <- 0
    for (i in 1:niters) {
        curr.qtiles <- qtiles[i]
        curr.qrates <- qrates[(count+1):(count+curr.qtiles)]
        curr.xcoord <- xcoord[(count+1):(count+curr.qtiles)]
        curr.ycoord <- ycoord[(count+1):(count+curr.qtiles)]
        centers <- cbind(curr.xcoord,curr.ycoord)
        distances <- rdist(centers,marks)
        closest <- apply(distances,2,which.min)
        qmrks <- matrix(curr.qrates[closest],nxmrks,nymrks,byrow=FALSE)
        qmrks <- qmrks - mean(qmrks)
        Zqvar <- Zqvar + (qmrks - Zqmu)^2
        count <- count + curr.qtiles
    }
    return(list(Zqvar=Zqvar,niters=niters))
}
mcmc.qrates0 <- function(mcmcpath,dimns,Zqmu,Zqvar,longlat,
                         add.map=FALSE,col.map=map.color,lwd.map=1,
                         add.grid=TRUE,col.grid=grid.color,lwd.grid=1,
                         add.samples=TRUE,col.samples="black",pch.samples=19,cex.samples=1,max.cex.samples=2,
                         add.outline=TRUE,col.outline="red",lwd.outline=3) {
    ipmap <- scan(paste(mcmcpath,'/ipmap.txt',sep=''),what=numeric(),quiet=TRUE)
    coord <- scan(paste(mcmcpath,'/demes.txt',sep=''),what=numeric(),quiet=TRUE)
    edges <- scan(paste(mcmcpath,'/edges.txt',sep=''),what=numeric(),quiet=TRUE)
    outer <- scan(paste(mcmcpath,'/outer.txt',sep=''),what=numeric(),quiet=TRUE)
    coord <- matrix(coord,ncol=2,byrow=TRUE)
    edges <- matrix(edges,ncol=6,byrow=TRUE)
    outer <- matrix(outer,ncol=2,byrow=TRUE)

    if (!longlat) {
        coord <- coord[,c(2,1)]
        outer <- outer[,c(2,1)]
    }

    nv <- nrow(edges)
    nn <- ncol(edges)
    nodes <- 1:nv
    sizes <- table(ipmap)
    demes <- as.numeric(names(sizes))
    sizes <- as.numeric(sizes)
    eems.levels <- center.qrates(Zqmu)

    layout(matrix(c(1,2),1,2,byrow=TRUE),c(1,1/4*(dimns$yspan/dimns$xspan)),1,respect=FALSE)

    par(plt = c(0.03,0.97,0.03,0.92),las=1,font.main=1,cex.main=1.5,xpd=TRUE)
    myfilled.contour(dimns$xmrks,dimns$ymrks,Zqmu,asp=1,
                     xlim=dimns$xrange,ylim=dimns$yrange,
                     col=eems.colors,levels=eems.levels,frame.plot=FALSE,
                     plot.axes = { xlab = ''; ylab = '';
                                   if (add.map) {
                                       map(database="worldHires",fill=FALSE,add=TRUE,col=col.map,lwd=lwd.map)
                                   }
                                   if (add.outline) {
                                       for (v in 2:nrow(outer)) {
                                           lines(outer[c(v-1,v),1],outer[c(v-1,v),2],col=col.outline,lwd=lwd.outline)
                                       }
                                   }
                                   if (add.grid) {
                                       for (a in 1:nv) {
                                       for (i in 1:nn) {
                                           b <- edges[a,i]
                                           if (b %in% nodes) {
                                               lines(coord[c(a,b),1],coord[c(a,b),2],col=col.grid,lwd=lwd.grid)
                                           }
                                       } }
                                   }
                                   if (add.samples) {
                                       points(coord[demes,1],coord[demes,2],
                                              col=col.samples,pch=pch.samples,
                                              cex=cex.samples+max.cex.samples*sizes/max(sizes))
                                   }
                               })
    title(main='Effective diversity rates q : posterior mean',line=0.5)
    par(plt = c(0.03,0.45,0.03,0.92),cex.axis=1.7,las=1)
    myfilled.legend(dimns$xmrks,dimns$ymrks,Zmrks,asp=1,
                    xlim=dimns$xrange,ylim=dimns$yrange,
                    colors=eems.colors,levels=eems.levels,
                    key.axes = axis(4,tick=FALSE,line=4,hadj=1),
                    key.title = mtext(expression(paste("e"["q"],sep="")),side=3,line=1,cex=1.5))
    
    min.Zqvar <- min(Zqvar)
    max.Zqvar <- max(Zqvar)
    Zqvar <- (Zqvar - min.Zqvar)/(max.Zqvar - min.Zqvar) 
    par(plt = c(0.03,0.97,0.03,0.92),las=1,font.main=1,cex.main=1.5,xpd=TRUE)
    myfilled.contour(dimns$xmrks,dimns$ymrks,Zqvar,asp=1,
                     xlim=dimns$xrange,ylim=dimns$yrange,
                     col=var.colors,levels=var.levels,frame.plot=FALSE,
                     plot.axes = { xlab = ''; ylab = '';
                                   if (add.map) {
                                       map(database="worldHires",fill=FALSE,add=TRUE,col=col.map,lwd=lwd.map)
                                   }
                                   if (add.outline) {
                                       for (v in 2:nrow(outer)) {
                                           lines(outer[c(v-1,v),1],outer[c(v-1,v),2],col=col.outline,lwd=lwd.outline)
                                       }
                                   }                                       
                                   if (add.grid) {
                                       for (a in 1:nv) {
                                       for (i in 1:nn) {
                                           b <- edges[a,i]
                                           if (b %in% nodes) {
                                               lines(coord[c(a,b),1],coord[c(a,b),2],col=col.grid,lwd=lwd.grid)
                                           }
                                       } }
                                   }
                                   if (add.samples) {
                                       points(coord[demes,1],coord[demes,2],
                                              col=col.samples,pch=pch.samples,
                                              cex=cex.samples+max.cex.samples*sizes/max(sizes))
                                   }
                               })
    title(main='Effective diversity rates q : posterior variance',line=0.5)
}
mcmc.qrates <- function(mcmcpath,dimns,longlat,
                        add.map=FALSE,col.map=map.color,lwd.map=1,
                        add.grid=TRUE,col.grid=grid.color,lwd.grid=1,
                        add.samples=TRUE,col.samples="black",pch.samples=19,cex.samples=1,max.cex.samples=2,
                        add.outline=TRUE,col.outline="red",lwd.outline=3) {
    print('Plotting effective diversity rates m : posterior mean and variance')
    mcmcpath1 <- character()
    for (file in mcmcpath) {
        if (file.exists(paste(file,'/mcmcqtiles.txt',sep=''))&&
            file.exists(paste(file,'/mcmcwcoord.txt',sep=''))&&
            file.exists(paste(file,'/mcmczcoord.txt',sep=''))&&
            file.exists(paste(file,'/mcmcqrates.txt',sep=''))) {
            mcmcpath1 <- c(mcmcpath1,file)
        }
    }
    mcmcpath <- mcmcpath1
    nsimnos <- length(mcmcpath)
    niters <- 0
    if (nsimnos==0) { return(0) }
    Zqmu <- matrix(0,dimns$nxmrks,dimns$nymrks)
    for (file in mcmcpath) {
        rslt <- standardize.qmrks(file,dimns,longlat)
        niters <- niters + rslt$niters
        Zqmu <- Zqmu + rslt$Zmrks
    }
    Zqmu <- Zqmu/niters
    Zqvar <- matrix(0,dimns$nxmrks,dimns$nymrks)
    for (file in mcmcpath) {
        rslt <- standardize.qmrks.var(file,dimns,Zqmu,longlat)
        Zqvar <- Zqvar + rslt$Zqvar
    }
    Zqvar <- Zqvar/(niters-1)
    mcmc.qrates0(mcmcpath[1],dimns,Zqmu,Zqvar,longlat,
                 add.map=add.map,col.map=col.map,lwd.map=lwd.map,
                 add.grid=add.grid,col.grid=col.grid,lwd.grid=lwd.grid,
                 add.samples=add.samples,col.samples=col.samples,
                 pch.samples=pch.samples,cex.samples=cex.samples,max.cex.samples=max.cex.samples,
                 add.outline=add.outline,col.outline=col.outline,lwd.outline=lwd.outline)
}
mcmc.mrates.voronoi <- function(mcmcpath,dimns,longlat,max.niter=10,
                                add.grid=TRUE,col.grid=grid.color,lwd.grid=1) {
    ## If there are multiple runs, pick the first one
    mcmcpath <- mcmcpath[1]
    print('Plotting Voronoi tessellation of Effective migration surface')
    print(mcmcpath)
    ipmap <- scan(paste(mcmcpath,'/ipmap.txt',sep=''),what=numeric(),quiet=TRUE)
    coord <- scan(paste(mcmcpath,'/demes.txt',sep=''),what=numeric(),quiet=TRUE)
    edges <- scan(paste(mcmcpath,'/edges.txt',sep=''),what=numeric(),quiet=TRUE)
    outer <- scan(paste(mcmcpath,'/outer.txt',sep=''),what=numeric(),quiet=TRUE)
    coord <- matrix(coord,ncol=2,byrow=TRUE)
    edges <- matrix(edges,ncol=6,byrow=TRUE)
    outer <- matrix(outer,ncol=2,byrow=TRUE)
    if (!longlat) {
        coord <- coord[,c(2,1)]
        outer <- outer[,c(2,1)]
    }
    nv <- nrow(edges)
    nn <- ncol(edges)
    nodes <- 1:nv
    sizes <- table(ipmap)
    demes <- as.numeric(names(sizes))
    sizes <- as.numeric(sizes)
    mrates <- scan(paste(mcmcpath,'/mcmcmrates.txt',sep=''),what=numeric(),quiet=TRUE)
    mtiles <- scan(paste(mcmcpath,'/mcmcmtiles.txt',sep=''),what=numeric(),quiet=TRUE)
    xcoord <- scan(paste(mcmcpath,'/mcmcxcoord.txt',sep=''),what=numeric(),quiet=TRUE)
    ycoord <- scan(paste(mcmcpath,'/mcmcycoord.txt',sep=''),what=numeric(),quiet=TRUE)
    if (!longlat) {
        tempic <- xcoord
        xcoord <- ycoord
        ycoord <- tempic
    }
    mrates <- log10(mrates)
    xrange <- dimns$xrange
    yrange <- dimns$yrange
    niter <- length(mtiles)
    niter <- min(niter,max.niter)
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
        par(plt=c(0.03,0.97,0.03,0.9),las=1,font.main=1,xpd=TRUE)
        plot(0,0,type="n",xlab="",ylab="",axes=FALSE,asp=1,
             xlim=dimns$xrange,ylim=dimns$yrange,
             main=paste('Effective migration rates m : iteration ',i,' (after thinning)',sep=''))
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
            if (add.grid) {
                for (a in 1:nv) {
                for (i in 1:nn) {
                    b <- edges[a,i]
                    if (b %in% nodes) {
                        lines(coord[c(a,b),1],coord[c(a,b),2],col=col.grid,lwd=lwd.grid)
                    }
                } }
            }
        }
        ## Plot the centers
        points(centers,pch=19)
        count <- count + curr.mtiles
    }
    return(list(colors=eems.colors,levels=eems.levels))    
}
mcmc.qrates.voronoi <- function(mcmcpath,dimns,longlat,max.niter=10,
                                add.grid=TRUE,col.grid=grid.color,lwd.grid=1) {
    ## If there are multiple runs, pick the first one
    mcmcpath <- mcmcpath[1]
    print('Plotting Voronoi tessellation of Effective diversity surface')
    print(mcmcpath)
    ipmap <- scan(paste(mcmcpath,'/ipmap.txt',sep=''),what=numeric(),quiet=TRUE)
    coord <- scan(paste(mcmcpath,'/demes.txt',sep=''),what=numeric(),quiet=TRUE)
    edges <- scan(paste(mcmcpath,'/edges.txt',sep=''),what=numeric(),quiet=TRUE)
    outer <- scan(paste(mcmcpath,'/outer.txt',sep=''),what=numeric(),quiet=TRUE)
    coord <- matrix(coord,ncol=2,byrow=TRUE)
    edges <- matrix(edges,ncol=6,byrow=TRUE)
    outer <- matrix(outer,ncol=2,byrow=TRUE)
    if (!longlat) {
        coord <- coord[,c(2,1)]
        outer <- outer[,c(2,1)]
    }
    nv <- nrow(edges)
    nn <- ncol(edges)
    nodes <- 1:nv
    sizes <- table(ipmap)
    demes <- as.numeric(names(sizes))
    sizes <- as.numeric(sizes)
    qrates <- scan(paste(mcmcpath,'/mcmcqrates.txt',sep=''),what=numeric(),quiet=TRUE)
    qtiles <- scan(paste(mcmcpath,'/mcmcqtiles.txt',sep=''),what=numeric(),quiet=TRUE)
    xcoord <- scan(paste(mcmcpath,'/mcmcwcoord.txt',sep=''),what=numeric(),quiet=TRUE)
    ycoord <- scan(paste(mcmcpath,'/mcmczcoord.txt',sep=''),what=numeric(),quiet=TRUE)
    if (!longlat) {
        tempic <- xcoord
        xcoord <- ycoord
        ycoord <- tempic
    }
    qrates <- log10(qrates)
    xrange <- dimns$xrange
    yrange <- dimns$yrange
    niter <- length(qtiles)
    niter <- min(niter,max.niter)
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
        par(plt=c(0.03,0.97,0.03,0.9),las=1,font.main=1,xpd=TRUE)
        plot(0,0,type="n",xlab="",ylab="",axes=FALSE,asp=1,
             xlim=dimns$xrange,ylim=dimns$yrange,
             main=paste('Effective diversity rates q : iteration ',i,' (after thinning)',sep=''))
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
            if (add.grid) {
                for (a in 1:nv) {
                for (i in 1:nn) {
                    b <- edges[a,i]
                    if (b %in% nodes) {
                        lines(coord[c(a,b),1],coord[c(a,b),2],col=col.grid,lwd=lwd.grid)
                    }
                } }
            }
        }
        ## Plot the centers
        points(centers,pch=19)
        count <- count + curr.qtiles
    }
    return(list(colors=eems.colors,levels=eems.levels))    
}
plot.logposterior <- function(mcmcpath) {
    print('Plotting log posterior')
    mcmcpath1 <- character()
    for (file in mcmcpath) {
        if (file.exists(paste(file,'/mcmcpilogl.txt',sep=''))) {
            mcmcpath1 <- c(mcmcpath1,file)
        }
    }
    mcmcpath <- mcmcpath1
    nsimnos <- length(mcmcpath)
    colors <- rep_len(1:8,length.out=nsimnos)
    ltypes <- rep_len(1:3,length.out=nsimnos)
    if (nsimnos==0) { return(0) }
    posteriors <- list()
    nIters <- NULL
    yrange <- NULL
    for (i in 1:length(mcmcpath)) {
        file <- mcmcpath[i]; print(file)
        pilogl <- scan(paste(file,'/mcmcpilogl.txt',sep=''),quiet=TRUE)
        pilogl <- matrix(pilogl,ncol=2,byrow=TRUE)
        posterior <- pilogl[,1] + pilogl[,2]
        posteriors[[i]] <- posterior
        yrange <- range(c(yrange,posterior))
        nIters <- max(nIters,length(posterior))
    }
    plot(c(1,nIters),yrange,type="n",xlab="iteration (after thinning)",ylab="log posterior")
    for (i in 1:length(mcmcpath)) {
        posterior <- posteriors[[i]]
        lines(1:length(posterior),posterior,col=colors[i],lty=ltypes[i],lwd=2)
    }
    legend("bottomleft",legend=1:nsimnos,col=colors[1:nsimnos],lty=ltypes[1:nsimnos],lwd=2,bty="n",horiz=TRUE)
}
dist.scatterplot <- function(mcmcpath,remove.singletons=TRUE) {
    print('Plotting average dissimilarities within and between demes')
    mcmcpath1 <- character()
    for (file in mcmcpath) {
        if (file.exists(paste(file,'/rdistJtDobsJ.txt',sep=''))&&
            file.exists(paste(file,'/rdistJtDhatJ.txt',sep=''))&&
            file.exists(paste(file,'/rdistoDemes.txt',sep=''))) {
            mcmcpath1 <- c(mcmcpath1,file)
        }
    }
    mcmcpath <- mcmcpath1
    nsimnos <- length(mcmcpath)
    if (nsimnos==0) { return(0) }
    ## list of observed demes, with number of samples taken from each
    ## each row specifies: x coordinate, y coordinate, n samples
    oDemes <- scan(paste(mcmcpath[1],'/rdistoDemes.txt',sep=''),quiet=TRUE)
    oDemes <- matrix(oDemes,ncol=3,byrow=TRUE)
    if (nsimnos>1) {
        for (simno in 2:nsimnos) {
            file <- mcmcpath[simno]; print(file)
            oDemes2 <- scan(paste(file,'/rdistoDemes.txt',sep=''),quiet=TRUE)
            oDemes2 <- matrix(oDemes2,ncol=3,byrow=TRUE)
            if ((length(oDemes2)!=length(oDemes))||
                (sum(oDemes2!=oDemes)>0)) {
                stop('dist.scatterplot: mcmc results for different population graphs.')
            }
        }
    }
    Sizes <- oDemes[,3]
    nPops <- length(Sizes)
    matSize <- matrix(Sizes,nPops,nPops)
    minSize <- pmin(matSize,t(matSize))
    JtDobsJ <- matrix(0,nPops,nPops)
    JtDhatJ <- matrix(0,nPops,nPops)
    for (file in mcmcpath) {
        print(file)
        JtDobsJ <- JtDobsJ + as.matrix(read.table(paste(file,'/rdistJtDobsJ.txt',sep=''),header=FALSE))
        JtDhatJ <- JtDhatJ + as.matrix(read.table(paste(file,'/rdistJtDhatJ.txt',sep=''),header=FALSE))
    }
    JtDobsJ <- JtDobsJ/nsimnos
    JtDhatJ <- JtDhatJ/nsimnos
    if (remove.singletons) {
        ## Remove singleton locations
        remove <- which(Sizes<=1)
        if (length(remove)) {
            JtDobsJ <- JtDobsJ[-remove,-remove]
            JtDhatJ <- JtDhatJ[-remove,-remove]
            minSize <- minSize[-remove,-remove]
            Sizes <- Sizes[-remove]
            nPops <- length(Sizes)
        }
    }
    if (nPops<2) {
        print('There is a single observed deme')
        return (0)
    }
    Wobs <- diag(JtDobsJ)
    What <- diag(JtDhatJ)
    ones <- matrix(1,nPops,1)
    Bobs <- JtDobsJ - (Wobs%*%t(ones) + ones%*%t(Wobs))/2
    Bhat <- JtDhatJ - (What%*%t(ones) + ones%*%t(What))/2
    ypts <- Bobs[upper.tri(Bobs,diag=FALSE)]
    xpts <- Bhat[upper.tri(Bhat,diag=FALSE)]
    cnts <- minSize[upper.tri(minSize,diag=FALSE)]
    par(font.main=1)
    plot(xpts,ypts,col=c("black","gray60")[1+1*(cnts==1)],
         xlab=expression(paste("Fitted dissimilarity between demes,  ",hat(D)[ab]," - (",hat(D)[aa],"+",hat(D)[bb],")/2",sep="")),
         ylab=expression(paste("Observed dissimilarity between demes,  ",D[ab]," - (",D[aa],"+",D[bb],")/2",sep="")))
    if (remove.singletons) {
        title(main="Dissimilarities between demes\nSingleton demes excluded from plot (not from EEMS)")
    } else {
        title(main="Dissimilarities between demes\nGray means a or b has a single individual sampled from")
    }
    ypts <- Wobs
    xpts <- What
    plot(xpts,ypts,col=c("black","gray60")[1+1*(Sizes==1)],
         xlab=expression(paste("Fitted dissimilarity within demes,  ",hat(D)[aa],sep="")),
         ylab=expression(paste("Observed dissimilarity within demes,  ",D[aa],sep="")))
    if (remove.singletons) {
        title(main="Dissimilarities within demes\nSingleton demes excluded from plot (not from EEMS)")
    } else {
        title(main="Dissimilarities within demes\nGray means a has a single individual sampled from")
    }
}
####################################################################################
## eems.plots takes three required arguments:
##   mcmcpath: one or several output directories (for the same dataset)
##   plotpath: filename of figures to generate
##     The following figures are created by default:
##     plotpath-mrates01.png: effective migration rates
##     plotpath-qrates01.png: effective diversity rates
##     plotpath-rist01/02.png: fitted vs observed distances
##   longlat (TRUE or FALSE): are the coordinates ordered longitude/latitude or not?
## eems.plots take a variety of optional arguments:
##   plot.width and plot.height: width and height of the graphics region (in inches)
##   res: resolution of png figures (generated with the bitmap function)
##   add.map: add 'worldHires' map (using the mapdata package)?
##   add.grid: add triangular population grid?
##   add.samples: add samples to their assigned location in the grid?
##   add.outline: add the habitat ring (as declared in the .outer file)?
##   col.map/col.grid/col.samples/col.outline: specify the colors
##   lwd.max/lwd.grid/lwd.outline: specify the line width
##   pch.samples: specify the character
##   cex.samples: specify the character size
##   max.cex.samples: some demes might be assigned more samples than others.
##     If max.cex.samples>0, then demes with more samples will also have bigger size.
##     If the sampling is uneven, then max.cex.samples>0 will underline this fact.
####################################################################################
eems.plots <- function(mcmcpath,plotpath,longlat,
                       plot.width=0,plot.height=0,res=600,
                       add.map=FALSE,col.map=map.color,lwd.map=1,
                       add.grid=TRUE,col.grid=grid.color,lwd.grid=1,
                       add.samples=TRUE,col.samples="black",pch.samples=19,cex.samples=1,max.cex.samples=2,
                       add.outline=TRUE,col.outline="red",lwd.outline=3) {
    
    if (add.map) { library(mapdata) }
    if (max.cex.samples<0) { max.cex.samples = 0 }
    
    mcmcpath1 <- character()
    for (file in mcmcpath) {
        if (file.exists(paste(file,'/ipmap.txt',sep=''))&&
            file.exists(paste(file,'/demes.txt',sep=''))&&
            file.exists(paste(file,'/edges.txt',sep=''))) {
            mcmcpath1 <- c(mcmcpath1,file)
        }
    }
    mcmcpath <- mcmcpath1
    nsimnos <- length(mcmcpath)
    if (nsimnos==0) { return(0) }
    
    dimns <- read.dimns(mcmcpath[1],longlat)
    if ((plot.height<=0)||(plot.width<=0)) {
        plot.height <- 1.1*dimns$yspan
        plot.width <- 1.25*dimns$xspan
    }

    bitmap(paste(plotpath,'-mrates%02d.png',sep=''),type='png16m',res=res,
           height=plot.height,width=plot.width,units='in')
    mcmc.mrates(mcmcpath,dimns,longlat,
                add.grid=add.grid,add.map=add.map,
                col.grid=col.grid,col.map=col.map,
                lwd.grid=lwd.grid,lwd.map=lwd.map,
                add.samples=add.samples,col.samples=col.samples,
                pch.samples=pch.samples,cex.samples=cex.samples,max.cex.samples=max.cex.samples,
                add.outline=add.outline,col.outline=col.outline,lwd.outline=lwd.outline)
    dev.off( )

    bitmap(paste(plotpath,'-qrates%02d.png',sep=''),type='png16m',res=res,
           height=plot.height,width=plot.width,units='in')
    mcmc.qrates(mcmcpath,dimns,longlat,
                add.grid=add.grid,add.map=add.map,
                col.grid=col.grid,col.map=col.map,
                lwd.grid=lwd.grid,lwd.map=lwd.map,
                add.samples=add.samples,col.samples=col.samples,
                pch.samples=pch.samples,cex.samples=cex.samples,max.cex.samples=max.cex.samples,
                add.outline=add.outline,col.outline=col.outline,lwd.outline=lwd.outline)
    dev.off( )

    plot.height <- 5
    plot.width <- 6
    bitmap(paste(plotpath,'-rdist%02d.png',sep=''),type='png16m',res=res,
           height=plot.height,width=plot.width,units='in')
    dist.scatterplot(mcmcpath)
    dev.off( )
    bitmap(paste(plotpath,'-pilogl%02d.png',sep=''),type='png16m',res=res,
           height=plot.height,width=plot.width,units='in')
    plot.logposterior(mcmcpath)
    dev.off( )
}
