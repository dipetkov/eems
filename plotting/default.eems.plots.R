
###### Required R packages ######

## The script requires several packages that handle geospatial data:
require(rgdal) ## spTransform
require(rgeos) ## gDifference
require(raster) ## filledContour
require(fields) ## the rdist/rdist.earth functions compute Euclidean/Great circle distances
##require(deldir) ## the deldir function calculates Voronoi tessellation
##require(rworldmap) ## the map database required to add geographic maps
##require(rworldxtra) ## and the high resolution world map

## The distance metric used by EEMS (Euclidean distance by default)
which.dist.metric <- function(mcmcpath) {
    dist.metric <- "euclidean"
    Lines = readLines(paste(mcmcpath,"/eemsrun.txt",sep=""))
    nLines = length(Lines)
    for (i in seq(nLines)) {
        ## Check for a line that specifies "distance = ..."
        str = gsub("\\s","",Lines[i])
        trm = strsplit(str,"distance=")[[1]]
        if (length(trm)==2) { dist.metric = trm[2] }
    }
    print(paste('Distance metric: ',dist.metric))
    return (dist.metric)
}
pairwise.dist <- function(x,y,dist.metric) {
    if (dist.metric=="greatcirc") {
        ## Great circle distance:
        return (rdist.earth(x,y)) ## miles=FALSE, R=6371
    } else {
        ## Euclidean distance:
        return (rdist(x,y))
    }
}

###### Define the default color schemes ######

## This is the default Dark Orange to Blue color scheme, with "white" as the midpoint color.
## It combines two color schemes from the dichromat package, which itself is based on
## a collection of color schemes for scientific data graphics:
## Light A and Bartlein PJ (2004). The End of the Rainbow? Color Schemes for Improved Data
## Graphics. EOS Transactions of the American Geophysical Union, 85(40), 385.
## Also see http://geog.uoregon.edu/datagraphics/color_scales.htm

## The default Dark Orange to Blue color scheme:
default.eems.colors <- function( ) {
    eems.colors <- c("#994000","#CC5800","#FF8F33","#FFAD66","#FFCA99","#FFE6CC", ## orange sequence
                     "#FFFFFF",                                                   ## white
                     "#CCFDFF","#99F8FF","#66F0FF","#33E4FF","#00AACC","#007A99") ## blue sequence
    return (eems.colors)
}
default.var.colors <- function( ) {
    x <- 1-1:13/(13+1)  ## There are 13 levels in the default scheme
    var.colors <- rgb(cbind(x,x,x))  ## gray sequence
    return (var.colors)
}
## Default color of the geographic map
default.map.color <- function( ) { return ("gray60") }
## Default color of the EEMS population grid
default.grid.color <- function( ) { return ("gray80") }

###### Define routines to normalize the rates before plotting them ######

mrates.levels <- function(Zvals) {
    minZ <- min(Zvals)
    maxZ <- max(Zvals)
    eems.colors <- default.eems.colors()
    numlevels <- length(eems.colors)
    ## First check whether the Z values are in the range (-2.5,+2.5), since
    ## it might be nice to have the migration rates always plotted on the same scale
    if ((minZ > -2.5)&&(maxZ < +2.5)) {
        minZ <- -2.5
        maxZ <- +2.5
        step <- (maxZ-minZ)/numlevels
        eems.levels <- seq(from=minZ,to=maxZ,by=step)
    } else {
        maxZ2 <- max(maxZ,-minZ)
        maxZ2 <- ceiling(maxZ2*1000)/1000
        minZ <- -maxZ2
        maxZ <- +maxZ2
        eems.levels <- seq(from=minZ,to=maxZ,length=numlevels+1)
    }
    return(eems.levels)
}
qrates.levels <- function(Zvals) {
    minZ <- min(Zvals)
    maxZ <- max(Zvals)
    eems.colors <- default.eems.colors()
    numlevels <- length(eems.colors)
    maxZ2 <- max(maxZ,-minZ)
    maxZ2 <- ceiling(maxZ2*1000)/1000
    minZ <- -maxZ2
    maxZ <- +maxZ2
    ## Currently, there is no default range for the diversity rates
    ## unless there is almost no variation, i.e., (maxZ-minZ)<0.001
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
    xlim <- range(outer[,1])
    ylim <- range(outer[,2])
    aspect <- (diff(ylim)/diff(xlim))/cos((mean(ylim) * pi)/180)
    ## Choose the number of interpolation in each direction
    if (is.null(nxmrks)&&is.null(nymrks)) {
        if (aspect>1) {
            nxmrks <- 100
            nymrks <- round(nxmrks*aspect)
        } else {
            nymrks <- 100
            nxmrks <- round(nymrks/aspect)
        }
    }
    ## The interpolation points are equally spaced
    xmrks <- seq(from=xlim[1],to=xlim[2],length=nxmrks)
    ymrks <- seq(from=ylim[1],to=ylim[2],length=nymrks)
    marks <- cbind(rep(xmrks,times=nymrks),rep(ymrks,each=nxmrks))
    dist.metric <- which.dist.metric(mcmcpath)
    return(list(nxmrks=nxmrks,xmrks=xmrks,xlim=xlim,xspan=diff(xlim),
                nymrks=nymrks,ymrks=ymrks,ylim=ylim,yspan=diff(ylim),
                marks=marks,aspect=aspect,
                outer=outer,dist.metric=dist.metric))
}
read.edges <- function(mcmcpath) {
    edges <- read.table(paste(mcmcpath,'/edges.txt',sep=''),colClasses=numeric())
    edges <- as.matrix(edges)
    ## Previously EEMS output the edges with one vertex per line and
    ## the six neighbors of each vertex listed in order
    ## Currently EEMS outputs the edges with one edge per line (as a
    ## pair of vertices) which is more general
    if (ncol(edges)==6) {
        ## Convert the old format to the new format
        edges0 <- edges
        edges <- matrix(0,nrow=sum(edges0>0),ncol=2)
        nv <- nrow(edges0)
        nn <- ncol(edges0)
        nodes <- 1:nv
        e <- 0
        for (a in 1:nv) {
        for (i in 1:nn) {
            b <- edges0[a,i]
            if (b %in% nodes) {
                e <- e + 1
                edges[e,1] = a
                edges[e,2] = b
            }
        } }
    }
    return(edges)
}
read.voronoi <- function(mcmcpath,longlat,is.mrates) {
    if (is.mrates) {
        rates <- scan(paste(mcmcpath,'/mcmcmrates.txt',sep=''),what=numeric(),quiet=TRUE)
        tiles <- scan(paste(mcmcpath,'/mcmcmtiles.txt',sep=''),what=numeric(),quiet=TRUE)
        xseed <- scan(paste(mcmcpath,'/mcmcxcoord.txt',sep=''),what=numeric(),quiet=TRUE)
        yseed <- scan(paste(mcmcpath,'/mcmcycoord.txt',sep=''),what=numeric(),quiet=TRUE)
    } else {
        rates <- scan(paste(mcmcpath,'/mcmcqrates.txt',sep=''),what=numeric(),quiet=TRUE)
        tiles <- scan(paste(mcmcpath,'/mcmcqtiles.txt',sep=''),what=numeric(),quiet=TRUE)
        xseed <- scan(paste(mcmcpath,'/mcmcwcoord.txt',sep=''),what=numeric(),quiet=TRUE)
        yseed <- scan(paste(mcmcpath,'/mcmczcoord.txt',sep=''),what=numeric(),quiet=TRUE)
    }
    if (!longlat) {
        tempi <- xseed
        xseed <- yseed
        yseed <- tempi
    }
    rates <- log10(rates)
    return(list(rates=rates,tiles=tiles,xseed=xseed,yseed=yseed))
}
compute.contour.vals <- function(dimns,seeds,rates,use.weighted.mean=TRUE) {
    ## Here 'seeds' stores the generator seeds of a Voronoi tessellation
    ## and 'rates' stores the log10-transformed rates of the tiles.
    ## If there are C seeds in the partition, then 'seeds' is a matrix
    ## with C rows and 2 columns and 'rates' is a vector with C elements
    distances <- pairwise.dist(dimns$marks,seeds,dimns$dist.metric)
    closest <- apply(distances,1,which.min)
    if (use.weighted.mean) {
        zvals <- matrix(rates[closest],dimns$nxmrks,dimns$nymrks,byrow=FALSE)
        zvals <- zvals - mean(zvals)
    } else {
        rates <- rates - mean(rates)
        zvals <- matrix(rates[closest],dimns$nxmrks,dimns$nymrks,byrow=FALSE)
    }
    return(zvals)
}
standardize.rates <- function(mcmcpath,dimns,longlat,is.mrates) {
    voronoi <- read.voronoi(mcmcpath,longlat,is.mrates)
    rates <- voronoi$rates
    tiles <- voronoi$tiles
    xseed <- voronoi$xseed
    yseed <- voronoi$yseed
    Zvals <- matrix(0,dimns$nxmrks,dimns$nymrks)
    niter <- length(tiles)
    count <- 0
    for (i in 1:niter) {
        now.tiles <- tiles[i]
        now.rates <- rates[(count+1):(count+now.tiles)]
        now.xseed <- xseed[(count+1):(count+now.tiles)]
        now.yseed <- yseed[(count+1):(count+now.tiles)]
        now.seeds <- cbind(now.xseed,now.yseed)
        zvals <- compute.contour.vals(dimns,now.seeds,now.rates)
        Zvals <- Zvals + zvals
        count <- count + now.tiles
    }
    return(list(Zvals=Zvals,niter=niter))
}
standardize.rates.var <- function(mcmcpath,dimns,Zmean,longlat,is.mrates) {
    voronoi <- read.voronoi(mcmcpath,longlat,is.mrates)
    rates <- voronoi$rates
    tiles <- voronoi$tiles
    xseed <- voronoi$xseed
    yseed <- voronoi$yseed
    Zvar <- matrix(0,dimns$nxmrks,dimns$nymrks)
    niter <- length(tiles)
    count <- 0
    for (i in 1:niter) {
        now.tiles <- tiles[i]
        now.rates <- rates[(count+1):(count+now.tiles)]
        now.xseed <- xseed[(count+1):(count+now.tiles)]
        now.yseed <- yseed[(count+1):(count+now.tiles)]
        now.seeds <- cbind(now.xseed,now.yseed)
        zvals <- compute.contour.vals(dimns,now.seeds,now.rates)
        Zvar <- Zvar + (zvals - Zmean)^2
        count <- count + now.tiles
    }
    return(list(Zvar=Zvar,niter=niter))
}
filled.countour.axes <- function(mcmcpath,longlat,plot.params) {
    ipmap <- scan(paste(mcmcpath,'/ipmap.txt',sep=''),what=numeric(),quiet=TRUE)
    demes <- scan(paste(mcmcpath,'/demes.txt',sep=''),what=numeric(),quiet=TRUE)
    outer <- scan(paste(mcmcpath,'/outer.txt',sep=''),what=numeric(),quiet=TRUE)
    demes <- matrix(demes,ncol=2,byrow=TRUE)
    outer <- matrix(outer,ncol=2,byrow=TRUE)
    edges <- read.edges(mcmcpath)
    if (!longlat) {
        demes <- demes[,c(2,1)]
        outer <- outer[,c(2,1)]
    }
    sizes <- table(ipmap)
    alpha <- as.numeric(names(sizes))
    sizes <- as.numeric(sizes)
    ## The filledContour fills a rectangular plot; now color the habitat exterior white
    boundary <- SpatialPolygons(list(Polygons(list(Polygon(outer, hole = FALSE)),"1")),
                                proj4string=CRS(plot.params$proj.in))
    boundary <- SpatialPolygonsDataFrame(boundary,data.frame(id="1"))
    boundary <- spTransform(boundary,CRS=CRS(plot.params$proj.out))
    exterior <- gDifference(gEnvelope(boundary),boundary)
    if (!is.null(exterior)) {
        plot(exterior,col="white",border="white",add=TRUE)
    }
    if (plot.params$add.map) {
        map <- getMap(resolution="high")
        map <- spTransform(map,CRS=CRS(plot.params$proj.out))
        map <- gIntersection(map,boundary,byid=TRUE)
        plot(map,col=NA,border=plot.params$col.map,lwd=plot.params$lwd.map,add=TRUE)
    }
    if (plot.params$add.grid) {
        segments <- list()
        for (e in 1:nrow(edges)) {
            segments[[e]] <- Line(cbind(demes[edges[e,],1],demes[edges[e,],2]))
        }
        segments <- SpatialLines(list(Lines(segments,ID="a")),proj4string=CRS(plot.params$proj.in))
        segments <- spTransform(segments,CRS=CRS(plot.params$proj.out))
        lines(segments,col=plot.params$col.grid,lwd=plot.params$lwd.grid)
    }
    if (plot.params$add.outline) {
        plot(boundary,col=NA,border=plot.params$col.outline,lwd=plot.params$lwd.outline,add=TRUE)
    }
    if (plot.params$add.samples) {
        samples <- SpatialPoints(cbind(demes[alpha,1],demes[alpha,2]),proj4string=CRS(plot.params$proj.in))
        samples <- spTransform(samples,CRS=CRS(plot.params$proj.out))
        points(samples,col=plot.params$col.samples,pch=plot.params$pch.samples,
               cex=plot.params$cex.samples+plot.params$max.cex.samples*sizes/max(sizes))
    }
}
one.eems.contour <- function(mcmcpath,dimns,Zmean,Zvar,longlat,plot.params,is.mrates) {
    eems.colors <- default.eems.colors( )
    if (is.mrates) {
        eems.levels <- mrates.levels(Zmean)
        main.title <- "Effective migration rates m : posterior mean"
        var.title <- "Effective migration rates m : posterior variance"
        key.title <- expression(paste("e"["m"],sep=""))
    } else {
        eems.levels <- qrates.levels(Zmean)
        main.title <- "Effective diversity rates q : posterior mean"
        var.title <- "Effective diversity rates q : posterior variance"
        key.title <- expression(paste("e"["q"],sep=""))
    }
    rr <- flip(raster(t(Zmean),
                      xmn=dimns$xlim[1],xmx=dimns$xlim[2],
                      ymn=dimns$ylim[1],ymx=dimns$ylim[2]),direction='y')
    projection(rr) <- CRS(plot.params$proj.in)
    filledContour(projectRaster(rr,crs=CRS(plot.params$proj.out)),
                  main=main.title,font.main=1,asp=1,
                  col=eems.colors,levels=eems.levels,frame.plot=FALSE,
                  key.axes = axis(4,tick=FALSE,hadj=1,line=3,cex.axis=1.5),
                  key.title = mtext(key.title,side=3,cex=1.5,font=1),
                  plot.axes = filled.countour.axes(mcmcpath,longlat,plot.params))
    
    min.Zvar <- min(Zvar)
    max.Zvar <- max(Zvar)
    Zvar <- (Zvar - min.Zvar)/(max.Zvar - min.Zvar)
    var.colors <- default.var.colors()
    var.levels <- seq(from=0,to=1,length.out=length(var.colors))

    rr <- flip(raster(t(Zvar),
                      xmn=dimns$xlim[1],xmx=dimns$xlim[2],
                      ymn=dimns$ylim[1],ymx=dimns$ylim[2]),direction='y')
    projection(rr) <- CRS(plot.params$proj.in)
    filledContour(projectRaster(rr,crs=CRS(plot.params$proj.out)),
                  main=var.title,font.main=1,asp=1,
                  col=var.colors,levels=var.levels,frame.plot=FALSE,
                  key.axes = axis(4,tick=FALSE,hadj=1,line=3,cex.axis=1.5),
                  key.title = mtext(key.title,side=3,cex=1.5),
                  plot.axes = filled.countour.axes(mcmcpath,longlat,plot.params))
}
average.eems.contours <- function(mcmcpath,dimns,longlat,plot.params,is.mrates) {
    mcmcpath1 <- character( )
    if (is.mrates) {
        print('Plotting effective migration rates m : posterior mean and variance')
        for (path in mcmcpath) {
            Files <- paste(path,c('/mcmcmtiles.txt','/mcmcmrates.txt',
                                  '/mcmcxcoord.txt','/mcmcycoord.txt'),sep='')
            if (min(file.exists(Files))) { mcmcpath1 <- c(mcmcpath1,path) }
        }
    } else {
        print('Plotting effective diversity rates q : posterior mean and variance')
        for (path in mcmcpath) {
            Files <- paste(path,c('/mcmcqtiles.txt','/mcmcqrates.txt',
                                  '/mcmcwcoord.txt','/mcmczcoord.txt'),sep='')
            if (min(file.exists(Files))) { mcmcpath1 <- c(mcmcpath1,path) }
        }
    }
    mcmcpath <- mcmcpath1
    nsimnos <- length(mcmcpath)
    niter <- 0
    if (nsimnos==0) { return(0) }
    Zmean <- matrix(0,dimns$nxmrks,dimns$nymrks)
    for (path in mcmcpath) {
        rslt <- standardize.rates(path,dimns,longlat,is.mrates)
        niter <- niter + rslt$niter
        Zmean <- Zmean + rslt$Zvals
    }
    Zmean <- Zmean/niter
    Zvar <- matrix(0,dimns$nxmrks,dimns$nymrks)
    for (path in mcmcpath) {
        rslt <- standardize.rates.var(path,dimns,Zmean,longlat,is.mrates)
        Zvar <- Zvar + rslt$Zvar
    }
    Zvar <- Zvar/(niter-1)
    one.eems.contour(mcmcpath[1],dimns,Zmean,Zvar,longlat,plot.params,is.mrates)
}
## If there are multiple runs, pick the first one and create at most max.niter Voronoi diagrams
## This function is mainly for testing purposes, so no need to create a plot for each iteration
voronoi.diagram <- function(mcmcpath,dimns,longlat,plot.params,is.mrates,max.niter=10) {
    require(deldir)
    mcmcpath <- mcmcpath[1]
    print('Plotting Voronoi tessellation of estimated effective rates')
    print(mcmcpath)
    voronoi <- read.voronoi(mcmcpath,longlat,is.mrates)
    rates <- voronoi$rates
    tiles <- voronoi$tiles
    xseed <- voronoi$xseed
    yseed <- voronoi$yseed
    if (is.mrates) {
        main.title <- 'Effective migration rates m'
        eems.levels <- mrates.levels(rates)
    } else {
        main.title <- 'Effective diversity rates q'
        eems.levels <- qrates.levels(rates)
    }
    count <- 0
    niter <- min(length(tiles),max.niter)
    eems.colors <- default.eems.colors( )
    for (i in 1:niter) {
        now.tiles <- tiles[i]
        now.rates <- rates[(count+1):(count+now.tiles)]
        now.xseed <- xseed[(count+1):(count+now.tiles)]
        now.yseed <- yseed[(count+1):(count+now.tiles)]
        ## Standardize the log-transformed rates, without taking into account
        ## the relative size of the tiles (this is hard to do without a Zvals grid)
        now.rates <- now.rates - mean(now.rates)
        L <- length(eems.levels)
        indices <- which(now.rates<eems.levels[1])
        now.rates[indices] <- 0.999*eems.levels[1]
        indices <- which(now.rates>eems.levels[L])
        now.rates[indices] <- 0.999*eems.levels[L]
        now.seeds <- cbind(now.xseed,now.yseed)
        now.colors <- character( )
        plot(0,0,type="n",axes=FALSE,xlab="",ylab="",xlim=dimns$xlim,ylim=dimns$ylim,asp=dimns$asp,
             main=paste(main.title,' : iteration ',i,' (after thinning)',sep=''))
        if (now.tiles==1) {
            ## There is only one tile
            tile.color <- eems.colors[round(L/2)]
            now.colors <- c(now.colors,tile.color)
            polygon(dimns$xlim,dimns$ylim,col=tile.color,border=FALSE)
        } else {
            ## Plot each tile in turn (as a polygon)
            Voronoi <- deldir(now.xseed,now.yseed,rw=c(dimns$xlim,dimns$ylim))
            tilelist <- tile.list(Voronoi)
            for (c in 1:now.tiles) {
                tile.color <- eems.colors[ max((1:L)[eems.levels<now.rates[c]]) ]
                now.colors <- c(now.colors,tile.color)
                polygon(tilelist[[c]]$x,tilelist[[c]]$y,col=tile.color,border=FALSE)
            }
            filled.countour.axes(mcmcpath,longlat,plot.params)
        }
        points(now.seeds,pch=4,col="red")
        count <- count + now.tiles
    }
    return(list(colors=eems.colors,levels=eems.levels))    
}
plot.logposterior <- function(mcmcpath) {
    print('Plotting log posterior')
    mcmcpath1 <- character()
    for (path in mcmcpath) {
        if (file.exists(paste(path,'/mcmcpilogl.txt',sep=''))) {
            mcmcpath1 <- c(mcmcpath1,path)
        }
    }
    mcmcpath <- mcmcpath1
    nsimnos <- length(mcmcpath)
    colors <- rep_len(1:8,length.out=nsimnos)
    ltypes <- rep_len(1:3,length.out=nsimnos)
    if (nsimnos==0) { return(0) }
    posteriors <- list()
    ylim <- NULL
    niter <- NULL
    for (i in 1:length(mcmcpath)) {
        path <- mcmcpath[i]; print(path)
        pilogl <- scan(paste(path,'/mcmcpilogl.txt',sep=''),quiet=TRUE)
        pilogl <- matrix(pilogl,ncol=2,byrow=TRUE)
        posterior <- pilogl[,1] + pilogl[,2]
        posteriors[[i]] <- posterior
        ylim <- range(c(ylim,posterior))
        niter <- max(niter,length(posterior))
    }
    plot(c(1,niter),ylim,type="n",xlab="iteration (after thinning)",ylab="log posterior")
    for (i in 1:length(mcmcpath)) {
        posterior <- posteriors[[i]]
        lines(1:length(posterior),posterior,col=colors[i],lty=ltypes[i],lwd=2)
    }
    legend("bottomleft",legend=1:nsimnos,col=colors[1:nsimnos],lty=ltypes[1:nsimnos],lwd=2,bty="n",horiz=TRUE)
}
dist.scatterplot <- function(mcmcpath,remove.singletons=TRUE) {
    print('Plotting average dissimilarities within and between demes')
    mcmcpath1 <- character()
    for (path in mcmcpath) {
        if (file.exists(paste(path,'/rdistJtDobsJ.txt',sep=''))&&
            file.exists(paste(path,'/rdistJtDhatJ.txt',sep=''))&&
            file.exists(paste(path,'/rdistoDemes.txt',sep=''))) {
            mcmcpath1 <- c(mcmcpath1,path)
        }
    }
    mcmcpath <- mcmcpath1
    nsimnos <- length(mcmcpath)
    if (nsimnos==0) { return(0) }
    ## List of observed demes, with number of samples taken collected
    ## Each row specifies: x coordinate, y coordinate, n samples
    oDemes <- scan(paste(mcmcpath[1],'/rdistoDemes.txt',sep=''),quiet=TRUE)
    oDemes <- matrix(oDemes,ncol=3,byrow=TRUE)
    if (nsimnos>1) {
        for (simno in 2:nsimnos) {
            path <- mcmcpath[simno]; print(path)
            oDemes2 <- scan(paste(path,'/rdistoDemes.txt',sep=''),quiet=TRUE)
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
    for (path in mcmcpath) {
        print(path)
        JtDobsJ <- JtDobsJ + as.matrix(read.table(paste(path,'/rdistJtDobsJ.txt',sep=''),header=FALSE))
        JtDhatJ <- JtDhatJ + as.matrix(read.table(paste(path,'/rdistJtDhatJ.txt',sep=''),header=FALSE))
    }
    JtDobsJ <- JtDobsJ/nsimnos
    JtDhatJ <- JtDhatJ/nsimnos
    if (remove.singletons) {
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
        print('Need at least two observed demes to plot pairwise differences')
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
    plot(xpts,ypts,col=c("black","gray60")[1+1*(cnts==1)],
         xlab=expression(paste("Fitted dissimilarity between demes,  ",hat(D)[ab]," - (",hat(D)[aa],"+",hat(D)[bb],")/2",sep="")),
         ylab=expression(paste("Observed dissimilarity between demes,  ",D[ab]," - (",D[aa],"+",D[bb],")/2",sep="")))
    if (remove.singletons) {
        title(main="Dissimilarities between demes\nSingleton demes, if any, excluded from plot (not from EEMS)")
    } else {
        title(main="Dissimilarities between demes\nGray means a or b has a single individual sampled from")
    }
    ypts <- Wobs
    xpts <- What
    plot(xpts,ypts,col=c("black","gray60")[1+1*(Sizes==1)],
         xlab=expression(paste("Fitted dissimilarity within demes,  ",hat(D)[aa],sep="")),
         ylab=expression(paste("Observed dissimilarity within demes,  ",D[aa],sep="")))
    if (remove.singletons) {
        title(main="Dissimilarities within demes\nSingleton demes, if any, excluded from plot (not from EEMS)")
    } else {
        title(main="Dissimilarities within demes\nGray means a has a single individual sampled from")
    }
}
## By default, all figures are saved as bitmap png images. However,
## it is straightforward to use another format (Here the alternative is pdf)
save.graphics <- function(plotpath,plot.height,plot.width,res=600,out.png=TRUE) {
    if (out.png) {
        bitmap(paste(plotpath,'%02d.png',sep=''),type='png16m',res=res,
               height=plot.height,width=plot.width,units='in')
    } else {
        pdf(paste(plotpath,'%02d.pdf',sep=''),
              height=plot.height,width=plot.width,onefile=FALSE)
    }
}
####################################################################################
## eems.plots takes three required arguments:
##   mcmcpath: one or several output directories (for the same dataset)
##   plotpath: filename of figures to generate
##     The following figures are created by default:
##     * plotpath-mrates01.png: effective migration rates
##     * plotpath-qrates01.png: effective diversity rates
##     * plotpath-rist01/02.png: fitted vs observed distances
##     * plotpath-pilogl01.png: trace of posterior probability
##   longlat (TRUE or FALSE): are the coordinates ordered longitude/latitude or not?
##
## eems.plots take a variety of optional arguments:
##   plot.width and plot.height: width and height of the graphics region (in inches)
##   plot.voronoi (TRUE or FALSE): plot a few posterior Voronoi diagrams?
##   add.map: add 'worldHires' map (using the rworldmap package)?
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
##   projection: specify the projection (the default is 'mercator')
####################################################################################
eems.plots <- function(mcmcpath,plotpath,longlat,plot.width=0,plot.height=0,out.png=TRUE,
                       add.map=FALSE,add.grid=TRUE,add.outline=TRUE,add.samples=TRUE,plot.voronoi=FALSE,
                       col.map=default.map.color(),col.grid=default.grid.color(),col.outline="white",col.samples="black",
                       lwd.map=1,lwd.grid=1,lwd.outline=2,pch.samples=19,cex.samples=1,max.cex.samples=2,
                       projection=NULL) {

    proj.in <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +no_defs"
    proj.out <- "+proj=merc"
    if (!is.null(projection)) { proj.out <- projection; }
    if (add.map) { require(rworldmap); require(rworldxtra); }
    if (max.cex.samples<0) { max.cex.samples <- 0; }
    
    plot.params <- list(add.map=add.map,add.grid=add.grid,add.outline=add.outline,add.samples=add.samples,
                        col.map=col.map,col.grid=col.grid,col.outline=col.outline,col.samples=col.samples,
                        lwd.map=lwd.map,lwd.grid=lwd.grid,lwd.outline=lwd.outline,pch.samples=pch.samples,
                        cex.samples=cex.samples,max.cex.samples=max.cex.samples,
                        proj.in=proj.in,proj.out=proj.out)

    mcmcpath1 <- character()
    for (path in mcmcpath) {
        if (file.exists(paste(path,'/ipmap.txt',sep=''))&&
            file.exists(paste(path,'/demes.txt',sep=''))&&
            file.exists(paste(path,'/edges.txt',sep=''))) {
            mcmcpath1 <- c(mcmcpath1,path)
        } else {
            print('The following EEMS output not found:')
            print(path)
        }
    }
    mcmcpath <- mcmcpath1
    nsimnos <- length(mcmcpath)
    if (nsimnos==0) { return(0) }

    print('Processing the following EEMS output:')
    print(mcmcpath)

    dimns <- read.dimns(mcmcpath[1],longlat)
    if ((plot.height<=0)||(plot.width<=0)) {
        plot.height <- min(dimns$yspan,12)
        plot.width <- min(dimns$xspan,12)
    }

    ## Plot filled contour of estimated effective migration rates
    save.graphics(paste(plotpath,'-mrates',sep=''),out.png=out.png,
                  plot.height=plot.height,plot.width=plot.width)
    average.eems.contours(mcmcpath,dimns,longlat,plot.params,is.mrates=TRUE)
    dev.off( )

    ## Plot fillet contour of estimated effective diversity rates
    save.graphics(paste(plotpath,'-qrates',sep=''),out.png=out.png,
                  plot.height=plot.height,plot.width=plot.width)
    average.eems.contours(mcmcpath,dimns,longlat,plot.params,is.mrates=FALSE)
    dev.off( )
    
    ## Plot Voronoi tessellations drawn from the posterior distributions on
    ## the migration and diversity rate parameters
    if (plot.voronoi) {
        save.graphics(paste(plotpath,'-mvoronoi',sep=''),out.png=out.png,
                      plot.height=plot.height,plot.width=plot.width)
        voronoi.diagram(mcmcpath,dimns,longlat,plot.params,is.mrates=TRUE)
        dev.off( )    
        save.graphics(paste(plotpath,'-qvoronoi',sep=''),out.png=out.png,
                      plot.height=plot.height,plot.width=plot.width)
        voronoi.diagram(mcmcpath,dimns,longlat,plot.params,is.mrates=FALSE)
        dev.off( )
    }

    plot.height <- 5
    plot.width <- 6

    ## Plot scatter plots of observed vs fitted genetic differences
    save.graphics(paste(plotpath,'-rdist',sep=''),out.png=out.png,
                  plot.height=plot.height,plot.width=plot.width)
    dist.scatterplot(mcmcpath)
    dev.off( )
    ## Plot trace plot of posterior probability to check convergence
    save.graphics(paste(plotpath,'-pilogl',sep=''),out.png=out.png,
                  plot.height=plot.height,plot.width=plot.width)
    plot.logposterior(mcmcpath)
    dev.off( )

}
