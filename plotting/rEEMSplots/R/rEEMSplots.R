
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

mrates.levels <- function(Zvals,standardize=TRUE) {
    minZ <- min(Zvals)
    maxZ <- max(Zvals)
    eems.colors <- default.eems.colors()
    numlevels <- length(eems.colors)
    maxZ2 <- max(maxZ,-minZ)
    maxZ2 <- ceiling(maxZ2*1000)/1000
    minZ <- -maxZ2
    maxZ <- +maxZ2
    eems.levels <- numeric()
    if (standardize) {
        ## First check whether the Z values are in the range (-2.5,+2.5), since
        ## it might be nice to have the migration rates always plotted on the same scale
        if ((minZ > -2.5)&&(maxZ < +2.5)) {
            minZ <- -2.5
            maxZ <- +2.5
            step <- (maxZ-minZ)/numlevels
            eems.levels <- seq(from=minZ,to=maxZ,by=step)
        } else {
            eems.levels <- seq(from=minZ,to=maxZ,length=numlevels+1)
        }
    } else {
        ## Otherwise, there is no default range for the diversity rates
        ## unless there is almost no variation, i.e., (maxZ-minZ)<0.001
        if ((maxZ - minZ)>0.1) {
            eems.levels <- seq(from=minZ,to=maxZ,length=numlevels+1)
        } else {
            eems.levels <- seq(from=minZ-0.1/2,to=maxZ+0.1/2,length=numlevels+1)
        }
    }
    return(eems.levels)
}
qrates.levels <- function(Zvals,standardize=TRUE) {
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
## The distance metric used by EEMS (Euclidean distance by default)
which.dist.metric <- function(mcmcpath) {
  dist.metric <- "euclidean"
  Lines = readLines(paste(mcmcpath,"/eemsrun.txt",sep=""))
  nLines = length(Lines)
  for (i in seq(nLines)) {
    s = gsub("\\s","",Lines[i])
    x = strsplit(s,"distance=")[[1]]
    if (length(x)==2) {
      dist.metric = x[2]
    }
  }
  return (dist.metric)
}
read.dimns <- function(mcmcpath,longlat,nxmrks=NULL,nymrks=NULL) {
    outer <- scan(paste(mcmcpath,'/outer.txt',sep=''),what=numeric(),quiet=TRUE)
    outer <- matrix(outer,ncol=2,byrow=TRUE)
    if (!longlat) {
        outer <- outer[,c(2,1)]
    }
    xlim <- range(outer[,1])
    ylim <- range(outer[,2])
    aspect <- abs((diff(ylim)/diff(xlim))/cos((mean(ylim) * pi)/180))
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
                marks=marks,nmrks=c(nxmrks,nymrks),aspect=aspect,
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
standardize.rates <- function(mcmcpath,dimns,longlat,is.mrates,standardize) {
    voronoi <- read.voronoi(mcmcpath,longlat,is.mrates)
    rates <- voronoi$rates
    tiles <- voronoi$tiles
    xseed <- voronoi$xseed
    yseed <- voronoi$yseed
    Zvals <- matrix(0,dimns$nxmrks,dimns$nymrks)
    niter <- length(tiles)
    if (standardize) {
        Zvals <- .Call("rEEMSplots__rcppstandardize_rates", PACKAGE = "rEEMSplots",
                       tiles, rates, xseed, yseed, dimns$marks, dimns$nmrks, dimns$dist.metric)
    } else {
        Zvals <- .Call("rEEMSplots__rcppnotstandardize_rates", PACKAGE = "rEEMSplots",
                       tiles, rates, xseed, yseed, dimns$marks, dimns$nmrks, dimns$dist.metric)
    }
    return(list(Zvals=Zvals,niter=niter))
}
filled.countour.axes <- function(mcmcpath,longlat,plot.params) {
    if (!is.null(plot.params$proj.in) && !is.null(plot.params$proj.out)) {
        filled.countour.axes.proj.known(mcmcpath,longlat,plot.params)
    } else {
        filled.countour.axes.proj.unknown(mcmcpath,longlat,plot.params)
    }
}
filled.countour.axes.proj.unknown <- function(mcmcpath,longlat,plot.params) {
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
    boundary <- SpatialPolygons(list(Polygons(list(Polygon(outer, hole = FALSE)),"1")))
    boundary <- SpatialPolygonsDataFrame(boundary,data.frame(id="1"))
    exterior <- gDifference(gEnvelope(boundary),boundary)
    if (!is.null(exterior)) {
        plot(exterior,col="white",border="white",add=TRUE)
    }
    if (plot.params$add.grid) {
        segments <- list()
        for (e in 1:nrow(edges)) {
            segments[[e]] <- Line(cbind(demes[edges[e,],1],demes[edges[e,],2]))
        }
        segments <- SpatialLines(list(Lines(segments,ID="a")))
        lines(segments,col=plot.params$col.grid,lwd=plot.params$lwd.grid)
    }
    if (plot.params$add.outline) {
        plot(boundary,col=NA,border=plot.params$col.outline,lwd=plot.params$lwd.outline,add=TRUE)
    }
    if (plot.params$add.samples) {
        samples <- SpatialPoints(cbind(demes[alpha,1],demes[alpha,2]))
        points(samples,col=plot.params$col.samples,pch=plot.params$pch.samples,
               cex=plot.params$cex.samples+plot.params$max.cex.samples*sizes/max(sizes))
    }
}
filled.countour.axes.proj.known <- function(mcmcpath,longlat,plot.params) {
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
one.eems.contour <- function(mcmcpath,dimns,Zmean,longlat,plot.params,is.mrates) {
    eems.colors <- default.eems.colors( )
    if (is.mrates) {
        eems.levels <- mrates.levels(Zmean,plot.params$standardize)
        main.title <- "Effective migration rates m : posterior mean"
        var.title <- "Effective migration rates m : posterior variance"
        key.title <- "m"
    } else {
        eems.levels <- qrates.levels(Zmean,plot.params$standardize)
        main.title <- "Effective diversity rates q : posterior mean"
        var.title <- "Effective diversity rates q : posterior variance"
        key.title <- "q"
    }
    rr <- flip(raster(t(Zmean),
                      xmn=dimns$xlim[1],xmx=dimns$xlim[2],
                      ymn=dimns$ylim[1],ymx=dimns$ylim[2]),direction='y')
    if (!is.null(plot.params$proj.in) && !is.null(plot.params$proj.out)) {
        projection(rr) <- CRS(plot.params$proj.in)
        rr <- projectRaster(rr,crs=CRS(plot.params$proj.out))
    }
    myfilledContour(rr,
                    main=main.title,font.main=1,asp=1,
                    col=eems.colors,levels=eems.levels,frame.plot=FALSE,
                    key.axes = axis(4,tick=FALSE,hadj=1,line=3,cex.axis=1.5),
                    key.title = mtext(key.title,side=3,cex=1.5,font=1),
                    plot.axes = filled.countour.axes(mcmcpath,longlat,plot.params))
}
average.eems.contours <- function(mcmcpath,dimns,longlat,plot.params,is.mrates) {
    mcmcpath1 <- character( )
    if (is.mrates) {
        writeLines('Plotting effective migration rates m : posterior mean and variance')
        for (path in mcmcpath) {
            Files <- paste(path,c('/mcmcmtiles.txt','/mcmcmrates.txt',
                                  '/mcmcxcoord.txt','/mcmcycoord.txt'),sep='')
            if (min(file.exists(Files))) { mcmcpath1 <- c(mcmcpath1,path) }
        }
    } else {
        writeLines('Plotting effective diversity rates q : posterior mean and variance')
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
        writeLines(path)
        rslt <- standardize.rates(path,dimns,longlat,is.mrates,plot.params$standardize)
        niter <- niter + rslt$niter
        Zmean <- Zmean + rslt$Zvals
    }
    Zmean <- Zmean/niter
    one.eems.contour(mcmcpath[1],dimns,Zmean,longlat,plot.params,is.mrates)
}
## If there are multiple runs, pick the first one and create at most max.niter Voronoi diagrams
## This function is mainly for testing purposes, so no need to create a plot for each iteration
voronoi.diagram <- function(mcmcpath,dimns,longlat,plot.params,is.mrates,max.niter=10) {
    mcmcpath <- mcmcpath[1]
    writeLines('Plotting Voronoi tessellation of estimated effective rates')
    writeLines(mcmcpath)
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
    writeLines('Plotting log posterior')
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
        path <- mcmcpath[i]; writeLines(path)
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
    writeLines('Plotting average dissimilarities within and between demes')
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
            path <- mcmcpath[simno]; writeLines(path)
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
        writeLines(path)
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
        writeLines('Need at least two observed demes to plot pairwise differences')
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
         xlab=expression(paste("Fitted dissimilarity between demes  ",hat(D)[ab]," - (",hat(D)[aa],"+",hat(D)[bb],")/2",sep="")),
         ylab=expression(paste("Observed dissimilarity between demes  ",D[ab]," - (",D[aa],"+",D[bb],")/2",sep="")))
    if (remove.singletons) {
        title(main="Dissimilarities between demes\nSingleton demes, if any, excluded from plot (not from EEMS)")
    } else {
        title(main="Dissimilarities between demes\nGray means a or b has a single individual sampled from")
    }
    ypts <- Wobs
    xpts <- What
    plot(xpts,ypts,col=c("black","gray60")[1+1*(Sizes==1)],
         xlab=expression(paste("Fitted dissimilarity within demes  ",hat(D)[aa],sep="")),
         ylab=expression(paste("Observed dissimilarity within demes  ",D[aa],sep="")))
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

###### Redefine the filled contour function (I don't like how the legend looks) ######

## Just copied over the 'filledContour' function in the 'raster' package
## For now, I have only changed how the legend/bar is plotted,
## the change itself is in the 'filled.contour' function
myfilledContour <- function (x, y = 1, maxpixels = 1e+05, ...) {
    if (nlayers(x) > 1) {
        y <- min(max(1, y), nlayers(x))
        x <- raster(x, y)
    }
    x <- sampleRegular(x, maxpixels, asRaster = TRUE, useGDAL = TRUE)
    X <- xFromCol(x, 1:ncol(x))
    Y <- yFromRow(x, nrow(x):1)
    Z <- t(matrix(getValues(x), ncol = x@ncols, byrow = TRUE)[nrow(x):1,])
    myfilled.contour(x = X, y = Y, z = Z, ...)
}
## I have only changed the legend/bar to remove the black border
myfilled.contour <- function (x = seq(0, 1, length.out = nrow(z)), y = seq(0, 1,length.out = ncol(z)), z,
                              xlim = range(x, finite = TRUE),
                              ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE),
                              levels = pretty(zlim, nlevels), nlevels = 20, color.palette = cm.colors,
                              col = color.palette(length(levels) - 1), plot.title, plot.axes,
                              key.title, key.axes, asp = NA, xaxs = "i", yaxs = "i", las = 1,
                              axes = TRUE, frame.plot = axes, ...) {
    if (missing(z)) {
        if (!missing(x)) {
            if (is.list(x)) {
                z <- x$z
                y <- x$y
                x <- x$x
            }
            else {
                z <- x
                x <- seq.int(0, 1, length.out = nrow(z))
            }
        }
        else stop("no 'z' matrix specified")
    }
    else if (is.list(x)) {
        y <- x$y
        x <- x$x
    }
    if (any(diff(x) <= 0) || any(diff(y) <= 0))
        stop("increasing 'x' and 'y' values expected")
    mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
    on.exit(par(par.orig))
    w <- (3 + mar.orig[2L]) * par("csi") * 2.54
    layout(matrix(c(2, 1), ncol = 2L), widths = c(1, lcm(w)))
    par(las = las)
    mar <- mar.orig
    mar[4L] <- mar[2L]
    mar[2L] <- 1
    par(mar = mar)
    plot.new()
    ## Change here: added 'border = NA' and removed 'box()'
    plot.window(xlim = c(0, 1), ylim = range(levels), xaxs = "i",yaxs = "i")
    rect(0, levels[-length(levels)], 1, levels[-1L], col = col, border = NA)
    if (missing(key.axes)) {
        if (axes)
            axis(4)
    }
    else key.axes
    ## box()
    if (!missing(key.title))
        key.title
    mar <- mar.orig
    mar[4L] <- 1
    par(mar = mar)
    plot.new()
    plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp)
    .filled.contour(x, y, z, levels, col)
    if (missing(plot.axes)) {
        if (axes) {
            title(main = "", xlab = "", ylab = "")
            Axis(x, side = 1)
            Axis(y, side = 2)
        }
    }
    else plot.axes
    if (frame.plot)
        box()
    if (missing(plot.title))
        title(...)
    else plot.title
    invisible()
}
load.required.package <- function(package) {
    available = suppressMessages(suppressWarnings(require(package,
        quietly = TRUE, character.only = TRUE, warn.conflicts = FALSE)))
    if (!available) {
        stop(paste("The option 'plot.voronoi=TRUE' requires the '",package,"' package",sep=""))
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
##   projection.in: specify the input projection 
##   projection.out: specify the output projection 
####################################################################################
eems.plots <- function(mcmcpath,plotpath,longlat,plot.width=0,plot.height=0,out.png=TRUE,
                       add.map=FALSE,add.grid=TRUE,add.outline=TRUE,add.samples=TRUE,plot.voronoi=FALSE,
                       col.map=default.map.color(),col.grid=default.grid.color(),col.outline="white",col.samples="black",
                       lwd.map=1,lwd.grid=1,lwd.outline=2,pch.samples=19,cex.samples=1,max.cex.samples=2,
                       projection.in=NULL,projection.out=NULL,standardize=TRUE) {

    if (!is.null(projection.in) && is.null(projection.out)) {
        projection.out = projection.in
    }
    writeLines(paste("Input projection: ",projection.in,"\n",
                     "Output projection: ",projection.out,sep=""))
    if (is.null(projection.in) && !is.null(projection.out)) {
        stop("Specify the input projection, projection.in")
    }
    if (add.map && (is.null(projection.in) || is.null(projection.out))) {
        stop(paste("To add a geographical map, specify the input and output projections.\n",
                   "For example, if the coordinates are longitude and latitude,\n",
                   "you can use the PROJ.4 string '+proj=loglat +datum=WGS84'",sep=""))
    }

    if (!is.null(projection.in) ||
        !is.null(projection.out)) {
        load.required.package('rgdal')
    }
    if (add.map) {
        load.required.package('rworldmap')
        load.required.package('rworldxtra')
    }
    if (plot.voronoi) {
        load.required.package('deldir')
    }
    
    plot.params <- list(add.map=add.map,add.grid=add.grid,add.outline=add.outline,add.samples=add.samples,
                        col.map=col.map,col.grid=col.grid,col.outline=col.outline,col.samples=col.samples,
                        lwd.map=lwd.map,lwd.grid=lwd.grid,lwd.outline=lwd.outline,pch.samples=pch.samples,
                        cex.samples=cex.samples,max.cex.samples=max.cex.samples,
                        proj.in=projection.in,proj.out=projection.out,standardize=standardize)

    mcmcpath1 <- character()
    for (path in mcmcpath) {
        if (file.exists(paste(path,'/ipmap.txt',sep=''))&&
            file.exists(paste(path,'/demes.txt',sep=''))&&
            file.exists(paste(path,'/edges.txt',sep=''))) {
            mcmcpath1 <- c(mcmcpath1,path)
        } else {
            writeLines('The following EEMS output not found:')
            writeLines(path)
        }
    }
    mcmcpath <- mcmcpath1
    nsimnos <- length(mcmcpath)
    if (nsimnos==0) { return(0) }

    writeLines('Processing the following EEMS output:')
    writeLines(mcmcpath)

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
        par(las=1,font.main=1)
        voronoi.diagram(mcmcpath,dimns,longlat,plot.params,is.mrates=TRUE)
        dev.off( )    
        save.graphics(paste(plotpath,'-qvoronoi',sep=''),out.png=out.png,
                      plot.height=plot.height,plot.width=plot.width)
        par(las=1,font.main=1)
        voronoi.diagram(mcmcpath,dimns,longlat,plot.params,is.mrates=FALSE)
        dev.off( )
    }

    plot.height <- 5
    plot.width <- 6

    ## Plot scatter plots of observed vs fitted genetic differences
    save.graphics(paste(plotpath,'-rdist',sep=''),out.png=out.png,
                  plot.height=plot.height,plot.width=plot.width)
    par(las=1,font.main=1)
    dist.scatterplot(mcmcpath)
    dev.off( )
    ## Plot trace plot of posterior probability to check convergence
    save.graphics(paste(plotpath,'-pilogl',sep=''),out.png=out.png,
                  plot.height=plot.height,plot.width=plot.width)
    par(las=0,font.main=1)
    plot.logposterior(mcmcpath)
    dev.off( )
    
}
