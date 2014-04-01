

library(dichromat) ## Defines color schemes appropriate for color-blindness ##


## Choose the Blue to Dark Orange color scheme (with 12 levels)
## Make the halfway level "white"
numlevels <- 13
colscheme <- colorschemes$BluetoDarkOrange.12
mycolors <- c(colscheme[1:(numlevels/2)],"white",colscheme[(1+numlevels/2):numlevels])
minZ <- -2.5
maxZ <- +2.5
step <- (maxZ-minZ)/numlevels
the.pts <- c(2.1,1.4,0.7,0,-0.7,-1.4,-2.1)
the.lab <- c(' 2.1',' 1.4',' 0.7','   0','-0.7','-1.4','-2.1')
mylevels <- seq(from=minZ,to=maxZ,by=step)

##############################################################

numcols <- 36
alpha <- 0.75
varcolors <- character(numcols)
varcolors[1]  <- rgb(0.835,0.804,0.584,alpha)
varcolors[2]  <- rgb(0.969,0.412,0.275,alpha)
varcolors[3]  <- rgb(1.000,0.784,0.133,alpha)
varcolors[4]  <- rgb(0.576,0.757,0.235,alpha)
varcolors[5]  <- rgb(0.220,0.447,0.447,alpha)
varcolors[6]  <- rgb(0.878,0.553,0.227,alpha)
varcolors[7]  <- rgb(0.357,0.302,0.702,alpha)
varcolors[8]  <- rgb(0.306,0.310,0.169,alpha)
varcolors[9]  <- rgb(0.808,0.176,0.475,alpha)
varcolors[10] <- rgb(0.259,0.702,0.725,alpha)
varcolors[11] <- rgb(0.333,0.200,0.000,alpha)
varcolors[12] <- rgb(0.439,0.545,0.855,alpha)
varcolors[13] <- rgb(0.784,0.808,0.275,alpha)
varcolors[14] <- rgb(0.627,0.000,0.086,alpha)
varcolors[15] <- rgb(0.204,0.275,0.561,alpha)
varcolors[16] <- rgb(0.800,0.800,0.800,alpha)
varcolors[17] <- rgb(0.576,0.149,0.502,alpha)
varcolors[18] <- rgb(0.400,0.400,0.400,alpha)
varcolors[19] <- rgb(0.310,0.710,0.824,alpha)
varcolors[20] <- rgb(0.306,0.306,0.306,alpha)
varcolors[21] <- rgb(0.643,0.949,0.984,alpha)
varcolors[22] <- rgb(0.125,0.255,0.737,alpha)
varcolors[23] <- rgb(0.937,0.745,0.306,alpha)
varcolors[24] <- rgb(0.169,0.525,0.310,alpha)
varcolors[25] <- rgb(0.796,0.267,0.255,alpha)
varcolors[26] <- rgb(0.518,0.510,0.208,alpha)
varcolors[27] <- rgb(0.200,0.200,0.200,alpha)
varcolors[28] <- rgb(0.025,0.373,0.208,alpha)
varcolors[29] <- rgb(0.576,0.329,0.224,alpha)
varcolors[30] <- rgb(0.200,0.600,0.600,alpha)
varcolors[31] <- rgb(0.922,0.937,0.196,alpha)
varcolors[32] <- rgb(0.259,0.702,0.725,alpha)
varcolors[33] <- rgb(0.576,0.329,0.000,alpha)
varcolors[34] <- rgb(0.800,0.600,0.200,alpha)
varcolors[35] <- rgb(0.973,0.706,0.078,alpha)
varcolors[36] <- rgb(0.600,0.298,0.098,alpha)

##############################################################

round2 <- function(x) { trunc(x+0.5) }
read.dimns <- function(datapath,dataset,xPop,yPop,nxmrks=NULL,nymrks=NULL) {
    
    dimns <- read.table(paste(datapath,'/',dataset,'.dimns',sep=''))
    
    xrange <- as.numeric(dimns[1,])
    yrange <- as.numeric(dimns[2,])
    xmidpt <- (xrange[2]+xrange[1])/2
    ymidpt <- (yrange[2]+yrange[1])/2
    xspan <- xrange[2]-xrange[1]
    yspan <- yrange[2]-yrange[1]

    if (is.null(nxmrks)&&is.null(nymrks)) {
        nxmrks <- 50
        nymrks <- round2(nxmrks*(yspan/xspan))
    }
    
    xmrks <- seq(xrange[1],xrange[2],length=nxmrks)
    ymrks <- seq(yrange[1],yrange[2],length=nymrks)

    return(list(xPop=xPop,xmrks=xmrks,xrange=xrange,xspan=xspan,
                yPop=yPop,ymrks=ymrks,yrange=yrange,yspan=yspan))
}    

