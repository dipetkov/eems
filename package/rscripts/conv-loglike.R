

rm(list=ls( ))

datapath <- '~/package/data'
plotpath <- '.'

source('src/myfilled.contour.R')
source('src/mydichromat.R')
library(fields) ## Computes Euclidean distances ##


q1 <- function(x) {
  return (quantile(x,0.1))
}
q5 <- function(x) {
  return (quantile(x,0.5))
}
q9 <- function(x) {
  return (quantile(x,0.9))
}

plotlikeli2 <- function(datapath,dataset,config,simnos) {

    filebase <- paste(datapath,'/',dataset,'-simno',simnos[1],'-',config,sep='')
    piloglike <- read.table(paste(filebase,'.mcmcpilogl',sep=''),header=FALSE)
    nsimno <- length(simnos)
    niters <- nrow(piloglike)
    posterior <- matrix(0,niters,nsimno)
    
    for (k in 1:nsimno) {
        filebase <- paste(datapath,'/',dataset,'-simno',simnos[k],'-',config,sep='')
        piloglike <- read.table(paste(filebase,'.mcmcpilogl',sep=''),header=FALSE)
        posterior[,k] <- piloglike[,1] + piloglike[,2]
    }
    
    matplot(posterior,type="l",col=varcolors,lwd=3,lty=1,
            xlab="iteration (after burn-in)",ylab="(log) posterior",
            main=paste(dataset,"-",config,sep=""),
            axes=FALSE,xlim=c(1,1.06*niters),bty="n")
    legend(1.01*niters,max(posterior),legend=simnos,col=varcolors,lwd=5,bty="n")
}


dataset <- 'hap-tribars-nIndiv300-nSites3000-gridSize12x8'
xPop <- 12
yPop <- 8
simnos <- 1:6


pdf(height=7,width=14,
    file=paste(plotpath,'/',dataset,
        '-g',xPop,'x',yPop,
        '-loglike.pdf',sep=''))
config = paste('g',xPop,'x',yPop,sep='')
plotlikeli2(datapath,dataset,config,simnos)
dev.off( )
