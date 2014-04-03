

## This is the Blue to Dark Orange scheme from the dichromat package,
## with "white" as the midpoint color
## The R package itself is based on a collection of color schemes for
## scientific data graphics:
## http://geog.uoregon.edu/datagraphics/color_scales.htm
## and 
## Light A and Bartlein PJ (2004). The End of the Rainbow? Color Schemes
## for Improved Data Graphics.
## EOS Transactions of the American Geophysical Union, 85(40), 385.

mycolors <- c("#1F8F99","#52C4CC","#99FAFF","#B2FCFF","#CCFEFF","#E6FFFF",
              "#FFFFFF","#FFE6CC","#FFCA99","#FFAD66","#FF8F33","#CC5800",
              "#994000")
numlevels <- length(mycolors)
minZ <- -2.5
maxZ <- +2.5
step <- (maxZ-minZ)/numlevels
the.pts <- c(  2.1 ,  1.4 ,  0.7 ,    0 , -0.7 , -1.4 , -2.1 )
the.lab <- c(' 2.1',' 1.4',' 0.7','   0','-0.7','-1.4','-2.1')
mylevels <- seq(from=minZ,to=maxZ,by=step)

## A collection of various colors

varcolors <- c("#D5CD95","#F76946","#FFC822","#93C13C","#387272","#E08D3A",
               "#5B4DB3","#4E4F2B","#CE2D79","#42B3B9","#553300","#708BDA",
               "#C8CE46","#A00016","#34468F","#CCCCCC","#932680","#666666",
               "#4FB5D2","#4E4E4E","#A4F2FB","#2041BC","#EFBE4E","#2B864F",
               "#CB4441","#848235","#333333","#065F35","#935439","#339999",
               "#EBEF32","#42B3B9","#935400","#CC9933","#F8B414","#994C19")
