

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
mycolors <- c("#994000","#CC5800","#FF8F33","#FFAD66","#FFCA99","#FFE6CC",
              "#FFFFFF",
              "#CCFDFF","#99F8FF","#66F0FF","#33E4FF","#00AACC","#007A99")


numlevels <- length(mycolors)
minZ <- -2.5
maxZ <- +2.5
step <- (maxZ-minZ)/numlevels
the.pts <- c(  2.1 ,  1.4 ,  0.7 ,    0 , -0.7 , -1.4 , -2.1 )
the.lab <- c(' 2.1',' 1.4',' 0.7','   0','-0.7','-1.4','-2.1')
mylevels <- seq(from=minZ,to=maxZ,by=step)


## Additional color schemes can be designed based on the 'dichromat' package
extra.color.schemes <- function(SchemeName='DarkOrangetoBlue2') {
    library(dichromat)
    ## 'DarkOrangetoBlue2' is the default color scheme
    if (SchemeName == 'BrowntoBlue') {
        colscheme <- colorschemes$BrowntoBlue.12
    } else if (SchemeName == 'BluetoBrown') {
        colscheme <- rev(colorschemes$BrowntoBlue.12)
    } else if (SchemeName == 'BluetoDarkOrange') {
        colscheme <- colorschemes$BluetoDarkOrange.12
    } else if (SchemeName == 'BluetoDarkOrange2') {
        Blues <- colorschemes$BrowntoBlue.12[12:7]
        DarkOranges <- colorschemes$BluetoDarkOrange.12[7:12]
        colscheme <- c(Blues,DarkOranges)
    } else if (SchemeName == 'DarkOrangetoBlue') {
        colscheme <- rev(colorschemes$BluetoDarkOrange.12)
    } else if (SchemeName == 'DarkOrangetoBlue2') {
        Blues <- colorschemes$BrowntoBlue.12[12:7]
        DarkOranges <- colorschemes$BluetoDarkOrange.12[7:12]
        colscheme <- rev(c(Blues,DarkOranges))
    } else if (SchemeName == 'DarkRedtoBlue') {
        colscheme <- colorschemes$DarkRedtoBlue.12
    } else if (SchemeName == 'BluetoDarkRed') {
        colscheme <- rev(colorschemes$DarkRedtoBlue.12)
    }
    mycolors <- c(colscheme[1:6],"#FFFFFF",colscheme[7:12])
    return(mycolors)
}
