
library(dichromat)

## Additional color schemes can be designed based on the 'dichromat' package
extra.color.schemes <- function(SchemeName='DarkOrangetoBlue2') {
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
    eems.colors <- c(colscheme[1:6],"#FFFFFF",colscheme[7:12])
    return(eems.colors)
}


## 'DarkOrangetoBlue2' is the default color scheme
eems.colors <- ('DarkOrangetoBlue2')
numlevels <- length(eems.colors)
