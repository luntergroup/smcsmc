calculate.quartiles <- function( frame, field="Ne", min=1000, max=100000, mint=30, maxt=10000 ) {

    options( stringsAsFactors=FALSE ) # seems not possible to give this as option to rbind
    new.df <- data.frame()
    factor.idx <- sapply( frame, is.factor )
    frame[ factor.idx ] <- lapply( frame[ factor.idx ], as.character )
    keys <- unique(data.frame(frame$type, frame$frm, frame$epoch, frame$aux_part_filt, frame$int_parameter, frame$np))
    for (i in 1:nrow(keys)) {
        data <- subset(frame, type==keys[i,1] & frm==keys[i,2] & epoch==keys[i,3] &
                              aux_part_filt==keys[i,4] & int_parameter==keys[i,5] & np==keys[i,6])
        if (field == "Ne") values <- data$ne else data <- values$rate
        q <- quantile(values)
        if (field == "Ne") data$ne <- q[3] else data$rate <- q[3]   # set to median
        new.row <- c( data[1,], pmin(max,pmax(min,q)))
        names(new.row) <- c( names(frame), "Q0","Q1","Q2","Q3","Q4" )
        new.df <- rbind( new.df, new.row)
    }
    ## make first epoch start at t0
    if (sum(new.df$start < mint) > 0)
        new.df[ new.df$start < mint, ]$start <- mint
    ## make last epoch end at t1
    endrows <- new.df[ new.df$epoch == max(new.df$epoch), ]
    endrows$start <- maxt
    new.df <- rbind( new.df, endrows )
    return(new.df)
}


## Multiple plot function
##
## ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
## - cols:   Number of columns in layout
## - layout: A matrix specifying the layout. If present, 'cols' is ignored.
##
## If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
## then plot 1 will go in the upper left, 2 will go in the upper right, and
## 3 will go all the way across the bottom.
##
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    library(grid)
    
    ## Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)

    numPlots = length(plots)
    
    ## If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
        ## Make the panel
        ## ncol: Number of columns of plots
        ## nrow: Number of rows needed, calculated from # of cols
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                         ncol = cols, nrow = ceiling(numPlots/cols))
    }

    if (numPlots==1) {
        print(plots[[1]])

    } else {
        ## Set up the page
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

        ## Make each plot, in the correct location
        for (i in 1:numPlots) {
            ## Get the i,j matrix positions of the regions that contain this subplot
            matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
            
            print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                            layout.pos.col = matchidx$col))
        }
    }
}

