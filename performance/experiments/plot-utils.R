calculate.quartiles <- function( frame, field="Ne", min=1000, max=100000, mint=30, maxt=10000 ) {

    options( stringsAsFactors=FALSE ) # seems not possible to give this as option to rbind
    new.df <- data.frame()
    factor.idx <- sapply( frame, is.factor )
    frame[ factor.idx ] <- lapply( frame[ factor.idx ], as.character )
    keys <- unique(data.frame(frame$type, frame$frm, frame$epoch))
    for (i in 1:nrow(keys)) {
        data <- subset(frame, type==keys[i,1] & frm==keys[i,2] & epoch==keys[i,3])
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
