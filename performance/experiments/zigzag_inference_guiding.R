library(ggplot2)

source("stat-stepribbon.R")

calculate.quartiles <- function( frame, field="Ne", min=1000, max=100000 ) {

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
    return(new.df)
}

plot.smcsmc <- function( data, g=30 ) {

    ##write.table(data, file="zigzag.dta", sep="\t")
    ##data <- read.table(file="zigzag.dta", header=TRUE)

    t0 <- 500
    t1 <- 2000000
    emstep <- max( data$iter )

    plot.data <- subset(data, iter==emstep & type=="Coal")
    quartiles <- calculate.quartiles( plot.data, min=100, max=1e5 )
    quartiles$Population <- paste( "Pop", quartiles$frm + 1, sep="" )
    ## make first epoch start at t0
    quartiles[ quartiles$start * g < t0, ]$start <- t0/g
    ## make last epoch end at t1
    endrows = quartiles[ quartiles$epoch == max(quartiles$epoch), ]
    endrows$start = t1 / g
    quartiles = rbind( quartiles, endrows )
    
    p <- ggplot(data=quartiles, aes(x=start*g))
    p <- p + scale_x_log10(limits=c(t0,t1),breaks=c(1000,10000,100000,1000000))
    p <- p + scale_y_log10(limits=c(1000,100000),breaks=c(2000,5000,10000,20000,50000))
    p <- p + geom_ribbon(data=subset(quartiles, Population=="Pop1"), aes(ymin=Q0, ymax=Q4),
                         stat="stepribbon", fill="blue", alpha=.15)
    p <- p + geom_ribbon(data=subset(quartiles, Population=="Pop1"), aes(ymin=Q1, ymax=Q3),
                         stat="stepribbon", fill="blue", alpha=.25)
    p <- p + geom_step(data=subset(quartiles, Population=="Pop1"), aes( y=Q2, x=start*g, colour=Population), lwd=1 )
    p <- p + xlab("years ago") + ylab("Ne") + theme(text = element_text(size=16), legend.position="top")
    #p <- p + geom_step(data=unidirmigr_true_parameters,aes(unidirmigr_times_in_years, unidirmigr_true_Ne_pop0, color="Truth for both"),lwd=1.1)
    p <- p + scale_colour_manual("",
                                 breaks = c("Pop1","Pop2","Truth for both"),
                                 values = c("Pop1"="blue","Pop2"="red","Truth for both"="black"))
    p
    ggsave(paste("zigzag_",emstep,"EMsteps.png",sep=""),
           width = 12, height = 8, units = "cm")
    
}

