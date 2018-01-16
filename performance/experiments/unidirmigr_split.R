library(ggplot2)
library(grid)
library(gridExtra)

source("stat-stepribbon.R")
source("plot-utils.R")

plot.ribbon <- function( p, data, col="blue", g=30, lwd=1 ) {
    p <- p + geom_ribbon(data=data, aes(ymin=Q0, ymax=Q4),
                     stat="stepribbon", fill=col, alpha=.15)
    p <- p + geom_ribbon(data=data, aes(ymin=Q1, ymax=Q3),
                         stat="stepribbon", fill=col, alpha=.25)
    p <- p + geom_step(data=data, aes( y=Q2, x=start*g, colour=Population), lwd=lwd )
    return (p)
}

plot.truth <- function( p, x, y, pop, g=30, lwd=1,lty=1 ) {
    data <- data.frame(start=x,
                       end=c(x[-1],tail(x,1)),
                       Q2=y,
                       Population=pop)
    p <- p + geom_step(data=data, aes( y=Q2, x=start, colour=Population), lwd=lwd, lty=lty )
    return (p)
}


plot.model <- function( data, model, g=30 ) {

    migr_time <- 50e3
    split_time <- 100e3
    t0 <- 500
    t1 <- 1500000
    t1m = split_time
    N0 <- 10000
    emstep <- max( data$iter )

    plot.data <- subset(data, iter==emstep & type=="Coal" & str_parameter==model)
    quartiles <- calculate.quartiles( plot.data, min=100, max=1e5, mint=t0/g, maxt=t1/g )
    quartiles$Population <- paste( "Pop", quartiles$frm + 1, sep="" )

    data.pop1 = subset(quartiles, Population=="Pop1")
    data.pop2 = subset(quartiles, Population=="Pop2" & epoch < 12)
    data.pop2 = rbind( data.pop2, data.pop1[ data.pop1$epoch == 12, ])
    data.pop2$Population <- "Pop2"

    p <- ggplot(data=quartiles, aes(x=start*g))
    p <- p + scale_x_log10(limits=c(t0,t1),breaks=c(1000,10000,100000,1000000), labels=NULL)
    labels <- NULL
    if (model == 0) labels <- waiver()
    p <- p + scale_y_log10(limits=c(400,100000),breaks=c(500,1000,2000,5000,10000,20000,50000), labels=labels)
    p <- plot.ribbon( p, data.pop2, "darkcyan" )
    p <- plot.truth( p, c(500, migr_time, split_time), c(0.5*N0, 0.05*N0, 1*N0), "Pop2", lty="11", lwd=0.5 )
    p <- plot.ribbon( p, data.pop1, "blue" )
    p <- plot.truth( p, c(500, migr_time, split_time, t1), c(1*N0, 0.95*N0, 1*N0, 1*N0), "Pop1", lty="11", lwd=0.5 )

    if (model == 0) {
        p <- p + ylab("Ne")
    } else {
        p <- p + theme(axis.title.y = element_blank() )
    }
    p <- p + theme(text = element_text(size=16), legend.position="top")
    p <- p + theme(axis.title.x = element_blank())
    p <- p + scale_colour_manual("",
                                 breaks = c("Pop1","Pop2"),
                                 values = c("Pop1"="blue","Pop2"="darkcyan"))

    "migration"
    plot.data <- subset(data, iter==emstep & type=="Migr" & str_parameter==model)
    quartiles <- calculate.quartiles( plot.data, field="rate", mint=t0/g, maxt=t1m/g, miny=0, maxy=10, scale=4*N0 )
    quartiles$Population <- paste( "Pop", quartiles$frm + 1, sep="" )

    data.pop1 = subset(quartiles, Population=="Pop1")
    data.pop2 = subset(quartiles, Population=="Pop2" & epoch < 12)

    pm <- ggplot(data=quartiles, aes(x=start*g))
    pm <- pm + scale_x_log10(limits=c(t0,t1),breaks=c(1000,10000,100000,1000000))
    labels <- NULL
    if (model == 0) labels <- waiver()
    pm <- pm + scale_y_continuous(limits=c(0,5), labels=labels)
    pm <- plot.ribbon( pm, data.pop1, "blue" )
    pm <- plot.ribbon( pm, data.pop2, "darkcyan" )

    pm <- pm + xlab("years ago")
    if (model == 0) {
        pm <- pm + ylab("Migration rate") + theme(text = element_text(size=16), legend.position="none")
    } else {
        pm <- pm + theme(axis.title.y = element_blank() ) + theme(legend.position="none")
    }
    pm <- pm + scale_colour_manual("",
                                 breaks = c("Pop1","Pop2"),
                                 values = c("Pop1"="blue","Pop2"="darkcyan"))


    grob <- rbind(ggplotGrob(p), ggplotGrob(pm), size="last")
    return(grob)
}


plot.smcsmc <- function( data, g=30 ) {

    write.table(data, file="udm_split.dta", sep="\t")
    data <- read.table(file="udm_split.dta", header=TRUE)

    gr0 <- plot.model( data, 0 )
    gr1 <- plot.model( data, 1 )
    gr2 <- plot.model( data, 2 )
    gr3 <- plot.model( data, 3 )

    grid.draw( cbind( gr0,gr1,gr2,gr3,size="last" ))

    ## """this doesn't work..."""
    emstep <- max( data$iter )
    ggsave(paste("udm_split_",emstep,"EMsteps.png",sep=""),
           width = 12, height = 8, units = "cm")


}

