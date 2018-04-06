library(ggplot2)
library(grid)
library(gridExtra)

source("stat-stepribbon.R")
source("plot-utils.R")

plot.model <- function( data, popidx, g=30, minne=100 ) {

    split_time <- 385e3
    splitepoch <- 21
    t0 <- 1000
    t1 <- 1500000
    t1m = split_time
    N0 <- 10000
    maxne <- 3e5
    min.migr.rate <- 2e-2
    ##max.migr.rate <- 5e2   ## for log scale
    max.migr.rate <- 20      ## for continuous scale
    emstep <- max( data$iter )

    model <- 2  ## for now

    plot.data <- subset(data, iter==emstep & type=="Coal" & int_parameter==popidx)
    plot.data.migration <- subset(data, iter==emstep & type=="Migr" & int_parameter==popidx)
    
    quartiles <- calculate.quartiles( plot.data, min=minne, max=maxne, mint=t0/g, maxt=t1/g )
    quartiles$Population <- paste( "Pop", quartiles$frm + 1, sep="" )

    data.pop1 = subset(quartiles, Population=="Pop1")
    data.pop2 = subset(quartiles, Population=="Pop2" & epoch < splitepoch)
    data.pop2 = rbind( data.pop2, data.pop1[ data.pop1$epoch == splitepoch, ])
    data.pop2$Population <- "Pop2"

    p <- ggplot(data=quartiles, aes(x=start*g))
    p <- p + scale_x_log10(limits=c(t0,t1),breaks=c(1000,10000,100000,1000000), labels=NULL)
    p <- p + annotation_logticks(sides="b")
    labels <- NULL
    if (model == 2) labels <- waiver()
    p <- p + scale_y_log10(limits=c(minne,200000),breaks=c(500,1000,2000,5000,10000,20000,50000), labels=labels)
    p <- plot.ribbon( p, data.pop2, "darkcyan" )
    p <- plot.ribbon( p, data.pop1, "blue" )

    if (model == 2) {
        p <- p + theme(text = element_text(size=16)) + ylab("Ne")
    } else {
        p <- p + theme(axis.title.y = element_blank() )
    }
    p <- p + theme(text = element_text(size=16), legend.position="top")
    p <- p + theme(axis.title.x = element_blank())
    p <- p + scale_colour_manual("",
                                 breaks = c("Pop1","Pop2"),
                                 values = c("Pop1"="blue","Pop2"="darkcyan"))

    "migration"
    plot.data <- plot.data.migration
    quartiles <- calculate.quartiles( plot.data, field="rate", mint=t0/g, maxt=t1m/g, miny=0, maxy=100, scale=4*N0, add=min.migr.rate )
    quartiles$Population <- paste( "Pop", quartiles$frm + 1, sep="" )

    data.pop1 = subset(quartiles, Population=="Pop1")
    data.pop2 = subset(quartiles, Population=="Pop2" & epoch < splitepoch)

    pm <- ggplot(data=quartiles, aes(x=start*g))
    if (model == 2) {
        pm <- pm + ylab("Migration rate") + theme(text = element_text(size=16), legend.position="none")
    } else {
        pm <- pm + theme(axis.title.y = element_blank() ) + theme(legend.position="none")
    }
    pm <- pm + scale_x_log10(limits=c(t0,t1),breaks=c(1000,10000,100000,1000000))
    pm <- pm + annotation_logticks(sides="b")
    labels <- NULL
    if (model == 0) labels <- waiver()
    pm <- pm + scale_y_continuous(limits=c(min.migr.rate,max.migr.rate), labels=labels)
    pm <- plot.ribbon( pm, data.pop1, "blue", maxy=max.migr.rate )
    pm <- plot.ribbon( pm, data.pop2, "darkcyan", maxy=max.migr.rate )

    pm <- pm + xlab("years ago")
    pm <- pm + scale_colour_manual("",
                                 breaks = c("Pop1","Pop2"),
                                 values = c("Pop1"="blue","Pop2"="darkcyan"))


    grob <- rbind(ggplotGrob(p), ggplotGrob(pm), size="last")
    return(grob)
}



plot.model.allemsteps <- function( data, popidx, g=30, minne=100 ) {

    split_time <- 385e3
    splitepoch <- 21
    t0 <- 1000
    t1 <- 1500000
    t1m = split_time
    N0 <- 10000
    maxne <- 3e5
    min.migr.rate <- 2e-2
    ##max.migr.rate <- 5e2   ## for log scale
    max.migr.rate <- 20      ## for continuous scale
    emstep <- max( data$iter )

    plot.data <- subset(data, type=="Coal" & int_parameter==popidx)
    plot.data.migration <- subset(data, type=="Migr" & int_parameter==popidx)
    
    quartiles <- calculate.quartiles( plot.data, min=minne, max=maxne, mint=t0/g, maxt=t1/g )
    quartiles$Population <- paste( "Pop", quartiles$frm + 1, sep="" )

    data.pop1 = subset(quartiles, Population=="Pop1")
    data.pop2 = subset(quartiles, Population=="Pop2" & epoch < splitepoch)
    data.pop2 = rbind( data.pop2, data.pop1[ data.pop1$epoch == splitepoch, ])
    data.pop2$Population <- "Pop2"

    p <- ggplot(data=quartiles, aes(x=start*g, colour=iter))
    p <- p + scale_x_log10(limits=c(t0,t1),breaks=c(1000,10000,100000,1000000))
    p <- p + annotation_logticks(sides="b")
    labels <- waiver()
    p <- p + scale_y_log10(limits=c(minne,200000),breaks=c(500,1000,2000,5000,10000,20000,50000), labels=labels)
    p <- p + geom_step(data=data.pop2, aes( y=Q2, x=start*g, colour=iter, group=iter), lwd=1 )
    p <- p + geom_step(data=data.pop1, aes( y=Q2, x=start*g, colour=iter, group=iter), lwd=1 )    

    p <- p + theme(text = element_text(size=16)) + ylab("Ne")
    p <- p + theme(text = element_text(size=16), legend.position="top")
    p <- p + theme(axis.title.x = element_blank())
    p <- p + facet_grid(Population ~ .)
    p <- p + scale_colour_gradientn(colours=rainbow(4))

    "migration"
    plot.data <- plot.data.migration
    quartiles <- calculate.quartiles( plot.data, field="rate", mint=t0/g, maxt=t1m/g, miny=0, maxy=100, scale=4*N0, add=min.migr.rate )
    quartiles$Population <- paste( "Pop", quartiles$frm + 1, sep="" )

    data.pop1 = subset(quartiles, Population=="Pop1")
    data.pop2 = subset(quartiles, Population=="Pop2" & epoch < splitepoch)

    pm <- ggplot(data=quartiles, aes(x=start*g, colour=iter))
    pm <- pm + ylab("Migration rate") + theme(text = element_text(size=16), legend.position="none")
    pm <- pm + scale_x_log10(limits=c(t0,t1),breaks=c(1000,10000,100000,1000000))
    pm <- pm + annotation_logticks(sides="b")
    labels <- waiver()
    pm <- pm + scale_y_continuous(limits=c(min.migr.rate,max.migr.rate), labels=labels)
    pm <- pm + geom_step(data=data.pop2, aes( y=Q2, x=start*g, colour=iter, group=iter), lwd=1 )
    pm <- pm + geom_step(data=data.pop1, aes( y=Q2, x=start*g, colour=iter, group=iter), lwd=1 )
    pm <- pm + facet_grid(Population ~ .)

    pm <- pm + xlab("years ago")
    pm <- pm + scale_colour_gradientn(colours=rainbow(4))

    grob <- rbind(ggplotGrob(p), ggplotGrob(pm), size="last")
    return(grob)
}






data <- load.from.out( "ceuchb4/result.out", int_parameter=1 )
data <- load.from.out( "ceugih4/result.out", data=data, int_parameter=2 )
data <- load.from.out( "ceumxl4/result.out", data=data, int_parameter=3 )
data <- load.from.out( "ceutsi4/result.out", data=data, int_parameter=4 )
data <- load.from.out( "ceuyri4/result.out", data=data, int_parameter=5 )
data$Population <- c("CHB","GIH","MXL","TSI","YRI")[ data$int_parameter ]

gr1 <- plot.model( data, 5, minne=1000 )
grid.draw( gr1 )

gr2 <- plot.model.allemsteps( data, 5, minne=1000 )
grid.draw( gr2 )

gr1 <- plot.model( data, 1, minne=1000 )
gr2 <- plot.model( data, 2, minne=1000 )
gr3 <- plot.model( data, 3, minne=1000 )
gr4 <- plot.model( data, 4, minne=1000 )
gr5 <- plot.model( data, 5, minne=1000 )
gr <- cbind(gr1,gr2,gr3,gr4,gr5,size="last")
grid.draw(gr)


data2 <- load.from.out( "ceu4-initmigr0.2-iter14.out", int_parameter=1 )
convergence <- subset( data2, epoch>= 12 & epoch <= 17 & type=="Migr" & frm==0 )
qplot( data=convergence, x=iter, y=rate, color=epoch, ylab="migration rate", xlab="iteration", log="y", geom="point" )

emstep <- max( data$iter )
ggsave(paste("migration_",emstep,"EMsteps.png",sep=""),
       gr,
       width = 24, height = 12, units = "cm")



