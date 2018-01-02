library(ggplot2)

source("stat-stepribbon.R")
source("plot-utils.R")

plot.smcsmc <- function( data, truth, g=30 ) {

    write.table(data, file="zigzag.dta", sep="\t")
    data <- read.table(file="zigzag.dta", header=TRUE)

    t0 <- 500
    t1 <- 1500000
    emstep <- max( data$iter )

    plot.data <- subset(data, iter==emstep & type=="Coal")
    quartiles <- calculate.quartiles( plot.data, min=100, max=1e5, mint=t0/g, maxt=t1/g )
    quartiles$Population <- paste( "Pop", quartiles$frm + 1, sep="" )
    
    p <- ggplot(data=quartiles, aes(x=start*g))
    p <- p + scale_x_log10(limits=c(t0,t1),breaks=c(1000,10000,100000,1000000))
    p <- p + scale_y_log10(limits=c(1000,100000),breaks=c(2000,5000,10000,20000,50000))
    p <- p + geom_ribbon(data=subset(quartiles, Population=="Pop1"), aes(ymin=Q0, ymax=Q4),
                         stat="stepribbon", fill="blue", alpha=.15)
    p <- p + geom_ribbon(data=subset(quartiles, Population=="Pop1"), aes(ymin=Q1, ymax=Q3),
                         stat="stepribbon", fill="blue", alpha=.25)
    p <- p + geom_step(data=subset(quartiles, Population=="Pop1"), aes( y=Q2, x=start*g, colour=Population), lwd=1 )
    p <- p + geom_step(data=truth, aes( y=Ne, x=t*g, colour="Truth"), lwd=0.5 )
    p <- p + xlab("years ago") + ylab("Ne") + theme(text = element_text(size=16), legend.position="top")
    p <- p + scale_colour_manual("",
                                 breaks = c("Pop1","Truth"),
                                 values = c("Pop1"="blue","Truth"="black"))
    p
    ggsave(paste("zigzag_",emstep,"EMsteps.png",sep=""),
           width = 12, height = 8, units = "cm")
    
}

