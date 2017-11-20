library(ggplot2)

source("stat-stepribbon.R")
source("plot-utils.R")

plot.smcsmc <- function( data, truth, g=30 ) {

    #write.table(data, file="cps-apf.dta", sep="\t")
    #data <- read.table(file="cps-apf.dta", header=TRUE)

    t0 <- 5000
    t1 <- 1500000
    minne <- 5000
    maxne <- 1e6
    emstep <- max( data$iter )

    plot.data <- subset(data, iter==emstep & type=="Coal")
    quartiles <- calculate.quartiles( plot.data, min=minne, max=maxne, mint=t0/g, maxt=t1/g )
    quartiles$Population <- "Inferred"

    quartiles$phased <- factor( -quartiles$int_parameter, labels = c("Phased","Unphased") )
    quartiles$aux_part_filt <- factor( quartiles$aux_part_filt, labels = c("","APF"))
    p <- ggplot(data=quartiles, aes(x=start*g))
    p <- p + scale_x_log10(limits=c(t0,t1),breaks=c(1000,10000,100000,1000000), labels=c("1k","10k","100k","1M"))
    p <- p + scale_y_log10(limits=c(minne,maxne),breaks=c(10000,100000,1000000), labels=c("10,000","100,000","1,000,000"))
    p <- p + geom_ribbon(aes(ymin=Q0, ymax=Q4),
                         stat="stepribbon", fill="blue", alpha=.15)
    p <- p + geom_ribbon(aes(ymin=Q1, ymax=Q3),
                         stat="stepribbon", fill="blue", alpha=.25)
    p <- p + geom_step(aes( y=Q2, x=start*g, colour=Population), lwd=0.5 )
    p <- p + geom_step(data=truth, aes( y=Ne, x=t*g, colour="Truth"), lwd=0.25 )
    p <- p + xlab("years ago") + ylab("Ne") + theme(text = element_text(size=16), legend.position="top")
    p <- p + scale_colour_manual("",
                                 breaks = c("Inferred","Truth"),
                                 values = c("Inferred"="blue","Truth"="black"))
    p <- p + facet_grid( phased + aux_part_filt ~ np )
    p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 0, size=9))
    p <- p + theme(axis.text.y = element_text(size=9))    
    p
    ggsave("cps_apf_phasing.png",
           width = 15, height = 15, units = "cm")
}

