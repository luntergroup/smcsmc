library(ggplot2)
library(grid)
library(gridExtra)
library(extrafont)

##font_import() ## once

source("stat-stepribbon.R")
source("plot-utils.R")
source("plot-themes.R")

schiffels.ne <- function( t, n0=1 ) {
    T = c(0.000582262,0.00232905,0.00931619,0.0372648,0.149059,0.596236,1e99)
    e = c(0,          1318.18,   -329.546,  82.3865,  -20.5966,5.14916, 0)
    ne = n0
    i = 1
    t0 = 0
    while (t0 < t) {
        t1 = min( T[i], t )
        ne = ne * exp( -e[i] * (t1 - t0) )
        t0 = t1
        i = i+1
    }
    return (ne)
}

plot.smcsmc <- function() {

    g = 30
    data = NULL
    truth = NULL
    data <- read.table(file="data/zigzag.dta", header=TRUE)
    N0 = 14312
    unscaled.t =  10^(c(250:650)/100)/(g*4*N0) 
    truth = data.frame( unscaled.t * 4 * N0, N0 * apply( matrix(unscaled.t), 1, schiffels.ne ) )
    colnames(truth) <- c("start","Ne")
    truth = rbind( data.frame(truth, int_parameter=0), data.frame(truth, int_parameter=1) )

    t0 <- 500
    t1 <- 1500000
    nemax <- 1500000
    emstep <- max( data$iter )
    emstep=72

    ## calculate.quartiles uses int_parameter as key, so move vb/not vb into there
    data$int_parameter <- 0
    data$int_parameter[ grep("True", data$str_parameter) ] = 1

    plot.data <- subset(data, iter==emstep & type=="Coal" )
    quartiles <- calculate.quartiles( plot.data, min=100, max=nemax, mint=t0/g, maxt=t1/g )
    ##quartiles$Population <- paste( "Pop", quartiles$frm + 1, sep="" )
    quartiles$Population <- "Inferred"
    quartiles$Method <- "SEM"
    quartiles$Method[ quartiles$int_parameter == 1 ] <- "VB"
        
    p <- ggplot(data=quartiles, aes(x=start*g))
    ##    p <- p + scale_x_log10(limits=c(t0,t1),breaks=c(1000,10000,100000,1000000))
    ##    p <- p + scale_y_log10(limits=c(1000,nemax),breaks=c(2000,5000,10000,20000,50000,100000))
    ##p <- p + scale_x_log10(limits=c(t0,t1),breaks=c(1000,10000,100000,1000000), labels=c(expression(10^3),expression(10^4),expression(10^5),expression(10^6)))
    ##p <- p + scale_y_log10(limits=c(minne,maxne),breaks=c(10000,100000,1000000), labels=c(expression(10^4),expression(10^5),expression(10^6)))
    ##p <- p + scale_y_continuous(breaks=c(2000,5000,10000,20000,50000,100000),labels=c("2000","5000","10,000","20,000","50,000","100,000"))

    p <- p + coord_trans(x="log10",y="log10", limx=c(t0,t1),limy=c(900,nemax * 0.90))
    p <- p + scale_x_continuous(breaks=c(1000,10000,100000,1000000),labels=c(expression(10^3),expression(10^4),expression(10^5),expression(10^6)))
    p <- p + scale_y_continuous(breaks=c(1000,10000,100000,1000000),labels=c(expression(10^3),expression(10^4),expression(10^5),expression(10^6)))
    ##p <- p + scale_x_log10(limits=c(t0,t1),breaks=c(1000,10000,100000,1000000), labels=c(expression(10^3),expression(10^4),expression(10^5),expression(10^6)))
    ##p <- p + scale_y_log10(limits=c(1000,nemax * 0.97),breaks=c(1000,10000,100000,1000000), labels=c(expression(10^3),expression(10^4),expression(10^5),expression(10^6)))

    p <- p + geom_ribbon(data=subset(quartiles, Population=="Inferred"), aes(ymin=Q0, ymax=Q4),
                         stat="stepribbon", fill="blue", alpha=.15)
    p <- p + geom_ribbon(data=subset(quartiles, Population=="Inferred"), aes(ymin=Q1, ymax=Q3),
                         stat="stepribbon", fill="blue", alpha=.25)
    p <- p + geom_step(data=subset(quartiles, Population=="Inferred"), aes( y=Q2, x=start*g, colour=Population), lwd=1 )
    p <- p + geom_step(data=truth, aes( y=Ne, x=start*g, colour="Truth"), lwd=0.25 )
    p <- p + annotation_logticks(sides="bl", short=unit(0.05,"cm"), mid=unit(0.05,"cm"),long=unit(0.05,"cm"))    
    p <- p + facet_grid( Method ~ . )
    p <- p + xlab("time [years ago]") + ylab("Ne") + theme(text = element_text(size=16), legend.position="top")
    p <- p + scale_colour_manual("",
                                 breaks = c("Inferred","Truth"),
                                 values = c("Inferred"="blue","Truth"="black"))
    p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 0, size=9))
    p <- p + theme(axis.text.y = element_text(size=9))    
    ##p <- p + theme( panel.grid.minor = element_line(colour="#f0f0f0")
    p <- p + theme_Publication(base_size=10,base_family="Helvetica",legend.title=element_blank())
    p <- p + theme(axis.line.x=element_blank())
    p <- p + theme(axis.line.y=element_blank())

    p
    ggsave(paste("zigzag_",emstep,"EMsteps.pdf",sep=""),
           width = 11, height = 11, units = "cm")    
}

