library(ggplot2)
library(grid)
library(gridExtra)

source("stat-stepribbon.R")
source("plot-utils.R")

data <- load.from.out( "ceu4.out", int_parameter=1 )
data <- load.from.out( "chb4.out", data=data, int_parameter=2 )
data <- load.from.out( "yri4.out", data=data, int_parameter=3 )
data <- load.from.out( "yri4-vb/emiter11/chunkfinal.out", data=data, int_parameter=4 )
data <- load.from.out( "ceu4-vb/emiter11/chunkfinal.out", data=data, int_parameter=5 )
data$Population <- c("CEU","CHB","YRI","CEUvb","YRIvb")[ data$int_parameter ]

t0 <- 5000
t1 <- 2500000
g <- 29
N0 <- 14312
emstep <- max( data$iter )
emstep <- 11

plot.data <- subset(data, iter==emstep & type=="Coal")
quartiles <- calculate.quartiles( plot.data, min=100, max=1e7, mint=t0/g, maxt=t1/g )

p <- ggplot(data=quartiles, aes(x=start*g))
labels <- waiver()
p <- p + scale_x_log10(limits=c(t0,t1),breaks=c(1000,3000,10000,30000,100000,300000,1000000), labels=labels)
p <- p + scale_y_log10(limits=c(2000,1000000),breaks=c(2000,5000,10000,20000,50000,100000,200000,500000,1000000), labels=labels)
p <- plot.ribbon( p, subset(quartiles, Population=="CEU"), "darkcyan" )
p <- plot.ribbon( p, subset(quartiles, Population=="CHB"), "blue" )
p <- plot.ribbon( p, subset(quartiles, Population=="YRI"), "darkred" )

p <- p + ylab("Ne")
p <- p + theme(text = element_text(size=16), legend.position="top")
p <- p + theme(axis.title.x = element_blank())
p <- p + scale_colour_manual("",
                             breaks = c("CEU","CHB","YRI","CEUvb","YRIvb"),
                             values = c("CEU"="blue","CHB"="darkcyan","YRI"="darkred","CEUvb"="cyan","YRIvb"="red"))

ggsave(paste("singlepops_",emstep,"EMsteps.png",sep=""),
       p,
       width = 12, height = 8, units = "cm")
