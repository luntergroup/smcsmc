library(ggplot2)
library(grid)
library(gridExtra)
library(extrafont)

##font_import() ## once

source("stat-stepribbon.R")
source("plot-utils.R")
source("plot-themes.R")

talk = TRUE
pops=c("ceu4","chb4","yri4","ceu4-vb","chb4-vb","yri4-vb")
labels=pops
data=NULL
for (popidx in 1:(length(pops))) {
    pop = pops[popidx]
    data <- load.from.out( paste("data/",pop,"/result.out",sep=""), data=data, int_parameter=popidx )
}
data$Population <- pops[ data$int_parameter ]

t0 <- 3000
t1 <- 5000000
ne0 <- 3000
ne1 <- 1000000
if (talk) {
    ne0 <- 1000
    t1 <- 3000000
}
g <- 29
N0 <- 14312
tbreaks=c(1000,10000,100000,1000000)
nebreaks=c(1000,10000,100000,1000000)
emstep <- max( data$iter )

plot.data <- subset(data, iter==emstep & type=="Coal")
quartiles <- calculate.quartiles( plot.data, miny=ne0, maxy=ne1, mint=t0/g, maxt=t1/g )

p <- ggplot(data=quartiles, aes(x=start*g))
xlabels <- waiver()
ylabels <- waiver()
xlabels <- c(expression(10^3),expression(10^4),expression(10^5),expression(10^6))
ylabels <- c(expression(10^3),expression(10^4),expression(10^5),expression(10^6))
p <- p + scale_x_log10(limits=c(t0,t1),breaks=tbreaks, labels=xlabels)
p <- p + scale_y_log10(limits=c(ne0,ne1),breaks=nebreaks, labels=ylabels)
p <- plot.ribbon( p, subset(quartiles, Population=="ceu4-vb"), "darkcyan", g=g, lwd=0.5 )
p <- plot.ribbon( p, subset(quartiles, Population=="chb4-vb"), "blue", g=g, lwd=0.5 )
p <- plot.ribbon( p, subset(quartiles, Population=="yri4-vb"), "darkred", g=g, lwd=0.5 )

p <- p + ylab("Ne") + xlab("time [years ago]")
p <- p + theme(text = element_text(size=16), legend.position="top")
p <- p + theme(axis.title.x = element_blank())
p <- p + annotation_logticks(sides="bl", mid=unit(0.1,"cm"),long=unit(0.1,"cm"))
p <- p + scale_colour_manual("",
                             breaks = c("ceu4-vb","chb4-vb","yri4-vb"),
                             values = c("ceu4-vb"="blue","chb4-vb"="darkred","yri4-vb"="darkcyan"),
                             labels = c("CEU","CHB","YRI"))

if (talk) {  ## for talk
    ggsave(paste("singlepops_talk.pdf",sep=""),
       p,
       width = 20, height = 10, units = "cm")
}
    

# publication tweaks
p <- p + scale_colour_Publication( labels=c("CEU","CHB","YRI"))
p <- p + theme_Publication(base_size=10,base_family="Helvetica",
                                                        legend.direction=c("vertical"),
                                                        legend.position=c(.6,.8),
                                                        legend.title=element_blank(),
                                                        legend.key.size=0.4)

ggsave(paste("singlepops_",emstep,"EMsteps.pdf",sep=""),
       p,
       width = 9, height = 7, units = "cm")

##
## make plot showing effect of VB.  Incomparable runs; need to re-run
##

onlyvb = TRUE
quartiles2 = quartiles
quartiles2$vb = FALSE; quartiles2$vb[ grep("-vb",quartiles2$Population) ] =TRUE
linetype = "solid2"
if (onlyvb) {
    quartiles2 <- quartiles2[ quartiles2$vb , ]
    linetype = "solid"
}
quartiles2$int_parameter[ quartiles2$int_parameter > 3 ] = quartiles2$int_parameter[ quartiles2$int_parameter > 3 ] -3
quartiles2$Population <- pops[ quartiles2$int_parameter ]

p <- ggplot(data=quartiles2, aes(x=start*g))
xlabels <- c(expression(10^3),expression(10^4),expression(10^5),expression(10^6))
ylabels <- c(expression(10^3),expression(10^4),expression(10^5),expression(10^6))
p <- p + scale_x_log10(limits=c(t0,t1),breaks=tbreaks, labels=xlabels)
p <- p + scale_y_log10(limits=c(ne0,ne1),breaks=nebreaks, labels=ylabels)
p <- plot.ribbon( p, subset(quartiles2, Population=="ceu4"), "darkcyan", g=g, lwd=0.5, linetype=linetype )
p <- plot.ribbon( p, subset(quartiles2, Population=="chb4"), "blue", g=g, lwd=0.5, linetype=linetype )
p <- plot.ribbon( p, subset(quartiles2, Population=="yri4"), "darkred", g=g, lwd=0.5, linetype=linetype )

p <- p + ylab("Ne") + xlab("time [years ago]")
p <- p + theme(text = element_text(size=16), legend.position="top")
p <- p + theme(axis.title.x = element_blank())
p <- p + annotation_logticks(sides="bl", mid=unit(0.1,"cm"),long=unit(0.1,"cm"))


# publication tweaks
p <- p + scale_colour_Publication( labels=c("CEU","CHB","YRI"))
p <- p + scale_linetype_manual(breaks=c(TRUE,FALSE),values=c("11","solid"),labels=c("VB","EM"))
p <- p + theme_Publication(base_size=10,base_family="Helvetica",
                                                        legend.direction=c("vertical"),
                                                        legend.position=c(.6,.8),
                                                        legend.title=element_blank(),
                                                        legend.key.size=0.4)

ggsave(paste("singlepops_vb_em",emstep,"EMsteps.pdf",sep=""),
       p,
       width = 11, height = 11, units = "cm")
