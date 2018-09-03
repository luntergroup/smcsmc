library(ggplot2)
library(grid)
library(gridExtra)
library(extrafont)

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

msmcpops=c("ceu4.4","ceu4.2a","ceu4.1a","chb4.4","chb4.2a","chb4.1a","yri4.4","yri4.2a","yri4.1a")
msmclabels=c("ceu4-vb","ceu4-vb","ceu4-vb","chb4-vb","chb4-vb","chb4-vb","yri4-vb","yri4-vb","yri4-vb")
msmcdata = NULL
for (popidx in 1:length(msmcpops)) {
    pop = msmcpops[popidx]
    msmcdata <- load.msmc( paste("data/msmc/",pop,".final.txt.out",sep=""), data=msmcdata, int_parameter=popidx )
}
msmcdata$Population <- msmclabels[ msmcdata$int_parameter ]
msmcdata$Experiment <- msmcpops[ msmcdata$int_parameter ]
## remove duplicated rows for smoother plot lines
msmcdata = msmcdata[ !(msmcdata$idx %in% seq(11,39,2)), ]

msmcdata[ msmcdata$int_parameter == 3 & msmcdata$start ==0, ]$ne = NA
msmcdata[ msmcdata$int_parameter == 6 & msmcdata$start ==0, ]$ne = NA
msmcdata[ msmcdata$int_parameter == 9 & msmcdata$start ==0, ]$ne = NA

t0 <- 1700
t1 <- 5000000
ne0 <- 2000
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

p <- p + geom_step( data=subset(msmcdata, Experiment=="ceu4.4"), aes(x=start, y=ne, col=Population, linetype="msmc 4 diploid"))
p <- p + geom_step( data=subset(msmcdata, Experiment=="ceu4.2a"), aes(x=start, y=ne, col=Population, linetype="msmc 2 diploid"))
p <- p + geom_step( data=subset(msmcdata, Experiment=="ceu4.1a"), aes(x=start, y=ne, col=Population, linetype="msmc 1 diploid"))

p <- p + geom_step( data=subset(msmcdata, Experiment=="yri4.4"), aes(x=start, y=ne, col=Population, linetype="msmc 4 diploid"))
p <- p + geom_step( data=subset(msmcdata, Experiment=="yri4.2a"), aes(x=start, y=ne, col=Population, linetype="msmc 2 diploid"))
p <- p + geom_step( data=subset(msmcdata, Experiment=="yri4.1a"), aes(x=start, y=ne, col=Population, linetype="msmc 1 diploid"))

p <- p + geom_step( data=subset(msmcdata, Experiment=="chb4.4"), aes(x=start, y=ne, col=Population, linetype="msmc 4 diploid"))
p <- p + geom_step( data=subset(msmcdata, Experiment=="chb4.2a"), aes(x=start, y=ne, col=Population, linetype="msmc 2 diploid"))
p <- p + geom_step( data=subset(msmcdata, Experiment=="chb4.1a"), aes(x=start, y=ne, col=Population, linetype="msmc 1 diploid"))

p <- plot.ribbon( p, subset(quartiles, Population=="chb4-vb"), "blue", g=g, lwd=1 )
p <- plot.ribbon( p, subset(quartiles, Population=="ceu4-vb"), "darkcyan", g=g, lwd=1 )
p <- plot.ribbon( p, subset(quartiles, Population=="yri4-vb"), "darkred", g=g, lwd=1 )

p <- p + ylab("Ne") + xlab("time [years ago]")
p <- p + theme(text = element_text(size=16), legend.position="top")
p <- p + theme(axis.title.x = element_blank())
p <- p + annotation_logticks(sides="bl", mid=unit(0.1,"cm"),long=unit(0.1,"cm"))
p <- p + scale_colour_manual("",
                             breaks = c("ceu4-vb","chb4-vb","yri4-vb"),
                             values = c("ceu4-vb"="blue","chb4-vb"="darkred","yri4-vb"="darkcyan"),
                             labels = c("CEU","CHB","YRI"))
p <- p + scale_linetype_manual(values=c("msmc 4 diploid"="solid","msmc 2 diploid"="22","msmc 1 diploid"="11"),
                               breaks=c("msmc 4 diploid","msmc 2 diploid","msmc 1 diploid"),labels=c("msmc 4 diploid","msmc 2 diploid","msmc 1 diploid"))


if (talk) {  ## for talk
    ggsave(paste("singlepops_talk.pdf",sep=""),
       p,
       width = 20, height = 10, units = "cm")
}
    

# publication tweaks
p <- p + scale_colour_Publication( labels=c("CEU","CHB","YRI") )
p <- p + theme_Publication(base_size=10,base_family="Helvetica",
                                                        legend.direction=c("vertical"),
                                                        legend.position=c(.6,.8),
                                                        legend.title=element_blank(),
                                                        legend.key.size=0.4)

ggsave(paste("singlepops_",emstep,"EMsteps.pdf",sep=""),
       p,
       width = 9, height = 9, units = "cm")

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
