library(ggplot2)
library(grid)
library(gridExtra)
library(extrafont)

source("stat-stepribbon.R")
source("plot-utils.R")
source("plot-themes.R")

pops=c("sim-ceu-0","sim-yri-0")
labels=c("CEU","YRI")
data=NULL
for (popidx in 1:(length(pops))) {
    pop = pops[popidx]
    data <- load.from.out( paste("data/",pop,"/result.out",sep=""), data=data, int_parameter=popidx )
}
data$Population <- labels[ data$int_parameter ]

## todo::
msmcpops=c("simceu.4","simceu.2","simceu.1","simyri.4","simyri.2","simyri.1")
msmclabels=c("msmc CEU","msmc CEU","msmc CEU","msmc YRI","msmc YRI","msmc YRI")
msmcdata = NULL
for (popidx in 1:length(msmcpops)) {
    pop = msmcpops[popidx]
    msmcdata <- load.msmc( paste("data/msmc/",pop,".final.txt.out",sep=""), data=msmcdata, int_parameter=popidx )
}
msmcdata$Population <- msmclabels[ msmcdata$int_parameter ]
msmcdata$Experiment <- msmcpops[ msmcdata$int_parameter ]
## remove duplicated rows for smoother plot lines
msmcdata = msmcdata[ !(msmcdata$idx %in% seq(11,39,2)), ]


truene <- function( t, pop="ceu" ) {
    if (pop == "ceu") {
        expt = c(0.00120468, 0.0180702,0.180702, 1.084212, 2.710529)   #  2kya, 30kya, 300kya, 1.8Mya, 4.5Mya
        expe = c(214.8965, -14.15827,1.33255, -0.563414)
    } else {
        expt = c(0.00120468, 0.00542106, 0.0451755, 0.180702, 1.084212, 2.710529)
        expe = c(502.8635, 0,  -5.89189, 1.33255, -0.563414)
    }
    return( apply( array(t), 1, function(tt) return( truene.t( tt, expt, expe ) ) ) )
}

truene.t <- function( t, expt, expe ) {
    N0 = 14312
    g = 29
    t = t/(4*N0*g)
    ne = 1.5e5
    idx = 1
    if (t < expt[idx]) return (ne)
    t = t - expt[idx]
    while (idx < length(expt)) {
        dt = min( t, expt[idx+1]-expt[idx] )
        ne = ne * exp( -expe[idx] * dt )
        if (dt == t) return(ne)
        idx = idx + 1
        t = t - dt
    }
    return(ne)
}
    
tt = exp( seq( log(3e3), log(6e6), length=1000 ) )
truth = data.frame( t = tt, ceu = truene( tt, pop="ceu" ), yri = truene( tt, pop="yri" ) )


t0 <- 3000
t1 <- 5000000
ne0 <- 1700
ne1 <- 1000000
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

p <- p + geom_step( data=subset(msmcdata, Experiment=="simceu.4"), aes(x=start, y=ne, col=Population, linetype="4 diploid"))
p <- p + geom_step( data=subset(msmcdata, Experiment=="simceu.2"), aes(x=start, y=ne, col=Population, linetype="2 diploid"))
p <- p + geom_step( data=subset(msmcdata, Experiment=="simceu.1"), aes(x=start, y=ne, col=Population, linetype="1 diploid"))

p <- p + geom_step( data=subset(msmcdata, Experiment=="simyri.4"), aes(x=start, y=ne, col=Population, linetype="4 diploid"))
p <- p + geom_step( data=subset(msmcdata, Experiment=="simyri.2"), aes(x=start, y=ne, col=Population, linetype="2 diploid"))
p <- p + geom_step( data=subset(msmcdata, Experiment=="simyri.1"), aes(x=start, y=ne, col=Population, linetype="1 diploid"))

p <- plot.ribbon( p, subset(quartiles, Population=="CEU"), col="blue", g=g, lwd=1 )
#p <- plot.ribbon( p, subset(quartiles, Population=="ceu1"), "lightblue", g=g, lwd=0.5 )
p <- plot.ribbon( p, subset(quartiles, Population=="YRI"), col="red", g=g, lwd=1 )
#p <- plot.ribbon( p, subset(quartiles, Population=="yri1"), "darkred", g=g, lwd=0.5 )
p <- p + geom_step(data=truth, aes( y=ceu, x=t, colour="truth"), lwd=0.35 )
p <- p + geom_step(data=truth, aes( y=yri, x=t, colour="truth"), lwd=0.35 )

p <- p + ylab("Ne") + xlab("time [years ago]")
p <- p + theme(text = element_text(size=16), legend.position="top")
p <- p + theme(axis.title.x = element_blank())
p <- p + annotation_logticks(sides="bl", mid=unit(0.1,"cm"),long=unit(0.1,"cm"))
p <- p + scale_colour_manual(values = c("CEU"="#386cb0","YRI"="#fdb462","msmc CEU"="#c0e0ff","msmc YRI"="#ffe060","truth"="#000000"))
p <- p + scale_linetype_manual(values=c("4 diploid"="solid","2 diploid"="22","1 diploid"="11"),
                               breaks=c("4 diploid","2 diploid","1 diploid"),labels=c("4 diploid","2 diploid","1 diploid"))

p <- p + guides(color = guide_legend(override.aes = list(size=0.5)))

# publication tweaks
p <- p + theme_Publication(base_size=10,base_family="Helvetica",
                                                        legend.direction=c("vertical"),
                                                        legend.position=c(.6,.8),
                                                        legend.title=element_blank(),
                                                        legend.key.size=0.4)

ggsave(paste("sim_singlepops_",emstep,"EMsteps.pdf",sep=""),
       p,
       width = 9, height = 9, units = "cm")

