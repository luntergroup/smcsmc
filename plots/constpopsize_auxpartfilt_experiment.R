library(ggplot2)
library(grid)
library(gridExtra)
library(extrafont)

##font_import() ## once

source("stat-stepribbon.R")
source("plot-utils.R")
source("plot-themes.R")

data <- read.table(file="data/cps-apf.dta", header=TRUE)
truth <- data.frame( matrix( c((10^(c(25:65/10.0))/(30)), rep(10000,41)), ncol=2) )
names(truth) <- c("t","Ne")

talk = FALSE
g=30
t0 <- 3000
t1 <- 1500000
minne <- 3000
maxne <- 1.3e6
emstep <- max( data$iter )

if (talk) {
    plot.data <- subset(data, iter==emstep & type=="Coal" & int_parameter == 0)
} else {
    plot.data <- subset(data, iter==emstep & type=="Coal")
}

quartiles <- calculate.quartiles( plot.data, min=minne, max=maxne, mint=t0/g, maxt=t1/g )
quartiles$Population <- "Inferred"

quartiles$phased <- factor( quartiles$int_parameter, labels = c("unphased","phased") )
quartiles$aux_part_filt <- factor( quartiles$aux_part_filt, labels = c("","lookahead"))
quartiles$np_label <- factor( quartiles$np, labels=c("N=1,000","N=3,000","N=10,000","N=30,000"))

p <- ggplot(data=quartiles, aes(x=start*g))
p <- p + scale_x_log10(limits=c(t0,t1),breaks=c(1000,10000,100000,1000000), labels=c(expression(10^3),expression(10^4),expression(10^5),expression(10^6)))
p <- p + scale_y_log10(limits=c(minne,maxne),breaks=c(10000,100000,1000000), labels=c(expression(10^4),expression(10^5),expression(10^6)))
p <- p + geom_ribbon(aes(ymin=Q0, ymax=Q4), stat="stepribbon", fill="blue", alpha=.15)
p <- p + geom_ribbon(aes(ymin=Q1, ymax=Q3), stat="stepribbon", fill="blue", alpha=.25)
p <- p + geom_step(aes( y=Q2, x=start*g, colour=Population), lwd=0.5 )
p <- p + geom_step(data=truth, aes( y=Ne, x=t*g, colour="Truth"), lwd=0.25 )
p <- p + annotation_logticks(sides="bl", short=unit(0.05,"cm"), mid=unit(0.05,"cm"),long=unit(0.05,"cm"))
p <- p + xlab("time [years ago]")
p <- p + ylab("Ne")
p <- p + theme(text = element_text(size=16), legend.position="top")
#p <- p + theme(axis.title.x = element_blank())
p <- p + scale_colour_manual("",
                             breaks = c("Inferred","Truth"),
                             values = c("Inferred"="blue","Truth"="black"))
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 0, size=9))
p <- p + theme(axis.text.y = element_text(size=9))    
#p <- p + theme( panel.grid.minor = element_line(colour="#f0f0f0")
p <- p + theme_Publication(base_size=10,base_family="Helvetica",legend.title=element_blank())
p <- p + theme(axis.line.x=element_blank())
p <- p + theme(axis.line.y=element_blank())

if (talk) {
    p <- p + facet_grid( aux_part_filt ~ np )
    ggsave("cps_apf_phasing.pdf",
           width = 13, height = 7, units = "cm")
} else {
    p <- p + facet_grid( aux_part_filt + phased ~ np_label)
    ggsave("cps_apf_phasing.pdf",
           width = 11, height = 11, units = "cm")
}


# <- plot.ribbon( p, subset(quartiles, Population=="ceu4-vb"), "darkcyan", g=g, lwd=0.5 )
#p <- plot.ribbon( p, subset(quartiles, Population=="chb4-vb"), "blue", g=g, lwd=0.5 )
#p <- plot.ribbon( p, subset(quartiles, Population=="yri4-vb"), "darkred", g=g, lwd=0.5 )

#p <- p + scale_colour_manual("",
#                             breaks = c("ceu4-vb","chb4-vb","yri4-vb"),
#                             values = c("ceu4-vb"="blue","chb4-vb"="darkcyan","yri4-vb"="darkred"),
#                             labels = c("CEU","CHB","YRI"))

# publication tweaks
##p <- p + scale_colour_Publication( labels=c("inferred","truth"))
                           #legend.direction=c("vertical"),
                           #legend.position=c(.6,.8),
                           #legend.key.size=0.4)
