library(ggplot2)
library(grid)
library(gridExtra)
library(reshape)
library(extrafont)

source("stat-stepribbon.R")
source("plot-utils.R")
source("model.R")


pops=c("ceutsi4-vb","ceugih4-vb","ceumxl4-vb","ceuchb4-vb","ceuyri4-vb","ceutsi4","ceugih4","ceumxl4","ceuchb4","ceuyri4")
labels=c("TSI","GIH","MXL","CHB","YRI","TSI","GIH","MXL","CHB","YRI")
data=NULL
for (popidx in 1:(length(pops))) {
    pop = pops[popidx]
    data <- load.from.out( paste("data/",pop,"/result.out",sep=""), data=data, int_parameter=popidx )
}
data$Population <- pops[ data$int_parameter ]
data$Pop1 <- "CEU"
data$Pop2 <- labels[ data$int_parameter ]


plot.model <- function( data, popidx, g=29, minne=100, modeldata=NULL, p.modelvars=NULL, m.modelvars=NULL, splitepoch=21 ) {

    split_time <- 385e3
    t0 <- 3000
    t1 <- 1500000
    ne0 <- 1000
    ne1 <- 1000000
    t1m = split_time
    N0 <- 14312
    maxne <- 3e5
    min.migr.rate <- 2e-2
    emstep <- max( data[ data$int_parameter==popidx, ]$iter )
    ##max.migr.rate <- 5e2   ## for log scale
    max.migr.rate <- 10      ## for continuous scale
    nebreaks=c(1000,10000,100000,1000000)
    xlabels <- c(expression(10^3),expression(10^4),expression(10^5),expression(10^6))
    ylabels <- c(expression(10^3),expression(10^4),expression(10^5),expression(10^6))
    if (is.null(emstep)) emstep <- max( data$iter )

    model <- 2  ## for now

    plot.data <- subset(data, iter==emstep & type=="Coal" & int_parameter==popidx)
    plot.data.migration <- subset(data, iter==emstep & type=="Migr" & int_parameter==popidx)
    
    quartiles <- calculate.quartiles( plot.data, miny=ne0, maxy=ne1, mint=t0/g, maxt=t1/g )
    quartiles$Population <- paste( "Pop", quartiles$frm + 1, sep="" )

    data.pop1 = subset(quartiles, Population=="Pop1")
    data.pop2 = subset(quartiles, Population=="Pop2" & epoch < splitepoch)
    data.pop2 = rbind( data.pop2, data.pop1[ data.pop1$epoch == splitepoch, ])
    data.pop2$Population <- "Pop2"

    p <- ggplot(data=quartiles, aes(x=start*g))
    p <- p + scale_x_log10(limits=c(t0,t1),breaks=c(1000,10000,100000,1000000), labels=xlabels)
    p <- p + annotation_logticks(sides="b")
    labels <- NULL
    if (model == 2) labels <- waiver()
    p <- p + scale_y_log10(limits=c(ne0,ne1),breaks=nebreaks, labels=ylabels)
    p <- plot.ribbon( p, data.pop2, "darkcyan" )
    p <- plot.ribbon( p, data.pop1, "blue" )
    p <- p + annotation_logticks(sides="bl", mid=unit(0.1,"cm"),long=unit(0.1,"cm"))
    
    for (m in p.modelvars) {
        p <- p + geom_line(data=subset(modeldata, variable==m) , aes(x=sim.t, y=value))
    }

    if (model == 2) {
        p <- p + theme(text = element_text(size=16)) + ylab("Ne")
    } else {
        p <- p + theme(axis.title.y = element_blank() )
    }
    p <- p + theme(text = element_text(size=16), legend.position="top")
    p <- p + theme(axis.title.x = element_blank())
    p <- p + scale_colour_manual("",
                                 breaks = c("Pop1", "Pop2","model"),
                                 values = c("Pop1"="blue","Pop2"="darkcyan","model"="black"))
    "migration"
    plot.data <- plot.data.migration
    quartiles <- calculate.quartiles( plot.data, field="rate", mint=t0/g, maxt=t1m/g, miny=0, maxy=100, scale=4*N0/(4*N0*g*1e-6), add=min.migr.rate )
    quartiles$Population <- paste( "Pop", quartiles$frm + 1, sep="" )

    data.pop1 = subset(quartiles, Population=="Pop1")
    data.pop2 = subset(quartiles, Population=="Pop2" & epoch < splitepoch)

    pm <- ggplot(data=quartiles, aes(x=start*g))
    for (m in m.modelvars) {
        pm <- pm + geom_line(data=subset(modeldata,variable==m), aes(x=sim.t, y=pmin(max.migr.rate,pmax(min.migr.rate,value))))
    }
    if (model == 2) {
        pm <- pm + ylab("Migration rate (per My)") + theme(text = element_text(size=16), legend.position="none")
    } else {
        pm <- pm + theme(axis.title.y = element_blank() ) + theme(legend.position="none")
    }
    pm <- pm + scale_x_log10(limits=c(t0,t1),breaks=c(1000,10000,100000,1000000),labels=xlabels)
    pm <- pm + annotation_logticks(sides="b")
    labels <- waiver()
    if (model == 0) labels <- waiver()
    pm <- pm + scale_y_continuous(limits=c(min.migr.rate,max.migr.rate), labels=labels)
    pm <- plot.ribbon( pm, data.pop1, "blue", maxy=max.migr.rate )
    pm <- plot.ribbon( pm, data.pop2, "darkcyan", maxy=max.migr.rate )

    pm <- pm + xlab("years ago")
    pm <- pm + scale_colour_manual("",
                                 breaks = c("Pop1","Pop2","model"),
                                 values = c("Pop1"="blue","Pop2"="darkcyan", "model"="black"))


    grob <- rbind(ggplotGrob(p), ggplotGrob(pm), size="last")
    return(grob)
}



plot.model.allemsteps <- function( data, popidx, g=29, minne=2000, max.migr.rate=10, emstep=NULL ) {

    split_time <- 386e3
    splitepoch <- 21
    t0 <- 3000
    t1 <- 1500000
    t1m = split_time
    N0 <- 14312
    maxne <- 3e5
    min.migr.rate <- 2e-2
    ##max.migr.rate <- 5e2   ## for log scale
    max.migr.rate <- 10      ## for continuous scale
    if (is.null(emstep)) emstep <- max( data$iter )

    plot.data <- subset(data, type=="Coal" & int_parameter==popidx)
    plot.data.migration <- subset(data, type=="Migr" & int_parameter==popidx)
    
    quartiles <- calculate.quartiles( plot.data, miny=minne, maxy=maxne, mint=t0/g, maxt=t1/g )
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
    quartiles <- calculate.quartiles( plot.data, field="rate", mint=t0/g, maxt=t1m/g, miny=0, maxy=max.migr.rate, scale=4*N0/(4*N0*g*1e-6), add=min.migr.rate )
    quartiles$Population <- paste( "Pop", quartiles$frm + 1, sep="" )

    data.pop1 = subset(quartiles, Population=="Pop1")
    data.pop2 = subset(quartiles, Population=="Pop2" & epoch < splitepoch)

    pm <- ggplot(data=quartiles, aes(x=start*g, colour=iter))
    pm <- pm + ylab("Migration rate (per My)") + theme(text = element_text(size=16), legend.position="none")
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






data <- load.from.out( "data/ceuchb4-vb/result.out", int_parameter=1 )
data <- load.from.out( "data/ceugih4-vb/result.out", data=data, int_parameter=2 )
data <- load.from.out( "data/ceumxl4-vb/result.out", data=data, int_parameter=3 )
data <- load.from.out( "data/ceutsi4-vb/result.out", data=data, int_parameter=4 )
data <- load.from.out( "data/ceuyri4-vb/result.out", data=data, int_parameter=5 )
data <- load.from.out( "data/ceuyri4-simb-PE-vb/result.out", data=data, int_parameter=6 )
data <- load.from.out( "data/ceuyri4-simb-P-vb/result.out", data=data, int_parameter=7 )
data <- load.from.out( "data/ceuyri4-simb-PL-vb/result.out", data=data, int_parameter=8 )
data$Population <- c("CHB","GIH","MXL","TSI","YRI","YRIsim-PE","YRIsim-P","YRIsim-PL")[ data$int_parameter ]
data$Pop1 <- "CEU"
data$Pop2 <- data$Population

if (FALSE) {
    data <- load.from.out( "data/ceuyri4-vb/result.out", int_parameter=1 )
    data <- load.from.out( "data/ceuyri4-simb-P-vb/result.out", data=data, int_parameter=2 )
    data <- load.from.out( "data/ceuyri4-simb-PE-vb/result.out", data=data, int_parameter=3 )
    data <- load.from.out( "data/ceuyri4-simb-PL-vb/result.out", data=data, int_parameter=4 )
    data$Population <- c("YRI","P","PE","PL")[ data$int_parameter ]
    data <- data[ data$iter <= 7 , ]
}


g <- 29
sim.t <- exp( seq( log(133),log(385e3/g),length.out=100) )
N0 <- 14312
modeldata <- data.frame(sim.t= sim.t * g,
                        pop.ceu=ceu(sim.t / 1.12),
                        pop.yri=yri(sim.t / 1.12),
                        migr.ceu=apply(matrix(sim.t / 1.12),1,migr_ceu) * 4 * N0,
                        migr.yri.early=apply(matrix(sim.t / 1.12),1,function(x) migr_yri( x, -1 )) * 4 * N0,
                        migr.yri=apply(matrix(sim.t / 1.12),1,function(x) migr_yri(x, 0)) * 4 * N0,
                        migr.yri.late=apply(matrix(sim.t / 1.12),1,function(x) migr_yri(x, 1)) * 4 * N0 )
modeldata <- melt( modeldata, id=c("sim.t") )
                        


gr1 <- plot.model( data, 5, minne=2000)
grid.draw( gr1 )

ggsave("migr-ceuyri-vb.pdf",gr1, width=20, height=20, unit="cm")

gr1 <- plot.model( data, 6, minne=2000, splitepoch=16, modeldata=modeldata, p.modelvars=c("pop.ceu","pop.yri"), m.modelvars=c("migr.ceu","migr.yri.early") )
grid.draw( gr1 )

gr1 <- plot.model( data, 7, minne=2000, splitepoch=16, modeldata=modeldata, p.modelvars=c("pop.ceu","pop.yri"), m.modelvars=c("migr.ceu","migr.yri") )
grid.draw( gr1 )

gr1 <- plot.model( data, 8, minne=2000, splitepoch=16, modeldata=modeldata, p.modelvars=c("pop.ceu","pop.yri"), m.modelvars=c("migr.ceu","migr.yri.late") )
grid.draw( gr1 )


gr2 <- plot.model.allemsteps( data, 5, minne=2000, max.migr.rate=10 )
grid.draw(gr2)

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



