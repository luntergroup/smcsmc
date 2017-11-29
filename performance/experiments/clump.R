library(reshape2)
library(plyr)
library(monomvn)

## use per-clump intercept vector?
use.intercept = FALSE

## use clump data (as opposed to aggregate data)?
use.clumps = TRUE

## maximum epoch distance for explanatory variables
min.epoch = -8
max.epoch = 8

#############################################################################################################

## recast data, including the EM outcome column (Clump == -1)
data <- read.table("zz-clump.out",header=T)
recastvar <- function( data, var) { return(dcast( data=data, Iter + Clump ~ Type + From + Epoch, value.var=var,
                                                 subset=.(Type %in% c("Coal","Migr")))) }
data.recast <- recastvar( data, "Rate" )
data$sd <- sqrt( pmax(0,data$Count) ) / pmax(1.0,data$Opp)
data.recast.sd <- recastvar( data, "sd" )

clumps <- 1 + max(data.recast$Clump)
iters <- max(data.recast$Iter)

## rescale (by fixed scaling factor, relative standard deviation) each variable, to normalize
## variance so that equal penalties can be applied.  This is done per clump, but the data
## used is the aggregated data, so that an additional factor sqrt(clumps) is needed.
scaling0 <- unlist( 1.0 / (data.recast.sd[ data.recast$Iter == iters & data.recast$Clump == -1, c(-1,-2) ]))
if (use.clumps)  scaling0 <- scaling0 * sqrt(clumps)
scaling <- c(1, 1, scaling0)
data.recast <- as.data.frame( t(t(data.recast) * scaling ) )

## set up auxiliary response variables, to infer per-clumb b vectors
dim0 <- dim(data.recast)[2] - 2  ## do not include Iter and Clump
if (use.intercept & use.clumps) {
    vars <- paste("Clump", unique(subset(data,data$Clump >= 0)$Clump), sep="_")
    data.recast[vars] <- 0
    for (i in 1:length(vars)) {
        clump <- as.integer(unlist(strsplit(vars[i],"_"))[2])
        data.recast[ data.recast$Clump == clump ,vars[i] ] <- 1
    }
}

## select part of the data, and setup explanatory (x) and response (y) variables
start = 0

## for the response data, remove the EM outcome data ($Clump == -1), and remove
## iteration and Clump columns (columns 1 and 2)
if (use.clumps) {
    clumps.to.use <- 0:(clumps-1)
} else {
    clumps.to.use <- c(-1)
}
ydata <- data.recast[ data.recast$Iter > start & data.recast$Clump %in% clumps.to.use, 3:length(data.recast) ]

## for the explanatory data, only use the EM outcome data, but of the previous iteration
## start by selecting the right iterations and columns; but leave the Iter column in
## then replace the data columns by the EM outcome data of the correct iteration
## finally, remove Iter column
xdata <- data.recast[ data.recast$Iter >= start &
                      data.recast$Iter < iters &
                      data.recast$Clump %in% clumps.to.use, c(1,3:length(data.recast)) ]
for (iter in start:(iters-1)) {
    xdata[ xdata$Iter == iter, 2:(dim0+1) ] <- data.recast[ data.recast$Clump == -1 &
                                                            data.recast$Iter == iter, 3:(dim0+2) ]
}
xdata$Iter <- NULL

## fit the model
fit.mx <- data.frame()
fit.b <- c()
if (use.intercept & use.clumps) {
    maxvars <- clumps+dim0
    mprior = 0.01
    samples <- 30
    rd <- c(0.5,10)
} else {
    maxvars = dim0
    mprior = 0.1
    samples <- 50
    rd <- FALSE ### c(5,200) # first=shape, larger=broader; first/second = mean lambda ~ 1/nonzero terms
}
for ( idx in 1:dim0 ) {
    epoch <- idx - 1
    cat("Computing row for epoch ",epoch,"\n")
    epochs <- unlist(lapply(strsplit(names(xdata),"_"), function(x) as.integer(x[3])))
    keep.epochidx <- which( epochs >= epoch + min.epoch & epochs <= epoch + max.epoch )
    xdata.nulled <- xdata
    xdata.nulled[ , -keep.epochidx ] <- 0
    model <- blasso( as.matrix(xdata.nulled), ydata[,idx], lambda2=0, RJ=TRUE, mprior=mprior, rd=rd,
                     theta=0, icept=TRUE, normalize=FALSE, verb=1, M=maxvars, T=samples )
    fit.b <- c( fit.b, median( model$mu[ floor(samples/2):samples ] ) )
    fit.mx <- rbind( fit.mx, apply( model$beta[ floor(samples/2):samples, ], 2, mean ) )
    colnames(fit.mx) <- colnames(xdata)
}

## heatmap
sq.mx <- fit.mx[1:36,1:36]
colnames(sq.mx) <- 0:35
sq.mx$X <- 1:36
melted.sq.mx <- melt(sq.mx, id="X")
sq.mx$X <- NULL
library(ggplot2)
ggplot(data = melted.sq.mx, aes(x=X, y=variable, fill=abs(value))) + 
  geom_tile()

## compute final answer
final.answer <- function( iter.final ) {
    x.final <- c(0)
    y.final <- c(0)
    for (it in iter.final) {
        x.final <- x.final + unlist( data.recast[ data.recast$Iter == it-1 & data.recast$Clump == -1, c(-1,-2) ] )
        y.final <- y.final + unlist( data.recast[ data.recast$Iter == it & data.recast$Clump == -1, c(-1,-2) ] )
    }
    x.final <- x.final / length(iter.final)
    y.final <- y.final / length(iter.final)
    A <- matrix( unlist(fit.mx[1:dim0,1:dim0]), nrow=dim0)
    result <- solve( diag(dim0) - A, y.final - A %*% x.final )
    return (result / scaling0)
}

em.answer <- function( iter.final ) {
    y.final <- unlist( data.recast[ data.recast$Iter == iter.final & data.recast$Clump == -1, c(-1,-2) ] )
    return( y.final / scaling0 )
}

its <- c(3,5,10,iters-3,iters-2,iters-1,iters)
mt <- apply( matrix(its), 1, function(x) 0.5/final.answer(x) )
mt2 <- apply( matrix(its), 1, function(x) 0.5/em.answer(x) )
matplot( mt, type="b", pch=1, col=1:length(its), ylim=c(0,20000) )
matlines( mt2, type="l", pch=1, col=1:length(its), lty=3 )

lines( 0.5 / final.answer( c(21,22,26,27,28,29,30) ), type="l", lty=5, lwd=3 )
