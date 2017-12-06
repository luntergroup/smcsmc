library(reshape2)
library(plyr)
library(monomvn)

## use clump data (as opposed to aggregate data)?
use.clumps = TRUE

## maximum epoch distance for explanatory variables
min.epoch = -40
max.epoch = 40

## select part of the data, and setup explanatory (x) and response (y) variables
start = 0
end = 999

#############################################################################################################

## recast data, including the EM outcome column (Clump == -1)
infile <- "zz-clump.out"
infile <- "udm-apf-clump.out"
data <- read.table(infile,header=T)

## !! Clean (for udm-apf-clump.out) -- missing data on iteration 39 (why?)
data <- data[ data$Iter != 39, ]

## add standard error in measurement
data$sd <- sqrt( pmax(0,data$Count) ) / pmax(1.0,data$Opp)

## recast
recastvar <- function( data, var) { return(dcast( data=data, Iter + Clump ~ Type + From + Epoch, value.var=var,
                                                 subset=.(Type %in% c("Coal","Migr")))) }
data.recast <- recastvar( data, "Rate" )
data.recast.sd <- recastvar( data, "sd" )
data.recast.count <- recastvar( data, "Count" )

## remove variables that have little or no counts at end of iterations -- either pegged to 0 or converged to 0
clumps <- 1 + max(data.recast$Clump)
iters <- min( end, max(data.recast$Iter))
keepcols <- c(1, 2, 2 + which(data.recast.count[ data.recast$Clump == -1 & data.recast$Iter == iters , c(-1,-2) ] > 0.5))
data.recast <- data.recast[ , keepcols ]
data.recast.sd <- data.recast.sd[ , keepcols ]
dim0 <- dim(data.recast)[2] - 2  ## number of inferred variables; do not include Iter and Clump

## sanity check - are all data present?
mlt <- melt(data.recast, id.vars=c("Iter","Clump"))
if (sum(is.na(mlt$value))>0) {
    stop("Missing data! stopping")
}

## rescale (by fixed scaling factor, relative standard deviation) each variable, to normalize
## variance so that equal penalties can be applied.  This is done per clump, but the data
## used is the aggregated data, so that an additional factor sqrt(clumps) is needed.
scaling0 <- unlist( 1.0 / (data.recast.sd[ data.recast$Iter == iters & data.recast$Clump == -1, c(-1,-2) ]))
if (use.clumps)  scaling0 <- scaling0 * sqrt(clumps)
scaling <- c(1, 1, scaling0)
data.recast <- as.data.frame( t(t(data.recast) * scaling ) )

## for the response data, remove the EM outcome data ($Clump == -1), and remove
## iteration and Clump columns (columns 1 and 2)
if (use.clumps) {
    clumps.to.use <- 0:(clumps-1)
} else {
    clumps.to.use <- c(-1)
}
ydata <- data.recast[ data.recast$Iter > start &
                      data.recast$Iter <= iters &
                      data.recast$Clump %in% clumps.to.use, 3:length(data.recast) ]

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
fit.mx.sampled <- list()
maxvars = dim0
mprior = 0.1
samples <- 50
rd <- FALSE ### c(5,200) # first=shape, larger=broader; first/second = mean lambda ~ 1/nonzero terms
for ( idx in 1:dim0 ) {
    epochs <- unlist(lapply(strsplit(names(xdata),"_"), function(x) as.integer(x[3])))
    epoch <- epochs[idx]
    cat("Computing row ",idx," for epoch ",epoch," and measurement ",names(xdata)[idx],"\n")
    keep.epochidx <- which( epochs >= epoch + min.epoch & epochs <= epoch + max.epoch )
    xdata.nulled <- xdata
    xdata.nulled[ , -keep.epochidx ] <- 0
    model <- blasso( as.matrix(xdata.nulled), ydata[,idx], lambda2=0, RJ=TRUE, mprior=mprior, rd=rd,
                     theta=0, icept=TRUE, normalize=FALSE, verb=1, M=maxvars, T=samples )
    fit.b <- c( fit.b, median( model$mu[ floor(samples/2):samples ] ) )
    ## choosing the mean in fit.mx appears to give better results
    fit.mx <- rbind( fit.mx, apply( model$beta[ floor(samples/2):samples, ], 2, mean ) )
    fit.mx.sampled[[ idx ]] <- model$beta[ floor(samples/2):samples, ]
    colnames(fit.mx) <- colnames(xdata)
}

## heatmap
sq.mx <- fit.mx[1:dim0,1:dim0]
colnames(sq.mx) <- 1:dim0
sq.mx$X <- 1:dim0
melted.sq.mx <- melt(sq.mx, id="X")
sq.mx$X <- NULL
library(ggplot2)
ggplot(data = melted.sq.mx, aes(x=X, y=variable, fill=abs(value))) + 
  geom_tile()

## compute final answer
final.answer <- function( iter.final, sample=NULL ) {
    if (is.null(sample)) {
        A <- matrix( unlist(fit.mx[1:dim0,1:dim0]), nrow=dim0)
    } else {
        A <- matrix( unlist( lapply(fit.mx.sampled, function(mx) mx[sample,])), nrow=dim0, byrow=TRUE)
    }
    x.final <- c(0)
    y.final <- c(0)
    for (it in iter.final) {
        x.final <- x.final + unlist( data.recast[ data.recast$Iter == it-1 & data.recast$Clump == -1, c(-1,-2) ] )
        y.final <- y.final + unlist( data.recast[ data.recast$Iter == it & data.recast$Clump == -1, c(-1,-2) ] )
    }
    x.final <- x.final / length(iter.final)
    y.final <- y.final / length(iter.final)
    result <- solve( diag(dim0) - A, y.final - A %*% x.final )
    return (result / scaling0)
}

em.answer <- function( iter.final ) {
    y.final <- unlist( data.recast[ data.recast$Iter == iter.final & data.recast$Clump == -1, c(-1,-2) ] )
    return( y.final / scaling0 )
}


## plot the final EM result
plot( 0.5/em.answer( iters ), col="darkred", lwd=2, type="l", lty=1, ylim=c(0,20000) )
text( rep(100, length(xdata)), names(xdata), srt=90, cex=0.75, adj=c(0,0) )

its <- c(iters)
#mt <- apply( matrix(its), 1, function(x) 0.5/final.answer(x) )
#matplot( mt, type="b", pch=1, col=1:length(its), ylim=c(0,20000), lty=3 )
#mt2 <- apply( matrix(its), 1, function(x) 0.5/em.answer(x) )
#matlines( mt2, type="l", pch=1, col=1:length(its), lty=1, lwd=2 )

## the averaged answer seems pretty good
##useiters <- c(iters=10,iters-9,iters-8,iters-7,iters-6,iters-5,iters-4,iters-3,iters-2,iters-1,iters)
useiters <- c(1:iters)
lines( 0.5 / final.answer( useiters ), type="l", lty=5, lwd=3 )

## the resampled answer seems not better than the averaged one
answers <- c()
for (sample in 1:(dim(fit.mx.sampled[[1]])[1])) {
    answers <- c(answers, final.answer( useiters, sample=sample ))
}
real.final.answer <- apply( matrix( answers, ncol=dim0, byrow=TRUE ), 2, median )
lines( 0.5 / real.final.answer, type="l", lty=2, lwd=3, col="darkgreen" )
