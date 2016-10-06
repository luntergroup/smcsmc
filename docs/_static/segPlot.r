rm(list=ls())

png("segPlot.png", width = 800, height = 300, bg = "transparent")
plot(c(1,9), c(1,3),type="n", axes=FALSE, frame.plot=F, xlab = "", ylab="")
#grid(10, 5, col = "lightgray", lty = "dotted",      lwd = par("lwd"), equilogs = TRUE)

startx = 2
starty = 2
endx = 8
endy = 2
points( startx, starty, pch  = 1, cex = 3)
points( endx, endy, pch  = 16, cex = 3 )
lines( c(startx,endx), c(starty,endy), lty = 2, cex = 3 )
upy = downy = 0.5
leftx = rightx = 0.5

text( startx, starty-downy, label = "segment start", cex = 3 )
text( endx, endy-downy, label = "segment end", cex = 3 )
dev.off()
