rm ( list = ls() )

#./pf-ARG_prof -l 0 -Np 10 -lag 100000 -seed 1
# usage is constant at 2436kb, replace the number of particle to get the y array"
n = c(10, 100, 500, 1000, 3000, 5000, 10000)
y = c(2436, 2936, 7684, 13492, 36920, 60088, 117916)

pdf("usage.pdf")
lm1 = lm (y~n)
overhead = lm1$coefficient[1]
usage = lm1$coefficient[2]
plot(n,y, main = paste(usage,"kb per partile ", overhead," kb overhead") )
ylabel = ( "Memory ( kb )" )
xlabel = ( "Number of particles" )
abline (lm(y~n))
dev.off()