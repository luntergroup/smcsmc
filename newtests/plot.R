library(ggplot2)

df <- read.table("newtests/2sample1.recomb", header=T)
df$iter <- as.factor(df$iter)

part <- df[ (df$locus < 1000000) , ]

# opportunities are quite comparable between EM iterations:
ggplot( part, aes(locus, opp_per_nt/3e6, color=iter)) + geom_step() 

# the opportunities in a facet plot, with raw posterior recombination counts.  Recombinations change more between iters...
ggplot( part, aes(locus))  + geom_step(aes(y=X1+X2)) + geom_step(aes(y=opp_per_nt/3e6), colour="blue") + facet_wrap( ~iter)

# ...but once you smooth it out, clear similarities are visible between iterations...
ggplot( part, aes(x=locus,weight=(X1+X2))) + geom_density( adjust=0.02, fill="blue" ) + facet_wrap(~iter)

# ...although lots of detail is different at smaller scales
ggplot( part, aes(x=locus,weight=(X1+X2))) + geom_density( adjust=0.000001, fill="blue" ) + facet_wrap(~iter)




