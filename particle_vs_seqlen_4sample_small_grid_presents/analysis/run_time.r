rm ( list = ls () )
library(ggplot2)
particles_array = c ( 1000, 3000, 10000 ) 
seqlen_array = c ( 30000000, 120000000)
#particles_array = c ( 1000, 3000 ) 
#seqlen_array = c ( 30000000, 120000000, 300000000)
DATA_dir = "../simulation_runs/"

Time_data = data.frame( Particle = numeric(0), Seqlen = numeric(0), Runtime = numeric(0), se = numeric(0) )
entry = 1

for ( particle in particles_array ) {
    for ( seqlen in seqlen_array ) {
        current_job = c()
        for ( task in c(1:15) ) { 
            FILEname = paste ( DATA_dir, "Particle", particle, "Seqlen", format(seqlen, scientific = FALSE), "TASK", task, "_time.txt", sep = "")  
            print(FILEname)
            dummy_time = read.table (FILEname)$V2[2]/3600
            current_job = c(current_job, dummy_time)
        }
        Time_data[entry,]$Particle = particle 
        Time_data[entry,]$Seqlen   = seqlen
        Time_data[entry,]$Runtime  = mean(current_job)
        Time_data[entry,]$se  = sd(current_job)/sqrt(15)
        entry =  entry + 1

    }
}

Time_data$Particle <- factor(Time_data$Particle)
Time_data$Seqlen <- factor(Time_data$Seqlen)

pdf("time_vs_particle.pdf")
ggplot(Time_data, aes(x=Particle, y=Runtime, fill=Seqlen)) +
    ylab("Time / hours")+
    geom_bar(position=position_dodge(), stat="identity") +
    geom_errorbar(aes(ymin=Runtime-se, ymax=Runtime+se),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9))

dev.off()

pdf("time_vs_seqlen.pdf")
ggplot(Time_data, aes(x=Seqlen, y=Runtime, fill=Particle)) +
    ylab("Time / hours")+
    geom_bar(position=position_dodge(), stat="identity") +
    geom_errorbar(aes(ymin=Runtime-se, ymax=Runtime+se),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9))

dev.off()


library(scales)     # Need the scales package

pdf("time_vs_particle_log.pdf")
time_plot = ggplot(Time_data, aes(x=Particle, y=Runtime, fill=Seqlen)) + 
    ylab("Time / hours")+
    geom_bar(position=position_dodge(), stat="identity") +
    geom_errorbar(aes(ymin=Runtime-se, ymax=Runtime+se),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9))
time_plot + scale_y_log10()

dev.off()

pdf("time_vs_seqlen_log.pdf")
time_plot = ggplot(Time_data, aes(x=Seqlen, y=Runtime, fill=Particle)) + 
    ylab("Time / hours")+
    geom_bar(position=position_dodge(), stat="identity") +
    geom_errorbar(aes(ymin=Runtime-se, ymax=Runtime+se),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9))
time_plot + scale_y_log10()

dev.off()
