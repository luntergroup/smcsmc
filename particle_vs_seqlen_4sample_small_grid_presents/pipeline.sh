#!/bin/bash
JOB="particle_vs_seqlen_4sample_smaller_grid"

mkdir /well/gerton/joezhu/${JOB}
DATA=sim-1Samples4msdata1.seg

WTCHG_prj=bsg.prjc

short_queue="short.qc"
long_queue="long.qc"

pattern_tmax="-p \"14*1+8*2+4*4+2*8\" -tmax 8" #1*4+25*2+1*4+1*6 = 64
EM=40
number_of_replicates=15

seqlen_index_array=(0 1 2)
particle_index_array=(0 1 2)
seqlen_index_array=(0 1)
particle_index_array=(0 1)

particle_array=(1000 3000 10000)
seqlen_array=(30000000 120000000 300000000) # 30Mb, 100Mb, 300Mb
theta_array=(30000 120000 300000)
#rho_array=(6000 24000 60000) # correct recombination rate
rho_array=(7000 28000 70000) # Higher recombination rate


for particle_index in "${particle_index_array[@]}"
    do 
    for seqlen_index in "${seqlen_index_array[@]}"
        do
            #current_queue=$([ $(expr $seqlen_index + $particle_index) -gt 0 ] && echo ${long_queue} || echo ${short_queue})
            current_queue=${long_queue} 
            # or alternatively, use the following line
            #[ $(expr $seqlen_index + $particle_index) -gt 1 ] &&    current_queue=${long_queue} ||    current_queue=${short_queue}
            #echo $(expr $seqlen_index + $particle_index) ${current_queue}
            currentjob=Particle${particle_array[particle_index]}Seqlen${seqlen_array[seqlen_index]}_4sample_smaller_grid
            #echo "particles ${particle_array[particle_index]}, seqlen ${seqlen_array[seqlen_index]}, mutation rate ${theta_array[seqlen_index]}, recombination rate ${rho_array[seqlen_index]}"
script="#!/bin/bash
#$ -cwd
#$ -V
#$ -P ${WTCHG_prj} -q ${current_queue}
#$ -e ${JOB}${currentjob}_ErrFiles
#$ -o ${JOB}${currentjob}_OutFiles
#$ -N ${JOB}${currentjob}
#$ -t 1-${number_of_replicates}
#$ -j y
TASK=\$(expr \${SGE_TASK_ID} )
outprefix=/well/gerton/joezhu/${JOB}/${currentjob}TASK\${TASK}

{ time -p pf-ARG  -nsam 4 -Np ${particle_array[particle_index]} -EM ${EM} -t ${theta_array[seqlen_index]} -r ${rho_array[seqlen_index]} ${seqlen_array[seqlen_index]} ${pattern_tmax} -seg ${DATA} -o \${outprefix} -seed \${TASK} \${TASK} \${TASK} ;} 2> \${outprefix}_time.txt
"     
        echo "${script}" > ${currentjob}.sh
        
        chmod a+x ${currentjob}.sh
        qsub ${currentjob}.sh
        done
    done

