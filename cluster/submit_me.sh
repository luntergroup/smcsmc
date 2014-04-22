#!/bin/bash
#$ -cwd
#$ -V
#$ -P bsg.prjb -q short.qb
#$ -e ErrFiles
#$ -o OutFiles
#$ -N test1_pfARG 
#$ -t 1-10
#$ -j y

source test1_pfARG
rep=$(expr $SGE_TASK_ID )
rep=1

../utility/pop_struct.py ${Case} ${nsample} ${rep}
pf-ARG -EM ${EMsteps} -p ${pattern} -Np ${Nparticle} -l 0 -t 10000.0 -r 10000.0 1000000 -tmax 5.0 -lag 200 -vcf ${Case}Samples${nsample}msdata${rep}.vcf -o ${Case}Samples${nsample}msdata${rep} -seed ${rep}

