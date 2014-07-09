#!/bin/bash
#$ -cwd
#$ -V
#$ -P bsg.prjb -q short.qb
#$ -e testing2_ErrFiles
#$ -o testing2_OutFiles
#$ -N testing2
#$ -t 1-15
#$ -j y

TASK=$(expr ${SGE_TASK_ID} )
outprefix=${prefix}rep${TASK}
echo ${outprefix}

pf-ARG -l 0 -o ${outprefix}


