#!/bin/bash
#$ -cwd
#$ -V
#$ -P bsg.prjb -q long.qb
#$ -e testing_ErrFiles
#$ -o testing_OutFiles
#$ -N testing
#$ -t 1-15
#$ -j y

TASK=$(expr ${SGE_TASK_ID} )
outprefix=${prefix}rep${TASK}
echo ${outprefix}
