#!/bin/bash
#$ -cwd
#$ -V
#$ -P bsg.prjb -q short.qb
#$ -e test-12sample1000npErrFiles
#$ -o test-12sample1000npOutFiles
#$ -N test-12sample1000np
#$ -t 1-15
#$ -j y

rep=$(expr $SGE_TASK_ID )
data_dir="/well/bsg/joezhu/experiment4_data/"
out_dir="/well/bsg/joezhu/experiment4/"

nsam=2
Np=1000
EM=20
pattern="1*4+25*2+1*4+1*6"
tmax=2

cmd="-EM ${EM} -Np ${Np} -t 25000 -r 2500 5000000 -p ${pattern} -tmax ${tmax} -l 0 -seed ${rep}"
case="test-1"

prefix=${case}Samples${nsam}msdata${rep}
outprefix=${out_dir}${prefix}Np${Np}
pf-ARG ${cmd} -o ${outprefix}  -vcf ${data_dir}${prefix}.vcf -log 
