#!/bin/bash
#$ -cwd
#$ -V
#$ -P bsg.prjb -q long.qb
#$ -e CEUCHB_priorNe_ErrFiles
#$ -o CEUCHB_priorNe_OutFiles
#$ -N CEUCHB_priorNe
#$ -t 1-15
#$ -j y

rep=1
#rep=$(expr $SGE_TASK_ID )

source CEUCHB_priorNe_param.src

data=${fileprefix}".vcf"
echo ${data}
Np=1000
EM=15

prefix=CEUCHB_priorNe_${nsam}sampleNp${Np}
mkdir ${prefix}

outprefix=${prefix}rep${rep}

pattern="1*3+25*2+1*4+1*6"
tmax=8

pf-ARG -l 0 ${cmd} -EM ${EM} -Np ${Np} -o ${outprefix} -p ${pattern} -tmax ${tmax} -vcf ${data} -seed ${rep}

# should add pattern to this!
