#!/bin/bash
#$ -cwd
#$ -V
#$ -P bsg.prjb -q long.qb
#$ -e CEUCHB_priorNeNp1000_ErrFiles
#$ -o CEUCHB_priorNeNp1000_OutFiles
#$ -N CEUCHB_priorNeNp1000
#$ -t 1-15
#$ -j y

#rep=1
rep=$(expr $SGE_TASK_ID )

source CEUCHB_priorNe_param.src

data=${fileprefix}".vcf"
echo ${data}
Np=1000
EM=15

prefix=CEUCHB_priorNe_${nsam}sampleNp${Np}
#mkdir ${prefix}

outprefix=${prefix}rep${rep}

pattern="1*3+25*2+1*4+1*6"
tmax=8

infer_mig_pattern="-eM 0.02599 1 -eM 0.03736 1 -eM 0.05044 1 "
pf-ARG -l 0 ${cmd} ${split} ${infer_mig_pattern} ${pop_struct} -EM ${EM} -Np ${Np} -o ${outprefix} -vcf ${data} -seed ${rep}

# should add pattern to this!
