#!/bin/bash
#$ -cwd
#$ -V
#$ -P bsg.prjb -q long.qb
#$ -e CEU_CHROM1_2sam_10000npErrFiles
#$ -o CEU_CHROM1_2sam_10000npOutFiles
#$ -N CEU_CHROM1_2sam_10000np
#$ -t 1-15
#$ -j y


#rep=1
rep=$(expr $SGE_TASK_ID )

#nsam=2
Np=10000
EM=20
#EM=0
pattern="1*3+28*2+1*4+1*6"
tmax=8
cmd="-EM ${EM} -Np ${Np} -t 250000 -r 100000 250000000 -p ${pattern} -tmax ${tmax} -l 0 -seed ${rep} -online"

case="CEU_CHROM1_2sam"
#for rep in $(seq 1 1 15)
    #do
#./pop_struct.py ${case} ${nsam} ${rep}
    #done
prefix=${case}rep${rep}
outprefix=${prefix}Np${Np}
pf-ARG ${cmd} -filter 100 -missing 1000000 -o ${outprefix}  -vcf CEU_2sample_CHROM1.vcf -log 
