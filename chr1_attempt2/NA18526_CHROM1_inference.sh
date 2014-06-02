#!/bin/bash
#$ -cwd
#$ -V
#$ -P bsg.prjb -q short.qb
#$ -e NA18526_CHROM1_1000npErrFiles
#$ -o NA18526_CHROM1_1000npOutFiles
#$ -N NA18526_CHROM1_1000np
#$ -t 1-15
#$ -j y


#rep=1
rep=$(expr $SGE_TASK_ID )

#nsam=2
Np=1000
EM=20
#EM=0
pattern="1*4+25*2+1*4+1*6"
tmax=8
cmd="-EM ${EM} -Np ${Np} -t 250000 -r 100000 250000000 -p ${pattern} -tmax ${tmax} -l 0 -seed ${rep} -online"

case="NA18526_CHROM1"
#for rep in $(seq 1 1 15)
    #do
#./pop_struct.py ${case} ${nsam} ${rep}
    #done
prefix=${case}rep${rep}
outprefix=${prefix}Np${Np}
pf-ARG ${cmd} -filter 100 -missing 1000000 -o ${outprefix}  -vcf NA18526_CHROM1.vcf -log 
