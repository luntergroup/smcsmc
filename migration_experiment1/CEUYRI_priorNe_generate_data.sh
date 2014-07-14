#!/bin/bash
source CEUYRI_priorNe_param.src

mkdir ${dir}
positionfile=${fileprefix}position
segfile=${fileprefix}seg

scrm ${nsam} 1 ${cmd} ${assign_haploid} ${split} ${mig_pattern} ${pop_struct} -l 100000 -seed 1 > ${fileprefix}

grep 'positions' ${fileprefix} | sed -e 's/positions: //' >  ${positionfile}
tail -${nsam} ${fileprefix} > ${segfile}

echo "./ms2something.py vcf ${seqlen} ${positionfile} ${segfile} ${fileprefix}"
./ms2something.py vcf ${seqlen} ${positionfile} ${segfile} ${fileprefix}
 
