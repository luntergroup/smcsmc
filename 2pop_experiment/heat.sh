#!/bin/bash

namefiles(){
    positionfile=${fileprefix}"position"
    segfile=${fileprefix}"seg"
    tree_file=${fileprefix}"msTrees"
    tmrca_file=${fileprefix}"mstmrca"
    mschange_file=${fileprefix}"mschange"
}

msprocess(){
    grep 'positions' ${fileprefix} | sed -e "s/positions: //" >  ${positionfile}
    tail -4 ${fileprefix} > ${segfile}
    echo "./ms2something.py vcf ${seqlen} ${positionfile} ${segfile} ${fileprefix}"
    ./ms2something.py vcf ${seqlen} ${positionfile} ${segfile} ${fileprefix}
    grep ";" ${fileprefix} | sed -e "s/\[.*\]//g" >  ${tree_file}     
    hybrid-Lambda -gt ${tree_file} -tmrca ${tmrca_file}
    grep ";" ${fileprefix} | sed -e "s/\[//g" | sed -e "s/\].*;//g" >  ${mschange_file}
}

prefix=Heat2popMerge
rm ${prefix}*
theta=180
rho=30
seqlen=200001
for rep in $(seq 1 1 5)
    do 
    fileprefix=${prefix}${rep}
    namefiles
    mscmd=" 4 1 -t ${theta} -r ${rho} ${seqlen} -I 2 2 2 -ej 2.0 2 1 -seed ${rep} ${rep} ${rep} -p 10 -T "
    ms ${mscmd} > ${fileprefix}
    msprocess
    pf-ARG  -Np 100000 -l 0 -t ${theta} -r ${rho} ${seqlen} -I 2 2 2 -ej 2.0 2 1 -seed ${rep} -vcf ${fileprefix}.vcf -o ${fileprefix} -heat
    ./generate_heatmap.py ${fileprefix} ${seqlen}
    done



