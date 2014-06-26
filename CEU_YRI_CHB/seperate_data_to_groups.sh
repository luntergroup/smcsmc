#!/bin/bash

ORG_DATA_DIR="$HOME/Desktop/CEUYRICHB_reduced_data/"

# Afrian data YRI NA18501, NA18502
indv1="NA18501"
indv2="NA18502"

# European ancenstry (Utah) CEU NA06985, NA06994
indv3="NA06994"
indv4="NA12890"

# Chinese population CHB NA18525, NA18526
indv5="NA18525"
indv6="NA18526"

for chr in $(seq 1 1 1)
    do 
    VCFfile="${indv1}_${indv2}_${indv3}_${indv4}_${indv5}_${indv6}_CHROM"${chr}.vcf
    vcf-subset -c ${indv1},${indv2} ${ORG_DATA_DIR}${VCFfile} | awk '$10!="0|0"||$11!="0|0"' | awk '$10!="1|1"||$11!="1|1"' > ${ORG_DATA_DIR}"YRI_2sample_CHROM"${chr}.vcf
    vcf-subset -c ${indv3},${indv4} ${ORG_DATA_DIR}${VCFfile} | awk '$10!="0|0"||$11!="0|0"' | awk '$10!="1|1"||$11!="1|1"' > ${ORG_DATA_DIR}"CEU_2sample_CHROM"${chr}.vcf
    vcf-subset -c ${indv5},${indv6} ${ORG_DATA_DIR}${VCFfile} | awk '$10!="0|0"||$11!="0|0"' | awk '$10!="1|1"||$11!="1|1"' > ${ORG_DATA_DIR}"CHB_2sample_CHROM"${chr}.vcf
    done


