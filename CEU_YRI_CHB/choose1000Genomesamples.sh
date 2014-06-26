#!/bin/bash

ORG_DATA_DIR="$HOME/Desktop/1000genome_pase1_analysis_result_shapeit2/"
CEUYRICHB_DATA_DIR="$HOME/Desktop/CEUYRICHB_data/"
mkdir ${CEUYRICHB_DATA_DIR}

# Afrian data YRI NA18501, NA18502
indv1="NA18501"
indv2="NA18502"

# European ancenstry (Utah) CEU NA12890, NA06994
indv3="NA06994"
indv4="NA12890"

# Chinese population CHB NA18525, NA18526
indv5="NA18525"
indv6="NA18526"

for chr in $(seq 1 1 1)
    do 
    VCFfile="${indv1}_${indv2}_${indv3}_${indv4}_${indv5}_${indv6}_CHROM"${chr}.vcf
    vcftools  --gzvcf ${ORG_DATA_DIR}ALL.chr${chr}.SHAPEIT2_integrated_phase1_v3.20101123.snps_indels_svs.genotypes.all.vcf.gz \
    --indv ${indv1} --indv ${indv2} --indv ${indv3} --indv ${indv4} --indv ${indv5} --indv ${indv6} \
    --recode --stdout > ${CEUYRICHB_DATA_DIR}${VCFfile}
    done
