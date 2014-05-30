#!/bin/bash

ORG_DATA_DIR="$HOME/Desktop/1000genome_pase1_analysis_result_shapeit2/"
CEUYRI_DATA_DIR="$HOME/Desktop/CEUYRI_data/"
mkdir ${CEUYRI_DATA_DIR}

# Afrian data YRI NA18501, NA18502
indv1="NA18501"
indv2="NA18502"

# European ancenstry (Utah) CEU NA06985, NA06994
indv3="NA06994"
indv4="NA12890"

for chr in $(seq 1 1 22)
    do 
    VCFfile="${indv1}_${indv2}_${indv3}_${indv4}_CHROM"${chr}.vcf
    vcftools  --gzvcf ${ORG_DATA_DIR}ALL.chr${chr}.SHAPEIT2_integrated_phase1_v3.20101123.snps_indels_svs.genotypes.all.vcf.gz \
    --indv ${indv1} --indv ${indv2} --indv ${indv3} --indv ${indv4} \
    --recode --stdout > ${CEUYRI_DATA_DIR}${VCFfile}
    done
