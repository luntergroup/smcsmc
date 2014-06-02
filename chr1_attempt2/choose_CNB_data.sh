#!/bin/bash

ORG_DATA_DIR="$HOME/Desktop/1000genome_pase1_analysis_result_shapeit2/"
#
# Chinese data CHB NA18526
indv1="NA18526"

chr=1

VCFfile="${indv1}_CHROM"${chr}.vcf

vcftools  --gzvcf ${ORG_DATA_DIR}ALL.chr${chr}.SHAPEIT2_integrated_phase1_v3.20101123.snps_indels_svs.genotypes.all.vcf.gz \
--indv ${indv1} --recode --stdout | sed "/0\|0/d" | sed "/1\|1/d" > ${VCFfile}
