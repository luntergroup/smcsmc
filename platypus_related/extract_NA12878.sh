#!/bin/bash

BAMFILES="/well/1000G/data/NA12878/high_coverage_alignment/NA12878.mapped.ILLUMINA.bwa.CEU.high_coverage_pcr_free.20130906.bam"
REFFILE="human_g1k_v37.fasta"
OUTPUT=NA12878.vcf

python ~/Platypus_0.7.4/Platypus.py callVariants --output=${OUTPUT} --bamFiles=${BAMFILES} --regions=1 --nCPU=8 --genIndels=0 --refFile=${REFFILE} --outputRefCalls=1  
