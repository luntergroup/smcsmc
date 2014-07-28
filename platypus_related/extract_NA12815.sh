#!/bin/bash
BAMFILES="/well/1000G/data/NA12815/alignment/NA12815.mapped.ILLUMINA.bwa.CEU.low_coverage.20130415.bam"
REFFILE="/well/1000G/reference/human_g1k_v37.fasta.gz"
OUTPUT=NA12815.vcf

python ~/Platypus_0.7.4/Platypus.py callVariants --output=${OUTPUT} --bamFiles=${BAMFILES} --regions=1 --nCPU=4 --genIndels=0 --refFile=${REFFILE}
