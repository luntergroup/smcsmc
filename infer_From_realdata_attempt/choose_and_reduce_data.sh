#!/bin/bash
vcf-subset -c NA18501 NA18501_NA18502_NA06994_NA12890_CHROM1.vcf | sed "/0\|0/d" | sed "/1\|1/d" > NA18501_CHROM1.vcf
