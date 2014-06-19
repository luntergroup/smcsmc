#!/bin/bash

foo(){
data_dir="/well/bsg/joezhu/experiment4_data/"
for rep in $(seq 1 1 15)
    do
./pop_struct.py ${case} ${nsam} ${rep}
mv ${case}Samples${nsam}msdata${rep}.vcf ${data_dir}
rm ${case}Samples${nsam}msdata${rep}*
    done
}

nsam=2

case="open"
foo

case="close"
foo

case="wakeley_a"
foo

case="wakeley_b"
foo

case="test-1"
foo


nsam=4

case="open"
foo

case="close"
foo

case="wakeley_a"
foo

case="wakeley_b"
foo

case="test-1"
foo
