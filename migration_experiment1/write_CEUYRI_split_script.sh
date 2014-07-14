#!/bin/bash

for TASK in $(seq 1 1 15)
    do 
    echo "#!/bin/bash
source CEUYRI_priorNe_param.src

data=\${fileprefix}".vcf"
echo \${data}
Np=1000
EM=15

prefix=CEUYRI_priorNe_\${nsam}sampleNp\${Np}

pattern=\"1*3+25*2+1*4+1*6\"
tmax=8

infer_mig_pattern=\"-eM 0 1 -eM 0.030252 1 \"

outprefix=\${prefix}_Splitrep${TASK}
pf-ARG -l 0 \${cmd} \${assign_haploid} \${split} \${infer_mig_pattern} \${pop_struct} -EM \${EM} -Np \${Np} -o \${outprefix} -vcf \${data} -seed ${TASK} ${TASK} ${TASK} " > CEUYRI_split_script${TASK}.sh
    chmod a+x CEUYRI_split_script${TASK}.sh
    done
# should add pattern to this!
