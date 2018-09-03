pth=/well/gerton/Tasks/smcsmc-analyses/run-singlepops

# first, run plot-zigzag.sh on bender to update .dta
scp rescomp:/well/gerton/Tasks/smcsmc-analyses/smcsmc-bender/performance/experiments/zigzag.dta .

for run in simb-P simb-PE simb-PL ; do
    dir=ceuyri4-${run}-vb
    mkdir -p ${dir}
    cat $pth/${dir}/emiter*/chunkfinal.out | awk '$13==-1 || NR==1' > ${dir}/result.out
done

mkdir -p ceuyri4-vb-dephase-30k
cp $pth/ceuyri4-vb-dephase-30k/result.out ceuyri4-vb-dephase-30k/result.out

for pop in ceu4 chb4 yri4 ; do
    echo $pop
    # 15 iterations, no vb
    mkdir -p ${pop}
    awk '$13==-1 || NR==1' $pth/${pop}.out > ${pop}/result.out
    # 15 iterations, with vb
    mkdir -p ${pop}-vb
    cp $pth/${pop}-vb/result.out ${pop}-vb/result.out
done

for pop in ceutsi4 ceugih4 ceuchb4 ceumxl4 ceuyri4 ; do
    echo $pop
    # 25 iterations, no vb
    mkdir -p ${pop}
    cp $pth/${pop}/result.out ${pop}/result.out
    # 30 iterations, with vb
    mkdir -p ${pop}-vb
    #cp $pth/${pop}-vb/result.out ${pop}-vb/result.out
    cat $pth/${pop}-vb/emiter*/chunkfinal.out | awk '$13==-1 || NR==1' > ${pop}-vb/result.out
done


