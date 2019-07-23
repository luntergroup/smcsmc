#!/bin/bash

testObj=smcsmcSCRM
testObjDbg=smcsmcSCRM_dbg
source test_binaries.src

begin=1
end=1

#echo "Testing the program for different number of taxa  "
#for size in $(seq 2 7);do
 #test_smcsmc -nsam ${size} -Np 1000 -r 130 1000000 || exit 1
#done
#echo ""

#echo "Testing -EM "
#for size in $(seq 2 6);do
 #test_smcsmc -nsam ${size} -Np 100 -r 130 1000000 -EM 2 || exit 1
#done
#echo ""


function test_smcsmc_specific {
  recomb_seqlen="-r 130 1000000"
  #recomb_seqlen="-r 2 100"
  echo -n " ${testObj} -nsam $@ -Np 1 ${recomb_seqlen} "

  for i in `seq ${begin} ${end}`; do
    echo -n "."
    scrmCMD=" $@ 1 ${recomb_seqlen} -T -l 0 -seed ${i}"
    ./scrm ${scrmCMD} | grep ";" | sed -e "s/\[.*\]//g" > scrm_tmp_tree
    smcsmcCMD=" -nsam $@ -Np 1 ${recomb_seqlen} -seed ${i}"
    ./smcsmcSCRM ${smcsmcCMD} | grep ");" > smcsmc_tmp_tree
    diff <(head -1 scrm_tmp_tree) <(head -1 smcsmc_tmp_tree)
    if [ $? -ne 0 ]; then
      echo ""
      echo "Failed at random seed: " ${i}
      echo "Testing for initial tree in scrm and smcsmc."
      exit 1
    #else
      #echo ""
      #echo "Initial tree is ok"
    fi

    #echo ""

    diff <(tail -1 scrm_tmp_tree) <(tail -1 smcsmc_tmp_tree)
    if [ $? -ne 0 ]; then
      echo ""
      echo "Failed at random seed: " ${i}
      echo "Testing for Last tree in scrm and smcsmc."
      echo "./scrm ${scrmCMD}"
      echo "./smcsmcSCRM ${smcsmcCMD}"
      exit 1
    fi
  done
  echo ""
}

begin=1
end=10

echo "Testing smcsmcSCRM "
 for sample_size in $(seq 2 10); do
  test_smcsmc_specific ${sample_size}  || exit 1
 done
echo ""
