#!/bin/bash

function test_smcsmc {
  recomb_seqlen="-r 130 1000000"
  echo -n " smcsmcSCRM -nsam $@ ${recomb_seqlen} "
  
  for i in `seq ${begin} ${end}`; do
    echo -n "."
    ../scrm $@ 1 ${recomb_seqlen} -T -l 0 -seed ${i} | grep ";" | sed -e "s/\[.*\]//g" > scrm_tmp_tree
    ../smcsmcSCRM -nsam $@ -Np 1 ${recomb_seqlen} -seed ${i} | grep ";" > smcsmc_tmp_tree
    diff <(head -1 scrm_tmp_tree) <(head -1 smcsmc_tmp_tree)
	  if [ $? -ne 0 ]; then
      echo ""
      echo "Failed at random seed: " ${i}
      echo "Testing for initial tree in scrm and smcsmc."
      exit 1
    #else 
      #echo "first is ok"
    fi
    
    diff <(tail -1 scrm_tmp_tree) <(tail -1 smcsmc_tmp_tree)
	  if [ $? -ne 0 ]; then
      echo ""
      echo "Failed at random seed: " ${i}
      echo "Testing for Last tree in scrm and smcsmc."
      exit 1
    fi
  done
  echo ""
}


begin=1
end=100

echo "Testing smcsmcRECOMBOFF "
 for sample_size in $(seq 2 10); do 
  test_smcsmc ${sample_size}  || exit 1
 done
echo ""
