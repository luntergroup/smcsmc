#!/bin/bash

function test_smcsmc {
  echo -n " smcsmc $@ "
  for i in `seq ${begin} ${end}`; do
    echo -n "."

    ## Test using smcsmc self-checks
    #smcsmc_dbg $@ -seed $i > /dev/null 
    #if [ $? -ne 0 ]; then
      #echo ""
      #echo "Executing \"smcsmc_dbg $@ -seed $i\" failed."
      #echo "Debug Call: make -mj2 smcsmc_dbg && smcsmc_dbg $@ -seed $i 2>&1 | less"
      #exit 1
    #fi

	# Test using smcsmc self-checks
    ../smcsmcRECOMBOFF_dbg $@ -seed $i > /dev/null 
    if [ $? -ne 0 ]; then
      echo ""
      echo "Executing \"smcsmcRECOMBOFF_dbg $@ -seed $i\" failed."
      echo "Debug Call: make -mj2 smcsmcRECOMBOFF_dbg && ../smcsmcRECOMBOFF_dbg $@ -seed $i 2>&1 | less"
      exit 1
    fi

    # Test for memory leaks
    valgrind --error-exitcode=1 --leak-check=full -q ../smcsmcRECOMBOFF $@ -seed $i > /dev/null
    if [ $? -ne 0 ]; then
      echo ""
      echo "Valgrind check of \"smcsmcRECOMBOFF $@ -seed $i\" failed."
      exit 1
    fi
	
    # Test for updating recomb rate
    ../smcsmcRECOMBOFF $@ -seed $i -EM 2 | grep "recomb rate " | sed -e "s/^.*= //g" -e "s/ //g" > rate_tmp
    diff <(head -1 rate_tmp) <(tail -1 rate_tmp)
	  if [ $? -ne 0 ]; then
      echo ""
      echo "Failed at random seed: " ${i}
      echo "Testing for smcsmcRECOMBOFF updating recomb rate."
      exit 1
    fi
    
    # Test for single iteration smcsmc vs smcsmcRECOMBOFF
    ../smcsmcRECOMBOFF $@ -seed $i -o recomboff > /dev/null
    ../smcsmc          $@ -seed $i -o recombon  > /dev/null
    diff <(tail -n+4  recomboff.log | sed -e "/inferred/d") <(tail -n+4  recombon.log | sed -e "/inferred/d")
    if [ $? -ne 0 ]; then
      echo ""
      echo "Failed at random seed: " ${i}
      echo "Testing for smcsmcRECOMBOFF vs smcsmc."
      exit 1
    fi
    
  done
  echo " done."
}

begin=1
end=2

echo "Testing smcsmcRECOMBOFF "
 test_smcsmc -nsam 2 -Np 1000 -r 130 1000000 || exit 1
 test_smcsmc -nsam 3 -Np 1000 -r 130 1000000 || exit 1
 test_smcsmc -nsam 4 -Np 1000 -r 130 1000000 || exit 1
 test_smcsmc -nsam 5 -Np 1000 -r 130 1000000 || exit 1
echo ""
