#!/bin/bash

function test_smcsmc {
  echo -n " smcsmc $@ "
  for i in `seq 1 10`; do
    echo -n "."

    # Test using smcsmc self-checks
    ./smcsmc_dbg $@ -seed $i > /dev/null 
    if [ $? -ne 0 ]; then
      echo ""
      echo "Executing \"./smcsmc_dbg $@ -seed $i\" failed."
      echo "Debug Call: make -mj2 smcsmc_dbg && ./smcsmc_dbg $@ -seed $i 2>&1 | less"
      exit 1
    fi

    # Test for memory leaks
    valgrind --error-exitcode=1 --leak-check=full -q ./smcsmc $@ -seed $i > /dev/null
    if [ $? -ne 0 ]; then
      echo ""
      echo "Valgrind check of \"./smcsmc $@ -seed $i\" failed."
      exit 1
    fi

	# Test for lagging utility
	./smcsmc $@ -seed $i -o tmp > /dev/null
	./smcsmc $@ -seed $i -o tmplag -lag 10000 > /dev/null
	diff tmpNe tmplagNe
	if [ $? -ne 0 ]; then
      echo ""
      echo "Lagging utility is not working, no data with/without lag failed."
      exit 1
    fi
  done
  echo " done."
}

echo "Testing (Decide the case name later ... ) "
 test_smcsmc -nsam 2 -Np 1000 -r 130 1000000 || exit 1
 test_smcsmc -nsam 3 -Np 1000 -r 130 1000000 || exit 1
 test_smcsmc -nsam 4 -Np 1000 -r 130 1000000 || exit 1
 test_smcsmc -nsam 5 -Np 1000 -r 130 1000000 || exit 1
echo ""

echo "Testing -EM "
 test_smcsmc -nsam 2 -Np 1000 -r 130 1000000 -EM 2 || exit 1
 test_smcsmc -nsam 3 -Np 1000 -r 130 1000000 -EM 2 || exit 1
 test_smcsmc -nsam 4 -Np 1000 -r 130 1000000 -EM 2 || exit 1
 test_smcsmc -nsam 5 -Np 1000 -r 130 1000000 -EM 2 || exit 1
echo ""
