#!/bin/bash

source test_binaries.src

begin=1
end=1

echo "Testing (Decide the case name later ... ) "
for size in $(seq 2 30);do 
 test_smcsmc -nsam ${size} -Np 1000 -r 130 1000000 || exit 1
done
echo ""

echo "Testing -EM "
for size in $(seq 2 30);do 
 test_smcsmc -nsam ${size} -Np 1000 -r 130 1000000 -EM 2 || exit 1
done
echo ""


