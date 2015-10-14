#!/bin/bash
# Timeout at 24 individuals
#smcsmc -nsam 23 -Np 1000 -r 130 1000000 . done.
#smcsmc -nsam 24 -Np 1000 -r 130 1000000 .command cd tests; ./test_binary_sample_size.sh took more than 3 minutes since last output

source test_binaries.src

begin=1
end=10

echo "Testing the program for different number of taxa  "
for size in $(seq 2 7);do
 test_smcsmc -nsam ${size} -Np 1000 -r 130 1000000 || exit 1
done
echo ""

echo "Testing -EM "
for size in $(seq 2 6);do
 test_smcsmc -nsam ${size} -Np 100 -r 130 1000000 -EM 2 || exit 1
done
echo ""


