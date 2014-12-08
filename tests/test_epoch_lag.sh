#!/bin/bash

source test_epoch_lag.src

begin=1
end=10

echo "Testing epoch lagging "
 test_smcsmc -nsam 2 -Np 1000 -r 130 1000000 || exit 1
 test_smcsmc -nsam 3 -Np 1000 -r 130 1000000 || exit 1
 test_smcsmc -nsam 4 -Np 1000 -r 130 1000000 || exit 1
 test_smcsmc -nsam 5 -Np 1000 -r 130 1000000 || exit 1
echo ""
