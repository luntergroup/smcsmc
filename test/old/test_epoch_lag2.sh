#!/bin/bash

source test_epoch_lag.src

begin=1
end=10000

echo "Testing epoch lagging "
 test_smcsmc -nsam 4 -Np 1 -r 130 1000000 || exit 1
echo ""
