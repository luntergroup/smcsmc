#!/bin/bash
rm recomb_rate Ne
for rep in $(seq 1 1 10000)
	do 
	pf-ARG -l 0 -Np 1 -o trial${rep} -seed ${rep}
	grep "inferred recomb rate = " trial${rep}.log | sed -e "s/inferred recomb rate =//" >> recomb_rate_pruned
	tail -1 trial${rep}.log | sed -e "s/  (       0 ) |  //" >> Ne_pruned
	rm trial${rep}*
	done
