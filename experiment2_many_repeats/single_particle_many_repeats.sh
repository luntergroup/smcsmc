#!/bin/bash
rm recomb_rate Ne
for rep in $(seq 1 1 100000)
	do 
	pf-ARG -Np 1 -o trial${rep} -seed ${rep}
	grep "inferred recomb rate = " trial${rep}.log | sed -e "s/inferred recomb rate =//" >> recomb_rate
	tail -1 trial${rep}.log | sed -e "s/  (       0 ) |  //" >> Ne
	rm trial${rep}*
	done
