#!/bin/bash
rm recomb_rate_no_prune Ne_no_prune
for rep in $(seq 1 1 10000)
	do 
	pf-ARG -Np 1 -o trial_no_prune${rep} -seed ${rep}
	grep "inferred recomb rate = " trial_no_prune${rep}.log | sed -e "s/inferred recomb rate =//" >> recomb_rate_no_prune
	tail -1 trial_no_prune${rep}.log | sed -e "s/  (       0 ) |  //" >> Ne_no_prune
	rm trial_no_prune${rep}*
	done
