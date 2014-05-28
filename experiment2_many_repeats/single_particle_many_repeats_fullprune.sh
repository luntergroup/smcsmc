#!/bin/bash
rm recomb_rate_full_prune Ne_full_prune
for rep in $(seq 1 1 10000)
	do 
	pf-ARG -Np 1 -o trial_full_prune${rep} -seed ${rep}
	grep "inferred recomb rate = " trial_full_prune${rep}.log | sed -e "s/inferred recomb rate =//" >> recomb_rate_full_prune
	tail -1 trial_full_prune${rep}.log | sed -e "s/  (       0 ) |  //" >> Ne_full_prune
	rm trial_full_prune${rep}*
	done
