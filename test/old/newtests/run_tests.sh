#!/bin/bash

numProcess=$( grep -c ^processor /proc/cpuinfo )
maxTimeOut=1300

nosetests -v --processes=${numProcess} --process-timeout=${maxTimeOut} 2> nosetests.report
cat nosetests.report | ssh rescomp1 mail -s "bender_regression_test_report" joe.zhu@well.ox.ac.uk gerton.lunter@well.ox.ac.uk donnah@well.ox.ac.uk

# for specific test
# nosetests -v test_2sample.py
