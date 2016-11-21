#!/bin/bash

numProcess=$( grep -c ^processor /proc/cpuinfo )
maxTimeOut=1300

nosetests -v --processes=${numProcess} --process-timeout=${maxTimeOut}

# for specific test
# nosetests -v test_2sample.py
