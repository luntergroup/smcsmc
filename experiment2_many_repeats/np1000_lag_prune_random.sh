#!/bin/bash
#$ -cwd
#$ -V
#$ -P bsg.prjb -q short.qb
#$ -e ErrFiles
#$ -o OutFiles
#$ -N np1000_lag_prune_random
#$ -j y

replicates=100
EMsteps=20
Np=1000
lag=100000
pruning=50000
ESS=1

source nodata_pfARG_singlepop.src
