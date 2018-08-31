# note: also edit Makefile.am when compiling on rescomp

# tell configure which gcc to use:
module unload python
module unload gcc
module load gcc/5.4.0
module load python/2.7.10-gcc5.4.0

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/apps/well/boost/1.59.0-gcc5.4.0-py27/lib
#export PATH=/apps/well/python/2.7.10/bin:$PATH

# work around SGE not passing through LD_LIBRARY_PATH
export LLP=$LD_LIBRARY_PATH

