# note: also edit Makefile.am when compiling on rescomp

# tell configure which gcc to use:
module purge
module load gcc/5.4.0

# optimizations
#export CXXMARCH="-march nehalem"

# before execution:
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/apps/well/boost/1.59.0-gcc5.4.0-py27/lib
export PATH=/apps/well/python/2.7.10/bin:$PATH

# this is to work around SGE not passing through LD_LIBRARY_PATH
export LLP=$LD_LIBRARY_PATH

