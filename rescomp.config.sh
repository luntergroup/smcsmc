# to make things compile and run on rescomp

module unload python
module unload gcc
module load gcc/5.4.0
module load python/2.7.10-gcc5.4.0

export BOOST_ROOT=/apps/well/boost/1.59.0-gcc5.4.0-py27
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/users/lunter/gerton/usrlocal/lib/
export PATH=/apps/well/python/2.7.10/bin:$PATH

# this is to work around SGE not passing through LD_LIBRARY_PATH
export LLP=$LD_LIBRARY_PATH
