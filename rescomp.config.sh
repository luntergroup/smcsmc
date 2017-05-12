# note: also edit Makefile.am when compiling on rescomp

# tell configure which gcc to use:
export CXX=/apps/well/gcc/5.4.0/bin/gcc

# before execution:
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/users/mcvean/joezhu/gcc-4.8.3/lib64
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/apps/well/boost/1.55.0/lib
export PATH=/apps/well/python/2.7.10/bin:$PATH

# python packages
pip install --user sqlalchemy
