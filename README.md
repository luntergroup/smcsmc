# SMCSMC -- demographic inference using particle filters

## Installation
On Debian/Ubuntu based systems:
```bash
apt-get install build-essential autoconf autoconf-archive doxygen libcppunit-dev
./bootstrap
make
```

On Mac OS:
```bash
port install autoconf autoconf-archive doxygen cppunit
./bootstrap
make
```

## Experimental -- CGATReport

I had to install `virtualenv` and `libfreetype6-dev` before starting.  The packages `python-tk` and `libxft-dev` may also be required.
        
To install:        
```
        git submodule init
        git submodule update
        cd performance
        virtualenv --no-site-packages env
        source env/bin/activate
        pip install numpy
        pip install future
        pip install --upgrade setuptools
        cd CGATReport
        python setup.py install
```

To run, go to `performance/experiments` and run (one of) the scripts, to create the database.  Then,

```
        cd performance
        source env/bin/activate
        cd report
        make html
```
and point a browser to performance/report/_build/html/contents.html

                        
##Package status

package  | Build Status
-------- | -----------------
smcsmc   | [![Build Status](https://magnum.travis-ci.com/luntergroup/smcsmc.svg?branch=master)](https://magnum.travis-ci.com/luntergroup/smcsmc)
scrm     | [![Build Status](https://travis-ci.org/luntergroup/scrm_jz_stable_branch.svg?branch=smcsmcSCRM)](https://travis-ci.org/luntergroup/scrm_jz_stable_branch)



