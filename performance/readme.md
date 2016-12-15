# Performance reports

## Overview

A [CGATReport](https://www.cgat.org/downloads/public/CGATReport) document summarizing
the performance of `smcsmc`

For installation notes see the top-level readme.
        
## Usage

Experiments are defined by Python scripts in `experiments/`.  When run, these store
their results in a mysql database `experiments/experimentsdb`.

These results are accessed by CGATReport using scripts in `report/trackers/`.  The
document itself is `contents.rst`.  It uses [ReStructuredText](http://docutils.sourceforge.net/docs/user/rst/quickref.html)
See the [CGATReport Tutorials](https://www.cgat.org/downloads/public/CGATReport/documentation/Tutorials.html) for various
examples.

Before you use `cgatreport`, be sure to activate the virtual environment, by going to the `performance` directory and
running `source env/bin/activate` (virtualenv) or `source activate env` (conda/anaconda).  (Use `deactivate` or `source deactivate`
to exit a virtual python environment.)
        
Plots can be tested using `cgatreport-test`, for example like this:

```
cd report        
cgatreport-test --path=trackers -t TrackerConstPopSize.ConstpopsizeLengthdependence -r sb-box-plot -o yrange=4000,20000
```
                
To view the full document, type `make html` in the `report` directory and point a browser to `report/_build/html/contents.html`.

        
## Installation -- Linux

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

## Installation -- Mac

This worked for me.  I had [Anaconda](https://docs.continuum.io/anaconda/install) and [Homebrew](http://brew.sh/) installed.

```
        brew install freetype
        git submodule init
        git submodule update
        cd performance
        conda create --name env python
        source activate
        conda install numpy future setuptools
        cd CGATReport
        python setup.py install
        
