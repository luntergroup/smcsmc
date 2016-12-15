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

Plots can be tested using `cgatreport-test`, for example like this:

```        
cgatreport-test --path=trackers -t TrackerConstPopSize.ConstpopsizeLengthdependence -r sb-box-plot -o yrange=4000,20000
```
                
To view the document, type `make html` and point a browser to `report/_build/html/contents.html`.
