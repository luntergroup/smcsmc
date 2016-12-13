# Performance reports

## Overview

A [CGATReport](https://www.cgat.org/downloads/public/CGATReport) document summarizing
the performance of `smcsmc`

## Setup

Install using
```
pip install CGATReport
pip install future
sudo apt-get install python-tk
```

Experiments are defined by Python scripts in `experiments/`.  When run, these store
their results in a mysql database `experiments/experimentsdb`.

These results are accessed by CGATReport using scripts in `trackers/`.  The document
itself is `contents.rst`.

To view the document, type `make html` and point a browser to `_build/html/contents.html`.


