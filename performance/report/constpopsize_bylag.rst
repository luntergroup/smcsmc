.. Test documentation master file, created by
   sphinxreport-quickstart 

**********************************************************
Dependence of inferred Ne on lag
**********************************************************

Population size inference with variable lag (see `constpopsize_bylag_experiment.py`).
On bender with `-j 50`, the experiment takes a few hours.

Rationale: increasing the lag increases the amount of information extracted from the data, and improves convergence.
However it also increases the variance of parameter estimates, and increases memory usage and runtime, so a balance should be struck.
In addition it is possible that the choice of lag influences the amount of bias in the estimates.


Conclusion from the experiment
------------------------------

A lag of 1.0 (or possibly 2.0) is sufficient.  However there is little evidence for increased variance of estimates, so the
downside of increasing the lag is not strong either -- except possibly for increased runtime and memory usage, which I didn't look at in this experiment.

Another observation from this experiment is that the bias in the inferences is more pronounced than I originally thought for this model, because I didn't let
previous experiments run for sufficient EM iterations.  The model doesn't seem to have completely converged even after 30 iterations.

=========
Overview
=========

.. report:: TrackerGeneral.Experiment
   :render: table
   :tracks: constpopsize_3epochs_lagdependence

   Overview of experimental parameters


Convergence
===========

.. report:: TrackerConstPopSize.ConstpopsizeEMConvergence
   :tracker: name=constpopsize_3epochs_lagdependence
   :render: sb-box-plot
   :slices: T80000
   :yrange: 6000,10000
   :layout: row
   :width: 200       

   Speed of convergence for different lags (epoch T80000, fixed data; truth is not 10,000)

      
Variable data
=============

.. report:: TrackerConstPopSize.ConstpopsizeNe
   :tracker: name=constpopsize_3epochs_lagdependence,column=lag
   :render: sb-box-plot
   :layout: row
   :function: 10000         
   :mpl-figure: tight_layout=True
   :width: 300
   :groupby: none
   :tracks: variableData       
   :yrange: 6000,11000

   Inference of a 3-epoch model on 50 Mb of sequence with variable lag around default of 4.0.


   

Fixed data
===========

.. report:: TrackerConstPopSize.ConstpopsizeNe
   :tracker: name=constpopsize_3epochs_lagdependence,column=lag
   :render: sb-box-plot
   :layout: row
   :function: 10000         
   :mpl-figure: tight_layout=True
   :width: 300
   :groupby: none
   :tracks: fixedData       
   :yrange: 6000,11000

   Inference of a 3-epoch model on 50 Mb of sequence with variable lag around default of 4.0.


   

