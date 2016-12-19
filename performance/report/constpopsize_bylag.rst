.. Test documentation master file, created by
   sphinxreport-quickstart 

**********************************************************
Dependence of inferred Ne on lag
**********************************************************

Population size inference with variable lag (see `constpopsize_bylag_experiment.py`).

Rationale: increasing the lag increases the amount of information extracted from the data, and improves convergence.
However it also increases the variance of parameter estimates, so a balance should be struck.

In addition it is possible that the choice of lag influences the amount of bias in the estimates.

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
   :render: sb-box-plot
   :slices: T80000
   :yrange: 9000,10000
   :layout: row
   :width: 200       

   Speed of convergence for different lags (epoch T80000, fixed data; truth is not 10,000)

      
Variable data
=============

.. report:: TrackerConstPopSize.ConstpopsizeLagdependence
   :render: sb-box-plot
   :layout: row
   :function: 10000         
   :mpl-figure: tight_layout=True
   :width: 300
   :groupby: none
   :tracks: variableData       
   :yrange: 8000,10500

   Inference of a 3-epoch model on 10 Mb of sequence with variable lag around default of 4.0.

Conclusion.
   

Fixed data
===========

.. report:: TrackerConstPopSize.ConstpopsizeLagdependence
   :render: sb-box-plot
   :layout: row
   :function: 10000         
   :mpl-figure: tight_layout=True
   :width: 300
   :groupby: none
   :tracks: fixedData       
   :yrange: 8000,10500

   Inference of a 3-epoch model on 10 Mb of sequence with variable lag around default of 4.0.

Conclusion.
   

