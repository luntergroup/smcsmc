.. Test documentation master file, created by
   sphinxreport-quickstart 

**********************************************************
Dependence of inferred Ne on lag
**********************************************************

Population size inference with variable lag (see `constpopsize_bylag_experiment.py`).

=========
Overview
=========

.. report:: TrackerGeneral.Experiment
   :render: table
   :tracks: constpopsize_3epochs_lagdependence

   Overview of experimental parameters
                                 


Unique data
===========

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
   

