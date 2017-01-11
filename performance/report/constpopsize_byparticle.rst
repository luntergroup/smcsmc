.. Test documentation master file, created by
   sphinxreport-quickstart 

**********************************************************
Dependence of inferred Ne on number of particles
**********************************************************

Population size inference with varying numbers of particles.
(see `constpopsize_byparticle_experiment.py`).
Runtime on bender with `-j 50` is about 15 minutes.

=========
Overview
=========

.. report:: TrackerGeneral.Experiment
   :render: table
   :tracks: constpopsize_3epochs_particledependence

   Overview of experimental parameters
                                 


Variable data
==============

.. report:: TrackerConstPopSize.ConstpopsizeNe
   :tracker: name=constpopsize_3epochs_particledependence
   :render: sb-box-plot
   :layout: row
   :function: 10000         
   :mpl-figure: tight_layout=True
   :width: 300
   :tracks: variableData     
   :groupby: none
   :yrange: 8000,10500

   Inference of a 3-epoch model on 50 Mb of sequence with varying numbers of particles.

The bias is clearly due to undersampling, and increasing the number of particles mitigates the issue.



Fixed data
=============

.. report:: TrackerConstPopSize.ConstpopsizeNe
   :tracker: name=constpopsize_3epochs_particledependence
   :render: sb-box-plot
   :layout: row
   :function: 10000         
   :mpl-figure: tight_layout=True
   :width: 300
   :tracks: fixedData
   :groupby: none
   :yrange: 8000,10500

   Inference of a 3-epoch model on 50 Mb of sequence with varying numbers of particles, on fixed input data.




