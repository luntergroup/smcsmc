.. Test documentation master file, created by
   sphinxreport-quickstart 

**********************************************************
Dependence of inferred Ne on amount of number of particles
**********************************************************

Population size inference with varying numbers of particles.
These runs used 10 Mb of sequence, using 10 replicates, 5 EM iterations, and a constant population size of 10000.
(see `constpopsize_byparticle_experiment.py`)

=========
Overview
=========

.. report:: TrackerGeneral.Experiment
   :render: table
   :tracks: constpopsize_3epochs_particledependence

   Overview of experimental parameters
                                 


Unique data
===========

.. report:: TrackerConstPopSize.ConstpopsizeParticledependence
   :render: sb-box-plot
   :layout: row
   :function: 10000         
   :mpl-figure: tight_layout=True
   :width: 300
   :groupby: none
   :yrange: 8000,10500

   Inference of a 3-epoch model on 10 Mb of sequence with varying numbers of particles.

The bias is clearly due to undersampling, and increasing the number of particles mitigates the issue.
   

------------

.. report:: TrackerConstPopSize.ConstpopsizeParticledependence
   :render: table
         
   Raw data table

