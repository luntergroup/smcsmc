.. Test documentation master file, created by
   sphinxreport-quickstart 

*************************************************************
Dependence of inferred Ne on number of particles -- revisited
*************************************************************

Population size inference with varying number of particles (see `constpopsize_byparticles2_experiment.py`).

Rationale: The previous experiment investigating the dependence on lag showed that earlier experiments used
too few EM iterations, and increasing the EM iterations showed a stronger bias than was observed before.
This experiment investigates whether increasing the number of particles mitigates this bias.

=========
Overview
=========

.. report:: TrackerGeneral.Experiment
   :render: table
   :tracks: constpopsize_3epochs_particles2,constpopsize_3epochs_particles3

   Overview of experimental parameters


Convergence
===========

.. report:: TrackerConstPopSize.ConstpopsizeEMConvergence2
   :render: sb-box-plot
   :slices: T0
   :yrange: 10000,11200
   :layout: row
   :width: 200
   :function: 10000
              
   Speed of convergence for different numbers of particles (epoch T0)


.. report:: TrackerConstPopSize.ConstpopsizeEMConvergence2
   :render: sb-box-plot
   :slices: T40000
   :yrange: 8000,11000
   :layout: row
   :width: 200
   :function: 10000
              
   Speed of convergence for different numbers of particles (epoch T40000)


   
.. report:: TrackerConstPopSize.ConstpopsizeEMConvergence2
   :render: sb-box-plot
   :slices: T80000
   :yrange: 6000,10000
   :layout: row
   :width: 200
   :function: 10000
              
   Speed of convergence for different numbers of particles (epoch T80000)
   
      
Variable data
=============

.. report:: TrackerConstPopSize.ConstpopsizeParticles2dependence
   :render: sb-box-plot
   :layout: row
   :function: 10000         
   :mpl-figure: tight_layout=True
   :width: 300
   :groupby: none
   :tracks: variableData
   :yrange: 6000,11000

   Inference of a 3-epoch model on 50 Mb of sequence with variable numbers of particles

Variable data - no pilot bias
=============================

.. report:: TrackerConstPopSize.ConstpopsizeParticles3dependence
   :render: sb-box-plot
   :layout: row
   :function: 10000         
   :mpl-figure: tight_layout=True
   :width: 300
   :groupby: none
   :tracks: variableData
   :yrange: 6000,11000

   Inference of a 3-epoch model on 50 Mb of sequence with variable numbers of particles.  This run
   does not use pilot biasing towards early epochs.




