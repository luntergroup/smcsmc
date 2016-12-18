.. Test documentation master file, created by
   sphinxreport-quickstart 

*************************************************
Dependence of inferred Ne on amount of input data
*************************************************

Population size inference for varying amounts of input data.
These runs are with 100 particles, 10 replicates, 5 EM iterations, and a constant population size of 10000.
(see `constpopsize_bylength_experiment.py`)

=========
Overview
=========

.. report:: TrackerGeneral.Experiment
   :render: table
   :tracks: constpopsize_3epochs_lengthdependence

   Overview of experimental parameters
                                 

No data
=======

.. report:: TrackerConstPopSize.ConstpopsizeLengthdependence
   :render: sb-box-plot
   :layout: row      
   :function: 10000
   :yrange: 4000,20000
   :tracks: missingData
   :groupby: none      
   :width: 300
   :mpl-figure: tight_layout=True

   Inference results in the absence of data.
                
..  groupby none causes 3 plots to be created, one for each slice (epoch)
..  mpl-figure passes the option 'tight_layout=True' to the matplotlib figure command, stopping x-labels from being chopped off
                
When there is no data, there is no particle degeneracy, and inferred Ne values have no bias and very little variance.

Fixed data
==========

.. report:: TrackerConstPopSize.ConstpopsizeLengthdependence
   :render: sb-box-plot
   :layout: row      
   :function: 10000
   :yrange: 4000,20000
   :tracks: fixedData
   :groupby: none      
   :width: 300
   :mpl-figure: tight_layout=True

   Results from running inference over different segments of a fixed input sequence.
   The blue line indicates the true population size.

When population sizes are inferred from data, the values have more variance because of particle degeneracy.
The means differ for the different length because of stochasticity in the coalescent process
that was used to simulate the data.


Variable data
=============

.. report:: TrackerConstPopSize.ConstpopsizeLengthdependence
   :render: sb-box-plot
   :layout: row      
   :function: 10000
   :yrange: 4000,20000
   :tracks: variableData
   :groupby: none      
   :width: 300
   :mpl-figure: tight_layout=True

   Resuls from running inference on different segments of different sequences, all simulated
   from the same model.

The variances in Ne estimates are now much higher due to the stochasticity
of the coalescent process, which dominates the variance due to the inference process itself.
This variance decreases with increasing amounts of input data.

It is also clear that except for the 0-40000 epoch, there is substantial bias in the estimates,
particularly for the >80000 epoch.  This bias does not appear to reduce with increasing amounts
of input data.

------------

.. report:: TrackerConstPopSize.ConstpopsizeLengthdependence
   :render: table
         
   Raw data table

