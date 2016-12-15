.. Test documentation master file, created by
   sphinxreport-quickstart 

=========
 Contents
=========

Population size inference for varying amounts of input data

.. report:: TrackerConstPopSize.ConstpopsizeLengthdependence
   :render: sb-box-plot
   :function: 10000
   :yrange: 4000,20000
   :track: missingData
   :width: 200

   When there is no data, there is no particle degeneracy, and inferred Ne values have very little variance.

.. report:: TrackerConstPopSize.ConstpopsizeLengthdependence
   :render: sb-box-plot
   :function: 10000
   :yrange: 4000,20000
   :track: fixedData
   :width: 200

   Results from running inference over different segments of a fixed input sequence.
   The inferred Ne values have now a considerable variance, because of particle degeneracy.
   The means differ for the different length because of stochasticity in the coalescent process
   that was used to simulate the data.

.. report:: TrackerConstPopSize.ConstpopsizeLengthdependence
   :render: sb-box-plot
   :function: 10000
   :yrange: 4000,20000
   :track: fixedData
   :width: 200

   Resuls from running inference on different segments of different sequences, all simulaed
   under the same model.  The variance is now even higher due to the stochasticity
   of the coalescent process now contributing to the variance of the estimates.

.. report:: TrackerConstPopSize.ConstpopsizeLengthdependence
   :render: table            

   Table of results
            
==================
Indices and Tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


