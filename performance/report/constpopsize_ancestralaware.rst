.. Test documentation master file, created by
   sphinxreport-quickstart 

**********************************************************
Dependence of inferred Ne on ancestry information
**********************************************************

Population size inference with variable ancestry information (see `constpopsize_ancestralaware_experiment.py`).
Here we are looking at a single population with 8 samples.

Rationale: Providing ancestral information should increase power (particularly for inferring migration, which
is not tested here). However a strict prior of 0=ancestral, 1=derived puts more pressure on the sampler,
as more trees will be completely incompatible with the data.


Conclusion from the experiment
------------------------------



=========
Overview
=========

.. report:: TrackerGeneral.Experiment
   :render: table
   :tracks: constpopsize_4epochs_ancestralawaredependence

   Overview of experimental parameters



Variable data
=============

.. report:: TrackerConstPopSize.ConstpopsizeNe
   :tracker: name=constpopsize_4epochs_ancestralawaredependence,column=ancestral_aware
   :render: sb-box-plot
   :layout: row
   :function: 10000         
   :mpl-figure: tight_layout=True
   :width: 300
   :groupby: none
   :tracks: variableData       
   :yrange: 6000,11000

   Inference of a 4-epoch model on 50 Mb of sequence with variable ancestry information.



Fixed data
===========

.. report:: TrackerConstPopSize.ConstpopsizeNe
   :tracker: name=constpopsize_4epochs_ancestralawaredependence,column=ancestral_aware
   :render: sb-box-plot
   :layout: row
   :function: 10000         
   :mpl-figure: tight_layout=True
   :width: 300
   :groupby: none
   :tracks: fixedData       
   :yrange: 6000,11000

   Inference of a 4-epoch model on 50 Mb of sequence with variable ancestry information.
