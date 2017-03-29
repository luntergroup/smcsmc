.. Test documentation master file, created by
   sphinxreport-quickstart 

***********************************************************************
Dependence of inferred Ne on lag fraction when NOT initialised at truth
***********************************************************************

Population size inference with variable lag fraction when Ne is not initialised at the truth (see `constpopsize_bylag_experiment.py`).


Conclusion from the experiment
------------------------------



=========
Overview
=========

.. report:: TrackerGeneral.Experiment
   :render: table
   :tracks: constpopsize_4epochs_falsestart_lagdependence

   Overview of experimental parameters


Results by Number of samples
============================

Here we are comparing lag fractions, using inferred pop size after the first iteration

.. report:: TrackerConstPopSize.ConstpopsizeLagdependenceFalsestart
   :tracker: name=constpopsize_4epochs_falsestart_lagdependence,column=lag
   :render: sb-box-plot
   :layout: row
   :function: 10000
   :mpl-figure: tight_layout=True
   :width: 300
   :groupby: none
   :tracks: 2s
   :yrange: 9000,17000

.. report:: TrackerConstPopSize.ConstpopsizeLagdependenceFalsestart
   :tracker: name=constpopsize_4epochs_falsestart_lagdependence,column=lag
   :render: sb-box-plot
   :layout: row
   :function: 10000
   :mpl-figure: tight_layout=True
   :width: 300
   :groupby: none
   :tracks: 4s
   :yrange: 9000,17000

.. report:: TrackerConstPopSize.ConstpopsizeLagdependenceFalsestart
   :tracker: name=constpopsize_4epochs_falsestart_lagdependence,column=lag
   :render: sb-box-plot
   :layout: row
   :function: 10000
   :mpl-figure: tight_layout=True
   :width: 300
   :groupby: none
   :tracks: 8s
   :yrange: 9000,17000

One EM step is not very informative, let's look at the inference after 30 iterations below

Results for last EM step by Number of samples
=============================================

Here we are comparing lag fractions, using inferred pop size after 30 EM iterations.

.. report:: TrackerConstPopSize.ConstpopsizeLagdependenceFalsestart_lastiter
   :tracker: name=constpopsize_4epochs_falsestart_lagdependence,column=lag
   :render: sb-box-plot
   :layout: row
   :function: 10000
   :mpl-figure: tight_layout=True
   :width: 300
   :groupby: none
   :tracks: 2s
   :yrange: 9000,17000

.. report:: TrackerConstPopSize.ConstpopsizeLagdependenceFalsestart_lastiter
   :tracker: name=constpopsize_4epochs_falsestart_lagdependence,column=lag
   :render: sb-box-plot
   :layout: row
   :function: 10000
   :mpl-figure: tight_layout=True
   :width: 300
   :groupby: none
   :tracks: 4s
   :yrange: 9000,17000

.. report:: TrackerConstPopSize.ConstpopsizeLagdependenceFalsestart_lastiter
   :tracker: name=constpopsize_4epochs_falsestart_lagdependence,column=lag
   :render: sb-box-plot
   :layout: row
   :function: 10000
   :mpl-figure: tight_layout=True
   :width: 300
   :groupby: none
   :tracks: 8s
   :yrange: 9000,17000

It is clear that a non-zero lag is necessary.
As for optimal choice of lag fraction, this seems to vary depending on the number of samples.
Note we are only using 500 particles without any pilots, so 4 and 8 samples are incredibly challenging.
I am currently running this analysis with 3k particles.

Note this is the last experiment run with the bug placing the SNP at the beginning of the segment.
The Np=3k version of this experiment has the fix, but these Np=500 results will need to be redone.
