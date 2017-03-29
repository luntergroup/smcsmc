.. Test documentation master file, created by
   sphinxreport-quickstart 

**********************************************************
Dependence of inferred Ne on lag fraction
**********************************************************

Population size inference with variable lag fraction (see `constpopsize_bylag_truestart_experiment.py`).


Conclusion from the experiment
------------------------------

This experiment shows that a non-zero lag is necessary to avoid divergence from the true value.
This does not conclusively show any lag fraction as optimal. In fact the divergence in the first
epoch is worrying for all choices. Is there a bug in the code? Running with a lot of particles
for the 2 sample case should help answer this question. If there is a bug, the fact that lag=0
does not diverge in one EM step indicates the bug is not in the event bookkeeping.

A similar experiment where Ne is not initialised at the true value is reported in the next section.


=========
Overview
=========

.. report:: TrackerGeneral.Experiment
   :render: table
   :tracks: constpopsize_4epochs_lagdependence

   Overview of experimental parameters


Results by Number of samples
============================

Here we are comparing lag fractions, using inferred pop size after the first iteration

.. report:: TrackerConstPopSize.ConstpopsizeLagdependence
   :tracker: name=constpopsize_4epochs_lagdependence,column=lag
   :render: sb-box-plot
   :layout: row
   :function: 10000
   :mpl-figure: tight_layout=True
   :width: 300
   :groupby: none
   :tracks: 2s
   :yrange: 9000,17000

.. report:: TrackerConstPopSize.ConstpopsizeLagdependence
   :tracker: name=constpopsize_4epochs_lagdependence,column=lag
   :render: sb-box-plot
   :layout: row
   :function: 10000
   :mpl-figure: tight_layout=True
   :width: 300
   :groupby: none
   :tracks: 4s
   :yrange: 9000,17000

.. report:: TrackerConstPopSize.ConstpopsizeLagdependence
   :tracker: name=constpopsize_4epochs_lagdependence,column=lag
   :render: sb-box-plot
   :layout: row
   :function: 10000
   :mpl-figure: tight_layout=True
   :width: 300
   :groupby: none
   :tracks: 8s
   :yrange: 9000,17000

Lag fraction = 0 does the best as without lag the events are driven solely by the sampler, which is the truth for this first iteration.
There is no observable difference in the non-zero lags for one iteration started at the truth.


Results for last EM step by Number of samples
=============================================

Here we are comparing lag fractions, using inferred pop size after 30 EM iterations.
This should show which lags can recover from a slight move away from the true value.

.. report:: TrackerConstPopSize.ConstpopsizeLagdependence_lastiter
   :tracker: name=constpopsize_4epochs_lagdependence,column=lag
   :render: sb-box-plot
   :layout: row
   :function: 10000
   :mpl-figure: tight_layout=True
   :width: 300
   :groupby: none
   :tracks: 2s
   :yrange: 9000,17000

.. report:: TrackerConstPopSize.ConstpopsizeLagdependence_lastiter
   :tracker: name=constpopsize_4epochs_lagdependence,column=lag
   :render: sb-box-plot
   :layout: row
   :function: 10000
   :mpl-figure: tight_layout=True
   :width: 300
   :groupby: none
   :tracks: 4s
   :yrange: 9000,17000

.. report:: TrackerConstPopSize.ConstpopsizeLagdependence_lastiter
   :tracker: name=constpopsize_4epochs_lagdependence,column=lag
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
I am running the false start version of this experiment with both 500 and 3k particles to address this.
