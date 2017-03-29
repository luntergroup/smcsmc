.. Test documentation master file, created by
   sphinxreport-quickstart 

**********************************************************
Dependence of inferred Ne on phasing information
**********************************************************

Population size inference with variable phasing (see `constpopsize_phasedependence_experiment.py`).


Conclusion from the experiment
------------------------------

As seen before we get better results when using unphased data. 
However both cases are diverging from the truth. Bug or 4 samples too challenging for Np?

I don't know why the multiple iterations Np=5k plots are not being generated correctly...

=========
Overview
=========

.. report:: TrackerGeneral.Experiment
   :render: table
   :tracks: constpopsize_4epochs_pahsedependence

   Overview of experimental parameters


4samples by Number of particles
======================

Here we are comparing unphased vs phased after the first iteration with 4 samples

.. report:: TrackerConstPopSize.ConstpopsizePhasedependence_bynp
   :tracker: name=constpopsize_4epochs_pahsedependence,column=phased
   :render: sb-box-plot
   :layout: row
   :function: 10000
   :mpl-figure: tight_layout=True
   :width: 300
   :groupby: none
   :tracks: Np100
   :yrange: 9000,17000

.. report:: TrackerConstPopSize.ConstpopsizePhasedependence_bynp
   :tracker: name=constpopsize_4epochs_pahsedependence,column=phased
   :render: sb-box-plot
   :layout: row
   :function: 10000
   :mpl-figure: tight_layout=True
   :width: 300
   :groupby: none
   :tracks: Np500
   :yrange: 9000,17000

.. report:: TrackerConstPopSize.ConstpopsizePhasedependence_bynp
   :tracker: name=constpopsize_4epochs_pahsedependence,column=phased
   :render: sb-box-plot
   :layout: row
   :function: 10000
   :mpl-figure: tight_layout=True
   :width: 300
   :groupby: none
   :tracks: Np1000
   :yrange: 9000,17000

.. report:: TrackerConstPopSize.ConstpopsizePhasedependence_bynp
   :tracker: name=constpopsize_4epochs_pahsedependence,column=phased
   :render: sb-box-plot
   :layout: row
   :function: 10000
   :mpl-figure: tight_layout=True
   :width: 300
   :groupby: none
   :tracks: Np5000
   :yrange: 9000,17000



8 samples by Number of particles
================================

Here we are comparing unphased and phased after the first iteration with 8 samples

.. report:: TrackerConstPopSize.ConstpopsizePhasedependence_bynp
   :tracker: name=constpopsize_4epochs_pahsedependence,column=phased
   :render: sb-box-plot
   :layout: row
   :function: 10000
   :mpl-figure: tight_layout=True
   :width: 300
   :groupby: none
   :tracks: Np100_8s
   :yrange: 9000,17000

.. report:: TrackerConstPopSize.ConstpopsizePhasedependence_bynp
   :tracker: name=constpopsize_4epochs_pahsedependence,column=phased
   :render: sb-box-plot
   :layout: row
   :function: 10000
   :mpl-figure: tight_layout=True
   :width: 300
   :groupby: none
   :tracks: Np500_8s
   :yrange: 9000,17000

.. report:: TrackerConstPopSize.ConstpopsizePhasedependence_bynp
   :tracker: name=constpopsize_4epochs_pahsedependence,column=phased
   :render: sb-box-plot
   :layout: row
   :function: 10000
   :mpl-figure: tight_layout=True
   :width: 300
   :groupby: none
   :tracks: Np1000_8s
   :yrange: 9000,17000

.. report:: TrackerConstPopSize.ConstpopsizePhasedependence_bynp
   :tracker: name=constpopsize_4epochs_pahsedependence,column=phased
   :render: sb-box-plot
   :layout: row
   :function: 10000
   :mpl-figure: tight_layout=True
   :width: 300
   :groupby: none
   :tracks: Np5000_8s
   :yrange: 9000,17000


4samples by Number of particles, show multiple iterations
=========================================================

Now look at the first and last iteration, to check the longterm behaviour

.. report:: TrackerConstPopSize.ConstpopsizePhasedependence_bynp_multiters
   :tracker: name=constpopsize_4epochs_pahsedependence,column=phased
   :render: sb-box-plot
   :layout: row
   :function: 10000
   :mpl-figure: tight_layout=True
   :width: 300
   :groupby: none
   :tracks: Np100
   :yrange: 9000,17000

.. report:: TrackerConstPopSize.ConstpopsizePhasedependence_bynp_multiters
   :tracker: name=constpopsize_4epochs_pahsedependence,column=phased
   :render: sb-box-plot
   :layout: row
   :function: 10000
   :mpl-figure: tight_layout=True
   :width: 300
   :groupby: none
   :tracks: Np500
   :yrange: 9000,17000

.. report:: TrackerConstPopSize.ConstpopsizePhasedependence_bynp_multiters
   :tracker: name=constpopsize_4epochs_pahsedependence,column=phased
   :render: sb-box-plot
   :layout: row
   :function: 10000
   :mpl-figure: tight_layout=True
   :width: 300
   :groupby: none
   :tracks: Np1000
   :yrange: 9000,17000

.. report:: TrackerConstPopSize.ConstpopsizePhasedependence_bynp_multiters
   :tracker: name=constpopsize_4epochs_pahsedependence,column=phased
   :render: sb-box-plot
   :layout: row
   :function: 10000
   :mpl-figure: tight_layout=True
   :width: 300
   :groupby: none
   :tracks: Np5000
   :yrange: 9000,17000

It appears that the advantage of not phasing is due to slower divergence from the truth, but it is still diverging. Bug?


8samples by Number of particles, show multiple iterations
=========================================================

Now look at the first and last iteration, to check the longterm behaviour

.. report:: TrackerConstPopSize.ConstpopsizePhasedependence_bynp_multiters
   :tracker: name=constpopsize_4epochs_pahsedependence,column=phased
   :render: sb-box-plot
   :layout: row
   :function: 10000
   :mpl-figure: tight_layout=True
   :width: 300
   :groupby: none
   :tracks: Np100_8s
   :yrange: 9000,17000

.. report:: TrackerConstPopSize.ConstpopsizePhasedependence_bynp_multiters
   :tracker: name=constpopsize_4epochs_pahsedependence,column=phased
   :render: sb-box-plot
   :layout: row
   :function: 10000
   :mpl-figure: tight_layout=True
   :width: 300
   :groupby: none
   :tracks: Np500_8s
   :yrange: 9000,17000

.. report:: TrackerConstPopSize.ConstpopsizePhasedependence_bynp_multiters
   :tracker: name=constpopsize_4epochs_pahsedependence,column=phased
   :render: sb-box-plot
   :layout: row
   :function: 10000
   :mpl-figure: tight_layout=True
   :width: 300
   :groupby: none
   :tracks: Np1000_8s
   :yrange: 9000,17000

.. report:: TrackerConstPopSize.ConstpopsizePhasedependence_bynp_multiters
   :tracker: name=constpopsize_4epochs_pahsedependence,column=phased
   :render: sb-box-plot
   :layout: row
   :function: 10000
   :mpl-figure: tight_layout=True
   :width: 300
   :groupby: none
   :tracks: Np5000_8s
   :yrange: 9000,17000
