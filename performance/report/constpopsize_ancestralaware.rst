.. Test documentation master file, created by
   sphinxreport-quickstart 

**********************************************************
Dependence of inferred Ne on ancestry information
**********************************************************

Population size inference with variable ancestry information (see `constpopsize_ancestralaware_experiment.py`).

Rationale: Providing ancestral information should increase power (particularly for inferring migration, which
is not tested here). However a strict prior of 0=ancestral, 1=derived puts more pressure on the sampler,
as more trees will be completely incompatible with the data.


Conclusion from the experiment
------------------------------

Known ancestry seems to provide slight improvements in a single population model. These slight improvements may be
due to more frequent resampling leading to slower divergence away from the initial true values; this should be tested
by initialising parameters elsewhere.

I now need to test if known ancestry improves migration estimates in a two population model. In particular, known ancestry should
help to detect the direction of migration.


=========
Overview
=========

.. report:: TrackerGeneral.Experiment
   :render: table
   :tracks: constpopsize_4epochs_ancestralawaredependence

   Overview of experimental parameters


4samples by Number of particles
===============================

Here we are comparing no ancestry information vs known ancestry after the first iteration with 4 samples

.. report:: TrackerConstPopSize.ConstpopsizeAncestralawaredependence_bynp
   :tracker: name=constpopsize_4epochs_ancestralawaredependence,column=ancestral_aware
   :render: sb-box-plot
   :layout: row
   :function: 10000
   :mpl-figure: tight_layout=True
   :width: 300
   :groupby: none
   :tracks: Np100
   :yrange: 9000,17000

.. report:: TrackerConstPopSize.ConstpopsizeAncestralawaredependence_bynp
   :tracker: name=constpopsize_4epochs_ancestralawaredependence,column=ancestral_aware
   :render: sb-box-plot
   :layout: row
   :function: 10000
   :mpl-figure: tight_layout=True
   :width: 300
   :groupby: none
   :tracks: Np500
   :yrange: 9000,17000

.. report:: TrackerConstPopSize.ConstpopsizeAncestralawaredependence_bynp
   :tracker: name=constpopsize_4epochs_ancestralawaredependence,column=ancestral_aware
   :render: sb-box-plot
   :layout: row
   :function: 10000
   :mpl-figure: tight_layout=True
   :width: 300
   :groupby: none
   :tracks: Np1000
   :yrange: 9000,17000

.. report:: TrackerConstPopSize.ConstpopsizeAncestralawaredependence_bynp
   :tracker: name=constpopsize_4epochs_ancestralawaredependence,column=ancestral_aware
   :render: sb-box-plot
   :layout: row
   :function: 10000
   :mpl-figure: tight_layout=True
   :width: 300
   :groupby: none
   :tracks: Np5000
   :yrange: 9000,17000

Surprisingly, for Np=100,500,1000 AA=True does more badly (relative to AA=False) for the experiments with more particles. I had expected that
AA=True would increase power but rely more on the sampler, and hence would do better in runs with more particles.
The additional experiment with Np=5000 shows my expectation holds for a more realistic Np. I should also try this
with biased sampling and eventually directed sampling, as those help our sampler too. 


8 samples by Number of particles
================================

Here we are comparing aa=False and aa=True after the first iteration with 8 samples

.. report:: TrackerConstPopSize.ConstpopsizeAncestralawaredependence_bynp
   :tracker: name=constpopsize_4epochs_ancestralawaredependence,column=ancestral_aware
   :render: sb-box-plot
   :layout: row
   :function: 10000
   :mpl-figure: tight_layout=True
   :width: 300
   :groupby: none
   :tracks: Np100_8s
   :yrange: 9000,17000

.. report:: TrackerConstPopSize.ConstpopsizeAncestralawaredependence_bynp
   :tracker: name=constpopsize_4epochs_ancestralawaredependence,column=ancestral_aware
   :render: sb-box-plot
   :layout: row
   :function: 10000
   :mpl-figure: tight_layout=True
   :width: 300
   :groupby: none
   :tracks: Np500_8s
   :yrange: 9000,17000

.. report:: TrackerConstPopSize.ConstpopsizeAncestralawaredependence_bynp
   :tracker: name=constpopsize_4epochs_ancestralawaredependence,column=ancestral_aware
   :render: sb-box-plot
   :layout: row
   :function: 10000
   :mpl-figure: tight_layout=True
   :width: 300
   :groupby: none
   :tracks: Np1000_8s
   :yrange: 9000,17000

.. report:: TrackerConstPopSize.ConstpopsizeAncestralawaredependence_bynp
   :tracker: name=constpopsize_4epochs_ancestralawaredependence,column=ancestral_aware
   :render: sb-box-plot
   :layout: row
   :function: 10000
   :mpl-figure: tight_layout=True
   :width: 300
   :groupby: none
   :tracks: Np5000_8s
   :yrange: 9000,17000

Strange result. When we have 8 samples, I expected the sampler to struggle when AA=True with few
particles. But it seems that AA=True does as well, if not better, than AA=False for 8 samples even with 100 particles. Could this be a result of
AA=True diverging more slowly than AA=False due to frequent resampling? To test this hypothesis, I would need
to intialize Ne at a different value or print the resampling file.


4samples by Number of particles, show multiple iterations
=========================================================

Now we'll look at the first iteration and the last iteration, to check the longterm behaviour

.. report:: TrackerConstPopSize.ConstpopsizeAncestralawaredependence_bynp_multiters
   :tracker: name=constpopsize_4epochs_ancestralawaredependence,column=ancestral_aware
   :render: sb-box-plot
   :layout: row
   :function: 10000
   :mpl-figure: tight_layout=True
   :width: 300
   :groupby: none
   :tracks: Np100
   :yrange: 9000,17000

.. report:: TrackerConstPopSize.ConstpopsizeAncestralawaredependence_bynp_multiters
   :tracker: name=constpopsize_4epochs_ancestralawaredependence,column=ancestral_aware
   :render: sb-box-plot
   :layout: row
   :function: 10000
   :mpl-figure: tight_layout=True
   :width: 300
   :groupby: none
   :tracks: Np500
   :yrange: 9000,17000

.. report:: TrackerConstPopSize.ConstpopsizeAncestralawaredependence_bynp_multiters
   :tracker: name=constpopsize_4epochs_ancestralawaredependence,column=ancestral_aware
   :render: sb-box-plot
   :layout: row
   :function: 10000
   :mpl-figure: tight_layout=True
   :width: 300
   :groupby: none
   :tracks: Np1000
   :yrange: 9000,17000

.. report:: TrackerConstPopSize.ConstpopsizeAncestralawaredependence_bynp_multiters
   :tracker: name=constpopsize_4epochs_ancestralawaredependence,column=ancestral_aware
   :render: sb-box-plot
   :layout: row
   :function: 10000
   :mpl-figure: tight_layout=True
   :width: 300
   :groupby: none
   :tracks: Np5000
   :yrange: 9000,17000

8samples by Number of particles, show multiple iterations
=========================================================

Now we'll look at the first iteration and the last iteration, to check the longterm behaviour

.. report:: TrackerConstPopSize.ConstpopsizeAncestralawaredependence_bynp_multiters
   :tracker: name=constpopsize_4epochs_ancestralawaredependence,column=ancestral_aware
   :render: sb-box-plot
   :layout: row
   :function: 10000
   :mpl-figure: tight_layout=True
   :width: 300
   :groupby: none
   :tracks: Np100_8s
   :yrange: 9000,17000

.. report:: TrackerConstPopSize.ConstpopsizeAncestralawaredependence_bynp_multiters
   :tracker: name=constpopsize_4epochs_ancestralawaredependence,column=ancestral_aware
   :render: sb-box-plot
   :layout: row
   :function: 10000
   :mpl-figure: tight_layout=True
   :width: 300
   :groupby: none
   :tracks: Np500_8s
   :yrange: 9000,17000

.. report:: TrackerConstPopSize.ConstpopsizeAncestralawaredependence_bynp_multiters
   :tracker: name=constpopsize_4epochs_ancestralawaredependence,column=ancestral_aware
   :render: sb-box-plot
   :layout: row
   :function: 10000
   :mpl-figure: tight_layout=True
   :width: 300
   :groupby: none
   :tracks: Np1000_8s
   :yrange: 9000,17000

.. report:: TrackerConstPopSize.ConstpopsizeAncestralawaredependence_bynp_multiters
   :tracker: name=constpopsize_4epochs_ancestralawaredependence,column=ancestral_aware
   :render: sb-box-plot
   :layout: row
   :function: 10000
   :mpl-figure: tight_layout=True
   :width: 300
   :groupby: none
   :tracks: Np5000_8s
   :yrange: 9000,17000

