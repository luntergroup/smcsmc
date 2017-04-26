.. Test documentation master file, created by
   sphinxreport-quickstart 

**********************************************************
Effect of focusing and guiding on likelihood and Ne
**********************************************************

To answer the question whether focusing (biasing towards recent recombinations) and guiding (biasing towards recombinations
at loci and in subtrees that have posterior support) improves sampling (see `constpopsize_guiding.py`).
.. On bender with `-j 50`, the experiment takes a few hours.

Rationale
---------
   
Rationale: important but rare events such as recent recombinations sampled infrequently, so can be missed when they do occur;
focusing oversamples them at the expense of slightly less dense sampling of ordinary events.  Inferring recombinations 
ignores information inferred at earlier EM iterations; guiding oversamples recombination events at loci and subtrees that
had posterior support at those iterations.

Conclusion
----------
Guiding improves the likelihood, and is equivalent to using 2.5x to 4x more particles without guiding.
The improvement is iterative, showing that useful recombination events are remembered, as intended.
These experiments do not show an advantage for focusing on the likelihood.  TODO: an experiment to look at parameter inferences.


------------------------------

=========
Overview
=========

.. report:: TrackerGeneral.Experiment
   :render: table
   :tracks: constpopsize_4epochs_guiding

   Overview of experimental parameters


Convergence
===========

.. report:: TrackerConstPopSize.ConstpopsizeLikelihood
   :tracker: name=constpopsize_4epochs_guiding
   :render: sb-box-plot
   :slices: F1.0      
   :yrange: -24500,-23500
   :function: -23600
   :layout: row
   :width: 400
   :mpl-figure: tight_layout=True

   Likelihood without (G=0) or with guiding (G=0.5; 5 iterations; parameters fixed); no focusing.

.. report:: TrackerConstPopSize.ConstpopsizeLikelihood
   :tracker: name=constpopsize_4epochs_guiding
   :render: sb-box-plot
   :slices: F2.0
   :yrange: -24500,-23500
   :function: -23600
   :layout: row
   :width: 400       
   :mpl-figure: tight_layout=True

   As above; focusing with bias_strength=2.0

.. report:: TrackerConstPopSize.ConstpopsizeEMConvergence
   :tracker: name=constpopsize_4epochs_guiding,track=str_parameter,type=LogL,value=rate,selector="np=5000"
   :render: sb-box-plot
   :tracks: guide0.5_bias1.0_mstepFalse,guide0.5_bias2.0_mstepFalse
   :yrange: -24500,-23500
   :function: -23600
   :layout: row
   :width: 400

   Speed of convergence of recombination guide, without and with focusing, with 5000 particles

.. report:: TrackerConstPopSize.ConstpopsizeEMConvergence
   :tracker: name=constpopsize_4epochs_guiding,track=str_parameter,type=LogL,value=rate,selector="np=2000"
   :render: sb-box-plot
   :tracks: guide0.5_bias1.0_mstepFalse,guide0.5_bias2.0_mstepFalse
   :yrange: -24500,-23500
   :function: -23600
   :layout: row
   :width: 400

   Same with 2000 particles
   
