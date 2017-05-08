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
These experiments do not show an advantage for focusing on the likelihood.

Initial experiments looking at the effect of focusing and guiding on parameter inferences were inconclusive -- I used
too little data (1 Mb).


------------------------------

=========
Overview
=========

.. report:: TrackerGeneral.Experiment
   :render: table
   :tracks: constpopsize_4epochs_guiding

   Overview of experimental parameters


Convergence of likelihood
=========================

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
   :tracker: name=constpopsize_4epochs_guiding,track=str_parameter,type=LogL,value=rate,selector="np=2000"
   :render: sb-box-plot
   :tracks: guide0.5_bias1.0_mstepFalse,guide0.5_bias2.0_mstepFalse
   :yrange: -24500,-23500
   :function: -23600
   :layout: row
   :width: 400

   Speed of convergence of recombination guide, without and with focusing, with 2000 particles

.. report:: TrackerConstPopSize.ConstpopsizeEMConvergence
   :tracker: name=constpopsize_4epochs_guiding,track=str_parameter,type=LogL,value=rate,selector="np=5000"
   :render: sb-box-plot
   :tracks: guide0.5_bias1.0_mstepFalse,guide0.5_bias2.0_mstepFalse
   :yrange: -24500,-23500
   :function: -23600
   :layout: row
   :width: 400

   Same with 5000 particles

.. report:: TrackerConstPopSize.ConstpopsizeEMConvergence
   :tracker: name=constpopsize_4epochs_guiding,track=str_parameter,type=LogL,value=rate,selector="np=10000"
   :render: sb-box-plot
   :tracks: guide0.5_bias1.0_mstepFalse,guide0.5_bias2.0_mstepFalse
   :yrange: -24500,-23500
   :function: -23600
   :layout: row
   :width: 400

   Same with 10000 particles


   
Convergence of parameter estimates
==================================


.. report:: TrackerConstPopSize.ConstpopsizeNe
   :tracker: name=constpopsize_4epochs_guiding,track=str_parameter
   :render: sb-box-plot
   :layout: row
   :function: 10000         
   :mpl-figure: tight_layout=True
   :width: 300
   :tracks: guide0.0_bias1.0_mstepTrue,guide0.0_bias2.0_mstepTrue,guide0.5_bias1.0_mstepTrue,guide0.5_bias2.0_mstepTrue
   :slices: T0      
   :groupby: none
   :yrange: 8000,10500

   Inference of population sizes, epoch T0
   

.. report:: TrackerConstPopSize.ConstpopsizeNe
   :tracker: name=constpopsize_4epochs_guiding,track=str_parameter
   :render: sb-box-plot
   :layout: row
   :function: 10000         
   :mpl-figure: tight_layout=True
   :width: 300
   :tracks: guide0.0_bias1.0_mstepTrue,guide0.0_bias2.0_mstepTrue,guide0.5_bias1.0_mstepTrue,guide0.5_bias2.0_mstepTrue
   :slices: T800
   :groupby: none
   :yrange: 8000,10500

   Inference of population sizes, epoch T800
   

.. report:: TrackerConstPopSize.ConstpopsizeNe
   :tracker: name=constpopsize_4epochs_guiding,track=str_parameter
   :render: sb-box-plot
   :layout: row
   :function: 10000         
   :mpl-figure: tight_layout=True
   :width: 300
   :tracks: guide0.0_bias1.0_mstepTrue,guide0.0_bias2.0_mstepTrue,guide0.5_bias1.0_mstepTrue,guide0.5_bias2.0_mstepTrue
   :slices: T4000
   :groupby: none
   :yrange: 8000,10500

   Inference of population sizes, epoch T4000
   

.. report:: TrackerConstPopSize.ConstpopsizeNe
   :tracker: name=constpopsize_4epochs_guiding,track=str_parameter
   :render: sb-box-plot
   :layout: row
   :function: 10000         
   :mpl-figure: tight_layout=True
   :width: 300
   :tracks: guide0.0_bias1.0_mstepTrue,guide0.0_bias2.0_mstepTrue,guide0.5_bias1.0_mstepTrue,guide0.5_bias2.0_mstepTrue
   :slices: T20000
   :groupby: none
   :yrange: 8000,10500

   Inference of population sizes, epoch T20000
   
      
