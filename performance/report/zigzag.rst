.. Test documentation master file, created by
   sphinxreport-quickstart 

**********************************************************************
Effect of focusing and guiding on likelihood and Ne - zigzag model
**********************************************************************

Investigate focusing and guiding on the bottleneck model, and particularly, investigate how well
inference works (see `bottleneck_inference_guiding.py`).

Rationale
---------
   
Conclusion
----------


=========
Overview
=========

.. report:: TrackerGeneral.Experiment
   :render: table
   :tracks: zigzag_inference_guiding

   Overview of experimental parameters


Convergence of likelihood
=========================

.. report:: TrackerConstPopSize.ConstpopsizeLikelihood
   :tracker: name=zigzag_inference_guiding,mstep=True
   :render: sb-box-plot
   :slices: F1.0
   :yrange: -2077000,-2040000
   :function: -2040700
   :layout: row
   :width: 400
   :mpl-figure: tight_layout=True

   Likelihood without (G=0) or with guiding (G=0.5; 5 iterations; parameters fixed); no focusing.

.. report:: TrackerConstPopSize.ConstpopsizeLikelihood
   :tracker: name=zigzag_inference_guiding,mstep=True
   :render: sb-box-plot
   :slices: F2.0
   :yrange: -2077000,-2040000
   :function: -2040700
   :layout: row
   :width: 400       
   :mpl-figure: tight_layout=True

   As above; focusing with bias_strength=2.0


Convergence of Ne estimates
===========================

.. note truth=scale,Ne(0),t1,a1,t2,a2,...,tn (series of exponential growth/declines)

.. report:: TrackerInferredNe.PlotNe
   :tracker: name=zigzag_inference_guiding,minne=1000,maxne=100000,mint=15,maxt=50000,bias=1.0,guide=0.0,truth="14312;5;0.000582262;1318.18;0.00232905;-329.546;0.00931619;82.3865;0.0372648;-20.5966;0.149059;5.14916;0.596236"
   :render: user
   :layout: row
   
   Inferred Ne values (unbiased, unguided)

.. report:: TrackerInferredNe.PlotNe
   :tracker: name=zigzag_inference_guiding,minne=1000,maxne=100000,mint=15,maxt=50000,bias=2.0,guide=0.0,truth="14312;5;0.000582262;1318.18;0.00232905;-329.546;0.00931619;82.3865;0.0372648;-20.5966;0.149059;5.14916;0.596236"
   :render: user
   :layout: row
   
   Inferred Ne values (bias 2.0, unguided)

.. report:: TrackerInferredNe.PlotNe
   :tracker: name=zigzag_inference_guiding,minne=1000,maxne=100000,mint=15,maxt=50000,bias=2.0,guide=0.0,showconv=True,truth="14312;5;0.000582262;1318.18;0.00232905;-329.546;0.00931619;82.3865;0.0372648;-20.5966;0.149059;5.14916;0.596236"
   :render: user
   :layout: row
   
   Convergence of Ne (bias 2.0, unguided)
            
