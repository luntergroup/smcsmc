.. Test documentation master file, created by
   sphinxreport-quickstart 

**********************************************************************
Effect of focusing and guiding on likelihood and Ne - bottleneck model
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
   :tracks: bottleneck_inference_guiding

   Overview of experimental parameters


Convergence of likelihood
=========================

.. report:: TrackerConstPopSize.ConstpopsizeLikelihood
   :tracker: name=bottleneck_inference_guiding,mstep=True
   :render: sb-box-plot
   :slices: F1.0
   :yrange: -380000,-372000
   :function: -372050
   :layout: row
   :width: 400
   :mpl-figure: tight_layout=True

   Likelihood without (G=0) or with guiding (G=0.5; 5 iterations; parameters fixed); no focusing.

.. report:: TrackerConstPopSize.ConstpopsizeLikelihood
   :tracker: name=bottleneck_inference_guiding,mstep=True
   :render: sb-box-plot
   :slices: F2.0
   :yrange: -380000,-372000
   :function: -372050
   :layout: row
   :width: 400       
   :mpl-figure: tight_layout=True

   As above; focusing with bias_strength=2.0


Convergence of Ne estimates
===========================

.. note truth=scale,Ne(0),t1,Ne(t1),t2,...,tN,Ne(tN)

.. report:: TrackerInferredNe.PlotNe
   :tracker: name=bottleneck_inference_guiding,minne=500,maxne=100000,mint=30,maxt=100000,bias=1.0,guide=0.0,truth="10000;1;0.01;0.1;0.06;1;0.2;0.5;1;1;2;2"
   :render: user
   :layout: row
   
   Inferred Ne values (unbiased, unguided)

.. report:: TrackerInferredNe.PlotNe
   :tracker: name=bottleneck_inference_guiding,minne=500,maxne=100000,mint=30,maxt=100000,bias=2.0,guide=0.0,truth="10000;1;0.01;0.1;0.06;1;0.2;0.5;1;1;2;2"
   :render: user
   :layout: row
   
   Inferred Ne values (bias 2.0, unguided)

.. report:: TrackerInferredNe.PlotNe
   :tracker: name=bottleneck_inference_guiding,minne=500,maxne=100000,mint=30,maxt=100000,bias=2.0,guide=0.0,showconv=True,truth="10000;1;0.01;0.1;0.06;1;0.2;0.5;1;1;2;2"
   :render: user
   :layout: row
   
   Convergence of Ne (bias 2.0, unguided)
            
