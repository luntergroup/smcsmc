.. Test documentation master file, created by
   sphinxreport-quickstart 

**********************************************************
Dependence of inferred Ne on ancestry information
**********************************************************

Population size inference with variable ancestry information (see `twopops_unidirmigr_ancestralaware_experiment.py`).

Rationale: Providing ancestral information should increase power, particularly for inferring migration. 
However a strict prior of 0=ancestral, 1=derived puts more pressure on the sampler,
as more trees will be completely incompatible with the data.


Conclusion from the experiment
------------------------------



=========
Overview
=========

.. report:: TrackerGeneral.Experiment
   :render: table
   :tracks: twopops_unidirmigr_falsestart_ancestralawaredependence

   Overview of experimental parameters


4samples by Number of particles
======================

Here we are comparing no ancestry information vs known ancestry after the first iteration with 4 samples

.. report:: TrackerUniDirMigr.TwopopsUnidirmigrAncestralawaredependence_bynp
   :tracker: name=twopops_unidirmigr_falsestart_ancestralawaredependence,column=ancestral_aware
   :render: sb-box-plot
   :layout: row
   :function: 10000
   :mpl-figure: tight_layout=True
   :width: 300
   :groupby: none
   :tracks: Np100
   :yrange: 9000,17000

.. report:: TrackerUniDirMigr.TwopopsUnidirmigrAncestralawaredependence_bynp
   :tracker: name=twopops_unidirmigr_falsestart_ancestralawaredependence,column=ancestral_aware
   :render: sb-box-plot
   :layout: row
   :function: 10000
   :mpl-figure: tight_layout=True
   :width: 300
   :groupby: none
   :tracks: Np500
   :yrange: 9000,17000

.. report:: TrackerUniDirMigr.TwopopsUnidirmigrAncestralawaredependence_bynp
   :tracker: name=twopops_unidirmigr_falsestart_ancestralawaredependence,column=ancestral_aware
   :render: sb-box-plot
   :layout: row
   :function: 10000
   :mpl-figure: tight_layout=True
   :width: 300
   :groupby: none
   :tracks: Np1000
   :yrange: 9000,17000




8 samples by Number of particles
======================

Here we are comparing aa=False and aa=True after the first iteration with 8 samples

.. report:: TrackerUniDirMigr.TwopopsUnidirmigrAncestralawaredependence_bynp
   :tracker: name=twopops_unidirmigr_falsestart_ancestralawaredependence,column=ancestral_aware
   :render: sb-box-plot
   :layout: row
   :function: 10000
   :mpl-figure: tight_layout=True
   :width: 300
   :groupby: none
   :tracks: Np100_8s
   :yrange: 9000,17000

.. report:: TrackerUniDirMigr.TwopopsUnidirmigrAncestralawaredependence_bynp
   :tracker: name=twopops_unidirmigr_falsestart_ancestralawaredependence,column=ancestral_aware
   :render: sb-box-plot
   :layout: row
   :function: 10000
   :mpl-figure: tight_layout=True
   :width: 300
   :groupby: none
   :tracks: Np500_8s
   :yrange: 9000,17000

.. report:: TrackerUniDirMigr.TwopopsUnidirmigrAncestralawaredependence_bynp
   :tracker: name=twopops_unidirmigr_falsestart_ancestralawaredependence,column=ancestral_aware
   :render: sb-box-plot
   :layout: row
   :function: 10000
   :mpl-figure: tight_layout=True
   :width: 300
   :groupby: none
   :tracks: Np1000_8s
   :yrange: 9000,17000


