Introduction
============


SMCSMC, short for **S**\ equential **M**\ onte **C**\ arlo for the **S**\ equentially **M**\ arkovian **C**\ oalescent, is a method for estimating ancestral population size and migration history from sequence data. It has several advantages over comparable methods, especially when inferring asymetric migration models. 

The method simulates many trees (:term:`particles<Particle>`) from a coalescent model and updates them using information about genetic variation along the sequence. Stochastic expectation maximization is used to update demographic parameters (population sizes, migration matrix) in iterative steps until convergence is reached. 

SMCSMC takes as input optionally phased sequencing data formatted as segments, and provides utilities for data conversion from common formats.


Installation
-----------

We *highly* recommend installing SMCSMC from :code:`conda`, as it comes packaged with all necessary dependencies. A seperate guide for manual compilation may be found `in the developer reference <https::github.com>`_\ .  `See here <https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html>`_ for a helpful guide to installing and using :code:`conda` to manage programs. 

First add both :code:`conda-forge` and :code:`terhorst` to your channel lists (if they are not there already), then install :code:`smcsmc`. 

.. code-block:: bash

   conda --add channel conda-forge
   conda --add channel terhorst

   conda install smcsmc


Basic Usage
----------

To use SMCSMC, start a :code:`python` session and import the :code:`smcsmc` module. As a part of the installation above, two binaries are installed into the :code:`conda-bin`, `smcsmc` (inference) and `scrm` (simulation). The front end, :code:`smcsmc` is a wrapper around these binaries providing convenient functions for data manipulation, conversion, plotting, and utilities surrounding the workflow for analysing sequences with SMCSMC. With a test seg file such as `this one <https://github.com>`_\ , the following will run a default session of SMCSMC.

.. code-block:: python

        import smcsmc

        test_args = {
                `seg`: `test_seg.seg`,
                `nsam`: 4
        }

        smcsmc.utils.run_smcsmc(test_args)
                        

Follow the Getting Started guide to become familiar with the basic structure and function of SMCSMC commands, then look at one of the tutorials for analysing simulated or real data. 

Other Methods
------------

SMCSMC is part of the `PopSim consortium <https://github.com/popgensims>`_\ , and we are actively involved in building a framework to standardize population genetic analyses. Part of this involves making it easy to run the same analysis with many different methods. We have built :code:`smcsmc` with this goal in mind. For the latest information about comparisons between different population genetic software, including :code:`smc++`, :code:`stairwayplot`, :code:`msmc`, and :code:`dadi/fastcoal`, check out the `PopSim analysis repository <https://github.com/popgensims/analysis>`_\ .

.. figure:: ../img/popsim.png
   :scale: 50 %
   :align: center
   
   Population history of a European-acting individual inferred from five replicates of the :code:`stdpopsim.homo_sapiens.GutenkunstThreePopOutOfAfrica` model of human history.


Citation
--------

If you use :code:`SMC2` in your work, please cite the following article:

        1. Henderson, D., Zhu, S. (Joe), & Lunter, G. (2018). Demographic inference using particle filters for continuous Markov jump processes. BioRxiv, 382218. https://doi.org/10.1101/382218
        2. Staab, P. R., Zhu, S., Metzler, D., & Lunter, G. (2015). scrm: efficiently simulating long sequences using the approximated coalescent with recombination. Bioinformatics, 31(10), 1680–1682. https://doi.org/10.1093/bioinformatics/btu861 


