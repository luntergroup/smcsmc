Getting Started
===============

This guide will familiarise your with the concepts necessary to run analyses in SMCSMC and take you through the basic aspects of an analysis using  toy data. You may also be interested in one of several tutorials using real and simulated data below. 

.. note::
        This tutorial is intended to be run on a personal computer, and values have been scaled accordingly. For analysis on real data, we expect that a user will have access to a compute cluster, and we provide guidance about setting up :code:`smcsmc` to function with common architechtures.  

Basic Concepts
--------------

SMCSMC is a particle filter, which means that a given number of :term:`particles<Particle>` are simulated, evaluated for an approximate likelihood, and resampled until demographic parameters have converged. Many of the options relate to the behaviour of this particle filter, and the remainder deal with the demographic model that you wish to infer. SMCSMC can infer complex demographic models involving several populations, though you should be weary of overspecifying your model. Run time increases drastically with the number of parameters needed to specify the model.   

Input Format
+++++++++++

Suppose we have the following data, saved into a file named :code:`toy.seg`. 

.. csv-table::
   :header: segment start, segment length, invariant, break, chr, genotype
   :widths: 5, 5, 5, 5, 5, 5

    1,       521,     T,       F,       1,       01\.0
    522,     2721,    T,       F,       1,       0111
    3243,    1758,    T,       F,       1,       10\/\/
    5001,    1296,    T,       F,       1,       0000
    6297,    1,       T,       F,       1,       \.\.\.\.
    6298,    4669,    T,       F,       1,       0110
    10967,   880,     T,       T,       1,       0100
    1,       708,     T,       F,       2,       1010


This gives the state of four haplotypes along the sequence, indicating segment starts and their lengths. If a segment is the last in a block before coordinates reset (i.e. when including multiple chromosomes in a single :code:`seg` file, it is called a genetic break. Genotypes in terms of either major/minor or ancestral/derived are given for each haplotype. Missing data is coded with a period (:code:`.`) whilst unphased variants are given with a forward slash (:code:`/`). 

.. note::
        This data format is similar to the input for `msmc <https://github.com/stschiff/msmc>`_ however, in the second column we require the length from the current segment to the next, rather than the number of called sites from the previous segment to the current one. 


Segment files are specified by the :code:`seg` flag, or if more than one are used, by the :code:`segs` flag. Globs are encouraged for readability. In the case of multiple files, genetic breaks are inferred from the beginnings and ends of files. Additionally, we tell :code:`smcsmc` the number of samples we expect with :code:`nsam`, so it knows to check the :code:`seg` file for formatting issues.

Basic Arguements
++++++++++++++++

In addition to the segment file, we need a few more pieces of information to kick off the particle filter.

**Particle Count:** We must specify a particle count (:code:`Np`) or use the default value of 1000. The number of particles refers to the number of individual ARGs which are simulated by :code:`SCRM`, higher numbers means that the model is more likely to converge on a reasonable answer, but increasing the particle count is computationally expensive. We recommend particle counts between 25 and 50 thousand for analysis. However, for exploratory work, smaller particle counts (a good starting point is 10 thousand) can be used to more effectively use computational resources.  In this guide, intended for use on a personal computer, we will use 10 particles.

**Number of Iterations:** We also explictly define the number of expectation maximization iterations that we want to use (:code:`EM`). A good place to start is 15 iterations. :code:`smcsmc` looks for output before starting, so should you be unsatisfied with the convergence after 15 iterations, simply run the exact same program with a higher number of iterations. :code:`smcsmc` will find the output for the iterations already run and simply continue from where it left off. 

**Demographic Parameters:** Several arguements are good to specify, especially when analysing real data as they help with convergence. Here we will use a mutation rate (:code:`mu`) of :code:`1.25e-8` and a recombination rate (:code:`rho`) of :code:`3e-9`. We give an effective population size (:code:`N0`) of 14312, and infer trees back to :code:`4*N0` generations with :code:`tmax`. We set this to 3.5, which is 1.2 Mya with an :code:`N0` of 10000. :code:`smcsmc` infers demographic parameters as discrete over intervals. We specify 31 intervals evenly spaced on the log scale with :code:`P 133 133032 31*1` giving the end times of the first and last epochs in generations, and the pattern for their generation. 

We also give the path to the output folder with :code:`o`.


Running :code:`smcsmc`
---------------------


.. tip:: 
        It is always a good idea to start a new :code:`conda` environment for a new analysis to ensure that there are no dependecy conflicts, though you can skip this step if you certain that your environment is correctly configured.

        .. code-block:: bash

                conda create --name smcsmc_tutorial
                conda activate smcsmc_tutorial

Once :code:`smcsmc` is installed, we can format the arguements detailed above into a dictionary. 

.. code-block:: python

   args = {
        'seg':    'test.seg',
        'nsam':   '4',
        'Np':     '10',
        'EM':     '1',
        'mu':     '1.25e-8',
        'rho':    'rho',
        'N0':     '10000',
        'tmax':   '3.5'
        'P':      '133 133032 31*1'
        'o':      'smcsmc_output'
   }

We directly use this dictionary with the :code:`run_smcsmc` command, which takes as its only arguement a dictionary of arguements.


.. code-block:: python

        import smcsmc
        smcsmc.run_smcsmc(args)

If your installation has been successful, then this will begin the process of parsing the input, merging any given :code:`seg` files, starting the particle filter, and iterating through the :code:`EM` steps requested.

Output
------

If your :code:`smcsmc` has run correctly, the resulting output directory will look something like this, with a seperate folder for each EM iteration, and a seperate file for each chunk, if this option has been used.


.. code-block:: bash
        
        output/
                emiterN/
                        chunkN.out
                        chunkN.stdout
                        chunk.stderr
                        chunkfinal.out
                merged.seg
                merged.map
                result.log
                result.out

If you are following this tutorial and are only using a single input :code:`seg` file, you will not see :code:`merged.seg` or :code:`merged.map` as there was no need to generate them. The output for each chunk is given, along with stdout and sterr, and results are aggregated over all chunks each epoch into :code:`chunkfinal.out`. The final epoch will be post processed into :code:`result.out`. Output and debugging information along with useful information to help interpret the results of your model are given in :code:`result.log`.

An example of a `results.out` file is given here:

.. code-block: bash

        EMstep EpochIndex EpochBegin   EpochEnd  EventType    FromPop      ToPop     Opportunity           Count            Rate              NE             ESS
             0          0   0.000000        Inf       Coal          0         -1     1.43803e+07         720.833     5.01264e-05     9974.779387               6
             0         -1   0.000000        Inf     Recomb         -1         -1      7.3339e+11         719.833     9.81516e-10              -1               6
             1          0   0.000000        Inf       Coal          0         -1     1.59149e+07         793.167      4.9838e-05    10032.496392               6
             1         -1   0.000000        Inf     Recomb         -1         -1     8.20343e+11         792.167     9.65653e-10              -1               6
             2          0   0.000000        Inf       Coal          0         -1     1.53825e+07             769      4.9992e-05    10001.597264               6
             2         -1   0.000000        Inf     Recomb         -1         -1     7.87734e+11             768     9.74948e-10              -1               6
             3          0   0.000000        Inf       Coal          0         -1     1.59971e+07         800.333       5.003e-05     9994.012474               6
             3         -1   0.000000        Inf     Recomb         -1         -1     8.21027e+11         799.333     9.73577e-10              -1               6
             4          0   0.000000        Inf       Coal          0         -1     1.57724e+07         771.833     4.89357e-05    10217.492772               6
             4         -1   0.000000        Inf     Recomb         -1         -1     7.97319e+11         770.833     9.66782e-10              -1               6
             5          0   0.000000        Inf       Coal          0         -1     1.61938e+07         796.333     4.91752e-05    10167.727706               6
             5         -1   0.000000        Inf     Recomb         -1         -1     8.21524e+11         795.333      9.6812e-10              -1               6
             6          0   0.000000        Inf       Coal          0         -1     1.49512e+07           752.5     5.03303e-05     9934.372706               6
             6         -1   0.000000        Inf     Recomb         -1         -1     7.64212e+11           751.5     9.83366e-10              -1               6
             7          0   0.000000        Inf       Coal          0         -1     1.68871e+07           846.5     5.01271e-05     9974.637421               6
             7         -1   0.000000        Inf     Recomb         -1         -1     8.43795e+11           845.5     1.00202e-09              -1               6
             8          0   0.000000        Inf       Coal          0         -1     1.63813e+07         816.333     4.98333e-05    10033.452317               6
             8         -1   0.000000        Inf     Recomb         -1         -1      8.2094e+11         815.333     9.93171e-10              -1               6
             9          0   0.000000        Inf       Coal          0         -1     1.72288e+07         862.333     5.00518e-05     9989.652610               6
             9         -1   0.000000        Inf     Recomb         -1         -1     8.82067e+11         861.333     9.76494e-10              -1               6
            10          0   0.000000        Inf       Coal          0         -1     1.48676e+07             751     5.05127e-05     9898.501864               6
            10         -1   0.000000        Inf     Recomb         -1         -1     7.67899e+11             750     9.76691e-10              -1               6



