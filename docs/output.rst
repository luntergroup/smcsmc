.. _sec-output:

==========================
Making sense of the output
==========================


************
Output files
************

``smcsmc`` outputs text files with user-specified prefix with flag **-o**.

*prefix*.log
    Log file records both ``smcsmc`` and ``scrm`` versions (git commits),
    input file paths, parameters used and population estimation at the final EM iteration.

*prefix*.out
    ``smcsmc`` records the inference results in the `*.out` file. Each row
    contains inference result of the population parameters including
    "population sizes", "migration rates" at different times and overall
    recombination rate.

* ``EMstep`` EM iteration.
* ``EpochIndex`` Epoch (time interval) index. ``-1`` for recombination rate inference.
* ``EpochBegin`` Epoch time begin (recent).
* ``EpochEnd`` Epoch time end (acient).
* ``EventType`` `Type of event`, ``Coal`` refers to coalescent event; ``Recomb`` refers to recombination event; ``Migr`` refers to migration event.
* ``FromPop`` Population index for coalescent rate inference, and migration (from) rate inference, 0-indexed. ``-1`` for recombination rate inference.
* ``ToPop`` Population index for migration (to) rate inference. 0-indexed. ``-1`` for recombination and coalescent rate inference.
* ``Opportunity`` Time (time times sequence for recombination events) span :math:`o` that `type of event` can happen.
* ``Count`` Number of `type of event` :math:`c` happened during the simulation.
* ``Rate`` Coalescent, recombination and migration rate, equal to :math:`c/o`.
* ``NE`` Population size. Calculated by :math:`o/2c`. ``-1`` for recombination and migration rate inference.
* ``ESS`` Effective sample size, .

For example:
::

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


******************************
Example of output interpretion
******************************

.. todo::
    Examples to come
