.. _sec-output:

==========================
Making sense of the output
==========================


************
Output files
************

``smcsmc`` outputs text files with user-specified prefix with flag **-o**.

*prefix*.log
    Log file records both ``smcsmc`` and ``scrm`` versions (git commits), input file paths, parameters used and population estimation at the final EM iteration.

*prefix* HIST
    Histrory, each iteration is seperated by =========, for example.

::

    =========
    RE	5.26809e-09
    ME	0	0
    ME	0.009672	0
    ME	0.025538	0
    ME	0.045624	0
    ME	0.07105	0
    ME	0.103238	0
    ME	0.143983	0
    ME	0.195563	0
    ME	0.260858	0
    ME	0.343515	0
    ME	0.44815	0
    ME	0.580607	0
    ME	0.748285	0
    ME	0.960548	0
    ME	1.22925	0
    ME	1.5694	0
    ME	2	0
    NE	0	0.172105
    NE	0.009672	0.118122
    NE	0.025538	0.0528732
    NE	0.045624	0.267541
    NE	0.07105	0.915965
    NE	0.103238	0.795306
    NE	0.143983	0.915542
    NE	0.195563	0.595508
    NE	0.260858	0.498153
    NE	0.343515	0.437207
    NE	0.44815	0.585799
    NE	0.580607	0.557354
    NE	0.748285	0.457153
    NE	0.960548	0.433858
    NE	1.22925	1.10365
    NE	1.5694	1.92134
    NE	2	2.49941

*prefix* Ne
    The final EM iteratio, last iteration of **prefix** HIST.

*prefix* Count
    The total count of oppurtunities and event count, used for recombination rate inference.

******************************
Example of output interpretion
******************************
