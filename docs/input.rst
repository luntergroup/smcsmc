.. _sec-input:

=============
How it works?
=============

******************************
Program parameters and options
******************************

Mostly used
-----------

-seg [*file*]
    File path of the input data file in seg format. Each line of the segment
    file represent a segment of the sequence. The first column marks the
    position of segment start. The second column records the segment length,
    which is the difference between the start and the end of a segment.
    The third column records T(rue) if it's invariant. The forth column
    records T(rue) if it is a genetic break, F(alse) otherwise. The fifth
    column shows the chromosome number. The last column records the genotype,
    and use \. to denote missing, and \/\/ to denote unphased.

.. image:: _static/segPlot.png
   :width: 524px
   :alt: segment data

For example, a segment input file looks like this:

.. csv-table::
   :header: segment start, segment length, invariant, genetic break, chromosome, genotype
   :widths: 5, 5, 5, 5, 5, 5

    1,       521,     T,       F,       1,       01\.0
    522,     2721,    T,       F,       1,       0111
    3243,    1758,    T,       F,       1,       10\/\/
    5001,    1296,    T,       F,       1,       0000
    6297,    1,       T,       F,       1,       \.\.\.\.
    6298,    4669,    T,       F,       1,       0110
    10967,   880,     T,       T,       1,       0100
    1,       708,     T,       F,       2,       1010


-Np [int]
    Number of particles (default value 1000).

-ESS [float]
    Fractional ESS threshold for resampling, in range of :math:`(0, 1]`, and 1
    impiles to use random likelihoods (default value 0.6).


-tmax [float]
    Maximum time, in unit of :math:`4N_0` (default value :math:`3`, i.e. :math:`12N_0`).

-p [string]
    Pattern of time segments. Given a maximum
    TMRCA specified by **-tmax** flag, in the :math:`4N_0` scale, we divide the time
    interval :math:`[0, T_{max}]` to n atomic time intervals, for example, default
    pattern is "3*1+2*3+4". Therefore :math:`n = 3\times1+2\times3+4 = 13`. We
    then block adjacent atomic time intervals, and assume they
    have equal population sizes. For the default pattern, the first 3 blocks span
    1 interval, then each of the 4th and 5th intervals spans over three
    atomic intervals, the last interval spans over 4 atomic intervals [Li2011]_.


-EM [int]
    Number of EM iterations (default value 20).

-o [string]
    Prefix for output files

-online\
    Perform online EM

-xr [int]
    Epoch or epoch range to exclude from recombination EM (1-based, closed)

-xc [int]
    Epoch or epoch range (e.g. 1-10) to exclude from coalescent/migration EM

-log\
    Generate \*.log file.


You may also try
----------------

-heat\
    Generate \*TMRCA and \*WEIGHT for heatmap.

-v\
    Display timestamp and git versions.

.. _sec-eg:

***************************
Example of data exploration
***************************

.. [Li2011] Li, H. and R. Durbin (2011). Inference of human population history from individual whole-genome sequences. *Nature* 475, 493–496.
