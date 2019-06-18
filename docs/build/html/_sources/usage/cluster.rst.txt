.. _cluster:

Clustering Computing 
====================

:code:`smcsmc` uses a lot of compute time. To ammeliorate this, it is *highly* recommend that users should make use of compute clusters, whether academic or otherwise. In this growing guide, we show how to leverage various compute cluster environments to best use :code:`smcsmc`. Please help us to fill in this section if you use a different cluster computing system.


SGE / :code:`qsub`
-----------------

The code ships with native support for SGE clusters and usage is relatively straightforward. Procedure for PBS clusters should be almost identical, and only minor tweaking of :code:`model.py` will give identical functionality.

.. tip::

        It's a good idea to run :code:`smcsmc` in a :code:`tmux` or :code:`screen` session on the head node of your compute cluster. This will ensure that if your connection is interupted the program will not stopped. 

1. If you wish to use a queue other than your default, place a :code:`qsub.conf` file **in the same directory that you want to conduct your analysis.** 

2. Specify a number of chunks in your input arguements. A good number to start with for typical human sequences is 100. 

3. Tell :code:`smcsmc` that you wish to spawn cluster jobs by specifying the :code:`c` option. Note that this overrides the use of multithreading. If you wish to use specific parameters for your cluster session, the :code:`C` flag can be used to provide specific directives and will take priority over **all** instructions given in :code:`qsub.conf`.

That's it! 

As an example, the following would be **added** to any other input arguements that you wish to use:

.. code-block:: python

   cluster_args = {
      'c':      '',
      'chunks': '100',
      'C':      '-P your.prc -q short.prj'
   }

For more information about how input arguements are handled by :code:`smcsmc` see 
