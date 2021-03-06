.. _args:

Input Arguments
================

This section describes the general structure of :code:`smcsmc` arguments and provides information on how they are interpreted. 

General Structure
-----------------

The user has a single entry point into :code:`smcsmc`, :code:`smcsmc.run_smcsmc`. This function takes a single main argument.  

Input arguments are always formatted into a `dict <https://developers.google.com/edu/python/dict-files>`_ and require both a key and a value. The key is always the name of the argument, and the value is generally the value of the arguement. In some cases they differ, and it is important to understand when this is the case. All values should be given as strings unless otherwise noted. This is simply a convenience to avoid complicated post processing.

.. code-block:: python
   :emphasize-lines: 3-7
   
   import smcsmc

   args = {
      'seg':    'test_seg.seg'
      'nsam':   '4'

   }

   smcsmc.run_smcsmc(args)

The arguments given entirely determine the :code:`smcsmc` inference. 

Processing of Arguments
------------------------

There are three main kinds of arguments to :code:`smcsmc`. 

* **Key value pairings**: These are the most common arguments, and require both a key and value. 

  + e.g. :code:`{'chunk': '100'}`

* **Boolean arguments**: To input a boolean, the name of the argument should be given and an empty string passed as its value. This will be picked up in processing and treated correctly.

  + e.g. :code:`{'no-infer-recomb': ''}`

* **Mulitple arguments of the same name**: Some arguements to :code:`smcsmc` are passed directly to :code:`SCRM` and act as psuedo-:code:`ms` code. In this case, pass arguements having the same name as a vector and :code:`smcsmc` will process them into the correct number of identically named arguements.

  + e.g. :code:`{'eM': ['0.0092 10', '0 5']}`

.. todo::

        This is not yet implemented.

:code:`smcsmc` will also understand that an argument passed with :code:`None` as a value will be removed entirely. This is a convenience for reusing input.  

Arguments
-------------------------------------------

:code:`[*]` arguments are required. :code:`[+]` arguments are optional, but one amungst the group is required.

The following are arguments that define properties of the population you are studying. They are fixed and do not change.

.. csv-table::
   :file: ../../args/pop_args.txt
   :widths: 20, 20, 60
   :header-rows: 1

The following arguements describe initial values which will be inferred and updated during runtime.

.. csv-table:: 
   :file: ../../args/inferred_args.txt
   :widths: 20,20, 60
   :header-rows: 1

The following arguements define inference related options.

.. csv-table::
   :file: ../../args/inference.txt
   :widths: 20,20,60
   :header-rows: 1

These arguments define the behaviour of the parameter updates via stochastic EM or Variational Bayes.

.. csv-table::
   :file: ../../args/em.txt
   :widths: 20,20,60
   :header-rows: 1

These are general options.

.. csv-table::
   :file: ../../args/general.txt
   :widths: 20,20,60
   :header-rows: 1
