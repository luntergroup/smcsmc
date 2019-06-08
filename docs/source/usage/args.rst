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

:code:`smcsmc` will also understand that an argument passed with :code:`None` as a value will be removed entirely. This is a convenience for reusing input.  

Arguments
-------------------------------------------
