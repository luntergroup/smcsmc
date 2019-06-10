Simulated Data with msprime
===========================

Another popular framework for coalescent simulation is :code:`msprime`, a reimplimentation of Hudson's :code:`ms`. In this tutorial we show how to simulate data with :code:`msprime` under any model and reinfer demographic parmeters with :code:`smcsmc`.

Data Generation
---------------

This is not a guide for *using* :code:`msprime`, and if you are unfamiliar with syntax, the `documentation <https://msprime.readthedocs.io/en/stable>`_ is very helpful. Here we are essentially following the analysis implemented in the `PopSim analysis <https://github.com/popgensims/analysis>`_ comparing multiple tools for demographic inference.  

.. code-block:: python

        import msprime
        from stdpopsim import homo_sapiens

        # Here we use a three-population Out of Africa (OoA) model of human history.
        #
        # It is stored in the stdpopsim homo_sapiens class. 
        #
        # Your model could be anything.
        model = getattr(stdpopsim.homo_sapiens, "GutenkunstThreePopOutOfAfrica")

        # We are interested in four individuals from the second (European) population.
        #
        # We take these samples in the present day.
        samples = [msprime.Sample(population = 1, time = 0)] * 4

        # Perform the simulation
        # We pass all of the demographic events from the model as a dictionary to the simulation function.
        ts = msprime.simulate(
                samples = samples,
                mutation_rate = 1.25e-8,
                **model.asdict())

        # And spit out the tree file.
        ts.dump('msprime.tree')


We include functions to convert from tree sequence dumps to the :code:`seg` files necessary for :code:`smcsmc` analysis.

.. code-block:: python

        import smcsmc

        smcsmc.ts_to_seg('msprime.tree')

Which will write :code:`msprime.seg` in the current working directory. We now, fairly straightforwardly, analysing this seg data using the same procedure outlined in :ref:`getting_started`.

.. code-block:: python

        import smcsmc

        args = {
                'EM': 				'15',
                'Np': 				'10000',
                
                # Submission Parameters
                'chunks': 			'100',
                'c':				'',
                'no_infer_recomb': 	        '',

                # Other inference parameters
                'mu': 				'1.25e-9',
                'N0':				'14312',
                'rho':				'3e-9',
                'calibrate_lag':	        '1.0',
                'tmax':				'3.5',
                'alpha': 			'0',
                'apf': 				'2',
                'P': 				'133 133016 31*1',
                'VB':				'',
                'nsam':				'4',

                # Files
                'o':                            'run',
                'seg':                          'msprime.seg'
        }
