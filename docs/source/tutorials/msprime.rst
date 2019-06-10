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

.. todo::

        Export the ts_to_seg function outside the utils submodule. Also make the number of haplotypes optional.

.. code-block:: python

        import smcsmc

        smcsmc.ts_to_seg('msprime.tree')

Which will write :code:`msprime.seg` in the current working directory.
