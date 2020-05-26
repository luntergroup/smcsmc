Simons Global Diversity Panel
=============================

In this tutorial we will use :code:`smcsmc` to analyse the population history of a population from the `Simons Global Diversity Panel <https://www.simonsfoundation.org/simons-genome-diversity-project/>`_. 

We very arbitrarily chose to analyse the population history of the Mbuti, a South-Central African group with an ancient divergence from the rest of the continent. 

Downloading and Converting Data
--------------------------------

We will only download one chromosome of data for the sake of this tutorial.

First, make a new folder to hold your analysis of the Mbuti:

.. code-block:: sh

        mkdir smcsmc-mbuti
        cd smcsmc-mbuti

We will download a small chromosome from the phased release of the SGDP provided by the Reich Lab and save it as :code:`sgdp_phased_chr21.vcf.gz`:

.. code-block:: sh

        wget https://sharehost.hms.harvard.edu/genetics/reich_lab/sgdp/phased_data/PS2_multisample_public/cteam_extended.v4.PS2_phase.public.chr21.vcf.gz -O sgdp_phased_chr21.vcf.gz 

We then need to convert this sample from :code:`vcf` file format to the segments which :code:`smcsmc` requires. We can use :func:`smcsmc.vcf_to_seg`. 

.. code-block:: python

        import smcsmc

        smcsmc.vcf_to_seg(
                [('sgdp_phased_chr21.vcf.gz', 'S_Mbuti-1')],
                output = 'mbuti_chr21.seg.gz',
                chroms = [21])

.. warning::

        This can take a few minutes.

Now we have a list of all the segments in the Mbuti individual. 

Running the Model
------------------

We are going to follow a very similar procedure as the :ref:`getting_started` guide:

.. code-block:: python

        import os

        smcsmcpath = os.path.expandvars('${CONDA_PREFIX}/bin/smcsmc')

        args = {
                'EM':                           '15',
                'Np':                           '10000',

                # Submission Parameters
                'chunks':                       '100',
                'no_infer_recomb':              '',

                'smcsmcpath':                   smcsmcpath,

                # Other inference parameters
                'mu':                           '1.25e-8',
                'N0':                           '14312',
                'rho':                          '3e-9',
                'calibrate_lag':                '1.0',
                'tmax':                         '3.5',
                'alpha':                        '0',
                'apf':                          '2',
                'P':                            '133 133016 31*1',
                'VB':                           '',
                'nsam':                         '2',

                # Files
                'o':                            'run',
                'seg':                          'mbuti_chr21.seg.gz'
        }

If you are able, add in the :code:`c` flag to run on a cluster compute system.

Then simply run the model.

.. code-block:: python

        smcsmc.run_smcsmc(args)

.. todo:: 

        Run this model and include a plot of the results.
