Glossary
========

.. glossary:: 

   Particle
      The basic unit of a particle filter, in our case a particle represents an ancestral recombination graph simulated by :code:`SCRM`\ . Particles are simulated, and updated along the sequence by genetic information that they encounter. Between iterations, particles are resampled according to their approximate likelihood, and asymtotically approach the true posterior in the limit of resampling steps.
