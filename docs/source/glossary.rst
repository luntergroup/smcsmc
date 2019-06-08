Glossary
========

.. glossary:: 

   Demographic Parameters
      Effective population size and migration rates discretised over a given number of epochs. Within each of these epochs, the rates are assumed to be **constant**. 

   Rainbow Plot
        A plot used to assess the convergence of EM. The demographic parameters inferred from each epoch are plotted in a seperate colour to visually determine if significant improvement could be made by letting the model run more iterations. 

   Particle
      The basic unit of a particle filter, in our case a particle represents an ancestral recombination graph simulated by :code:`SCRM`\ . Particles are simulated, and updated along the sequence by genetic information that they encounter.  Particles are resampled once the effective sample size of the population reaches a given threshold and are weighted according to their approximate likelihood. Asymtotically, this form of resampling importance sampling approaches the true posterior. 

   Particle Filter
      Particle filters, or Sequential Monte Carlo (SMC) methods, are a class of algorithms which use sampling importance resampling (SIR) to generate samples from the posterior distribution of a latent variable. In our case, this is the set of trees along the sequence. The posterior is approximated using weighted random samples denoted as "particles" drawn from a known, tractable distribution. See `Tulsyan, Gopaluni, and Khare 2016 <https://www.sciencedirect.com/science/article/pii/S0098135416302769>`_ for a very readable review of the basic principles behind particle filters for inference.

   Sequentially Markovian Coalescent
      This may refer to one of a number of variations on McVean and Cardin 2005. The authors derive an approximation of the coalescent with recombination in which lineages with no overlapping ancestral material may not coalesce. This approximation leads to a Markovian state space of genealogies along the sequence, and allows for tractable inference. 

   Variational Bayes
      A method used to infer the posterior distribution of demographic parameters from a set of particles. Variation Bayes is quite similar to stochastic Expectation Maximation, except that it generates a full posterior distribution of parameters rather than a single most probable value.
