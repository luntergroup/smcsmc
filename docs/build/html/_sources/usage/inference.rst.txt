Overview of Inference Procedure
================================

:code:`smcsmc` uses a particle and several novel methodological advancements to estimate epoch specific demographic parameters from sequence. Here we outline the general inference procedure, however for those interested, the complete procedure is given in the preprint. 

We observe polymorphisms across each of :math:`n` sequences and model a continuous-time stochastic process with piecewise constant realisations :math:`x:[1,L) \rightarrow \mathbb{T}` where :math:`\mathbb{T}` is the state space of the Markov process. In our case, :math:`\mathbb{T}` represents the set of genealogical trees over :math:`n` sequences. Our goal is to sample from the posterior distribution of genealogical trees. 

We extend the notion of a particle filter to Markov jump processes which are continuous in both time (sequence) and space (the uncountable set :math:`\mathbb{T}`). In this case, we are able to use Radon-Nikodym derivates as opposed to ratios of probability distributions. 

Sampling Importance Resampling
--------------------------------

Particle filters build up a sample from the posterior sequentially using the notion of Monte Carlo approximation. 

Suppose we want to find the expectation of some function :math:`f` and some hidden variable :math:`z`. We would like to draw samples from the posterior distribution :math:`p(z_n | X_n)` where :math:`X_n = (x_1, \dots, x_n)` are observed. Using Monte Carlo approximation:

.. math::
        
        \mathbb{E}_{\color{red}{p}} \left[ f \left( z_n \right) \right] &= \int f(z_n) p \left( z_n | X_n \right) dz_n \\
                                &\approx \sum^L_{l=1} w^{(l)}_n f \left( z^{(l)}_n \right)

where the approximation holds as :math:`n \rightarrow \infty`. Sampling weights :math:`\left\{ w^{(l)}_n \right\}` are defined as 

.. math::
        
        w^{(l)}_n = \frac{p \left(X_n | z^{(l)}_n \right)}{ \sum^L_{m=1} p \left( x_n | z^{(m)}_n \right)}


and the set :math:`\left\{ z^{(l)}_n \right\}` along with their weights represent the posterior :math:`p(z_n|x_n)`. However, in many cases the true distribution :math:`p` is difficult to sample from, and importance sampling is used. In IS, samples are taken from a so-called "proposal" distribution :math:`q` and weights are adjusted accordingly such that

.. math::

        \mathbb{E}_{\color{red}{q}} \left[ f (z_n) \right] &= \int f(z_n) \frac{p(z_n | X_n) }{ q(z_n | X_n)} q(z_n | X_n) dz_n \\
                &\approx \sum^L_{l=1} \tilde{w}^{(l)}_n f \left( z^{(l)}_n \right)

With new weights 

.. math::

        \tilde{w}^{(l)}_n = w^{(l)}_n \frac{p(z_n | X_n) }{ q(z_n | X_n)} 

And remembering that in our case, the above case would use Radon-Nikodym derivatives as opposed to ratios of probability distributions; for brevity we have illustrated the discrete case.


Using a very similar procecure, we can sample from :math:`p(z_{s+1} | z_{1:s})` and adjust weights to approximate :math:`p(z_{s+1} | x_{s+1})`. See the relevant section in the preprint for more information on this step.

However, because samples from :math:`q` generally are not close to the target distribution :math:`p`, the proportion of particles which are close to the posterior's mode diminishes exponentially each iteration. To address this, a resampling step draws samples from the proposal distribution, assinging each a weight of :math:`\frac{1}{N}`. However, in order to avoid unnessarily inflating the variation at the current time, resampling is only performed when the effective sample size (ESS) drops below a certain threshold.  

The "waypoints", or positions at which we consider resampling, require some attention. Too few waypoints will increase the variation of the approximation, while too many will impact computational efficiency without any gains in accuracy. We derive a criteria that avoids particle degeneracy and find that if waypoints are considered at each observation and additional waypoints are added such that waypoints are never more than :math:`\frac{1}{\sqrt{2 \sigma^2}` then the ESS will not drop more than :math:`\sqrt{\frac{1}{e}}` between waypoints. This avoids degeneracy.

In this way, we sequentially build up an approximation of the posterior.

A Lookahead Likelihood
-----------------------

Not every particle is equally important for approximating the posterior. Upweighting particles which will be relevant in the future at the expense of those only relevant to the current time will increase the overall accuracy of the algorithm. We extend Pitt and Shephard's discrete time Auxiliary Particle Filter to our continuous-time case. 

We evaluate the likelihood (though this term is used loosely) of the current state :math:`x^{(i)}_s` at some future point which is not too far beyond the current point, though how far ahead will depend on how well this lookahead distribution fits the true distribution. 

We modify the above sampling algorithm to incorporate a second set of "resampling weights" which incorporate this lookahead likelihood. When resampling is performed, particles are weighted by this second, informed, weight, to better adapt to future variation.

In practice, we require a tractable approximate likelihood of future data conditioned on a given genealogy. Several simplifications are neccessary, and we primarily incorporate singletons and doubletons which are informative about topology near the tips of the genealogy.  We primarily look at the distance :math:`s_i` to the nearest future singleton for each sequence, and the mutually consistent cherries with their supporting doubletons. We derive an approximation of the likelihood of the current genealogy given these data. 

Parameter Inference
--------------------

Parameters may be inferred either by stochastic expectation maximization (SEM) or via Variational Bayes (VB). SEM maximizes the expected log likelihood over a posterior distribution of latent variables (as generated above) and gives a maximum a posteriori (MAP) estimate of parameter values :math:`\theta`. VB can substantially improve on SEM by iteratively estimating the posterior distirbution of :math:`\theta` rather than taking a point estimate. In practice, this distribution is not reported in the current implementation, however we find that using VB is useful in avoiding fixed-point solutions in SEM.
