import numpy as np
from lightning.priors.base import AnalyticPrior
from scipy.special import erfinv

class NormalConnection(AnalyticPrior):
    r'''Treats the distance :math:`\delta x = x_1 - x_2` as a Normal random variable with :math:`\delta x \sim N(\mu, \sigma)`

    Parameters
    ----------
    params : array-like, (2,)
        Mean and sigma of the normal distribution.
    target_name : str
        Name of region to target
    param_idx : int
        Index of parameter in region model.

    '''

    type = 'analytic'
    model_name = 'normal-connection'
    Nparams = 2
    param_names = ['mu', 'sigma']
    param_descr = ['Mean', 'Sigma']
    param_bounds = np.array([[-np.inf, np.inf],
                             [0, np.inf]])

    def __init__(self, params, target_name, param_idx):

        self.params = params # i.e., mu and sigma
        # These two parameters are used to look up the target values
        self.target_name = target_name # which region is the parameter connected to?
        self.param_idx = param_idx # what is the index of the parameter in that region (sorry)

    def evaluate(self, x1, x2):
        r'''Calculate prior probability.

        Returns an array with the same shape as ``x1`` that is
        equal to

        .. math::

            p = \frac{1}{\sigma \sqrt{2 \pi}} \exp[- (\delta x - \mu)^2 / \sigma^2)]

        where

        .. math::

            \delta x = x_1 - x_2

        '''

        assert (x1.shape == x2.shape)

        mu = self.params[0]
        sigma = self.params[1]

        dx = x1 - x2

        p = 1 / (sigma * np.sqrt(2 * np.pi)) * np.exp(-1 * ((dx - mu) / sigma)**2)

        return p

    def quantile(self, q):
        r'''Calculate quantile function:

        Return an array with the same shape as ``q`` that is
        equal to

        .. math::

            \delta x = \mu + \sigma \sqrt{2} {\rm erfinv}(2 q - 1)

        where

        .. math::

            \delta x = x_1 - x_2

        '''

        mu = self.params[0]
        sigma = self.params[1]

        dx = mu + sigma * np.sqrt(2) * erfinv(2 * q - 1)

        return dx

class UniformConnection(AnalyticPrior):
    r'''Treats the absolute distance :math:`dx = |x_1 - x_2|` as a Uniform random variable with `dx ~ U(0, D_{max})`, i.e., x2 can be at most Dmax from x1.

    Parameters
    ----------
    params : array-like, (2,)
        Lower and upper bound of the uniform distribution.
    target_name : str
        Name of region to target
    param_idx : int
        Index of parameter in region model.

    '''

    type = 'analytic'
    model_name = 'normal-connection'
    Nparams = 1
    param_names = ['Dmax']
    param_descr = ['Maximum distance']
    param_bounds = np.array([0, np.inf])

    def __init__(self, param, target_name, param_idx):

        self.params = param # i.e., mu and sigma
        # These two parameters are used to look up the target values
        self.target_name = target_name # which region is the parameter connected to?
        self.param_idx = param_idx # what is the index of the parameter in that region (sorry)

    def evaluate(self, x1, x2):
        r'''Calculate prior probability.

        Return an array with the same shape as x1 that's equal to :math:`1 / D_{max}`
        wherever :math:`dx` is in :math:`[a,b)` and 0 elsewhere.
        '''

        Dmax = self.params

        assert (x1.shape == x2.shape)

        dx = np.abs(x1 - x2)

        p = np.zeros_like(x1)
        p[dx < Dmax] = 1 / (Dmax)

        return p

    def quantile(self, q):
        '''Calculate quantile function.

        Return an array with the same shape as q that's equal to :math:`q D_{max}`.
        '''

        Dmax = self.param

        dx = q * Dmax

        return dx

class FixedConnection(AnalyticPrior):
    '''Dummy prior. We use this to hold a parameter value fixed to a parameter from another region.

    Note that this prior doesn't fully implement ``lightning.priors.AnalyticPrior`` and can't
    be used to sample

    Parameters
    ----------
    target_name : str
        Name of region to target
    param_idx : int
        Index of parameter in region model.
    '''

    type = 'analytic'
    model_name = 'fixed-connection'
    Nparams = 0
    param_names = None
    param_descr = None
    param_bounds = None

    def __init__(self, target_name, param_idx):

        # These two parameters are used to look up the target values
        self.target_name = target_name # which region is the parameter connected to?
        self.param_idx = param_idx # what is the index of the parameter in that region (sorry)

    def evaluate(self, x1, x2):
        '''Returns x1 == x2'''

        return np.array(x1 == x2, dtype='int')

    def quantile(self, q):
        '''Not defined'''
        return None
