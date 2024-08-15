import numpy as np
from lightning.priors.base import AnalyticPrior
from scipy.special import erfinv

class NormalConnection(AnalyticPrior):
    '''
    Treats the distance dx = x1 - x2 as a Normal random variable with dx ~ N(mu, sigma)
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
        '''
        Return an array with the same shape as ``x1`` that is
        equal to::

        p = 1 / [sigma * sqrt(2 * pi)] * exp[-1 * (x - mu)**2 / sigma**2)]

        '''

        assert (x1.shape == x2.shape)

        mu = self.params[0]
        sigma = self.params[1]

        dx = x1 - x2

        p = 1 / (sigma * np.sqrt(2 * np.pi)) * np.exp(-1 * ((dx - mu) / sigma)**2)

        return p

    def quantile(self, q):
        '''
        Return an array with the same shape as ``q`` that is
        equal to::

        dx = mu + sigma * sqrt(2) * erfinv(2 * q - 1)

        '''

        mu = self.params[0]
        sigma = self.params[1]

        dx = mu + sigma * np.sqrt(2) * erfinv(2 * q - 1)

        return dx

class UniformConnection(AnalyticPrior):
    '''
    Treats the absolute distance dx = |x1 - x2| as a Uniform random variable with dx ~ U(0, Dmax), i.e., x2 can be at
    most Dmax from x1.
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
        '''
        Return an array with the same shape as x1 that's equal to ``1 / Dmax``
        wherever x is in [a,b) and 0 elsewhere.
        '''

        Dmax = self.params

        assert (x1.shape == x2.shape)

        dx = np.abs(x1 - x2)

        p = np.zeros_like(x1)
        p[dx < Dmax] = 1 / (Dmax)

        return p

    def quantile(self, q):
        '''
        Return an array with the same shape as q that's equal to ``q * (Dmax)``.
        '''

        Dmax = self.param

        dx = q * Dmax

        return dx

class FixedConnection(AnalyticPrior):
    '''
    Dummy prior. We use this to hold a parameter value fixed to a parameter from another region.
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

        return np.array(x1 == x2, dtype='int')

    def quantile(self, q):

        return None
