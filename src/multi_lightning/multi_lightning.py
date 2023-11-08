import sys
from pprint import pprint
import numpy as np
import astropy.units as u
import emcee

class MultiLightning:
    '''A class interface to fit "multi-region" spectral energy distributions: SEDs where multiple regions
    are resolved in some, but not all of the available bandpasses.
    '''

    def __init__(self,
                 lgh,
                 flux_obs,
                 flux_unc,
                 model_unc=None,
                 Nregions=None,
                 reg_names=None
                 ):
        '''Constructor.

        Parameters
        ----------
        lgh : Lightning object or list of Nregions Lightning objects
            If a single Lightning object, the same model will be applied to each region (albeit with different parameters)
            and the `Nregions` keyword must be set. More generally, this can by a list of Lightning objects, one per region.
            Things like the redshift, luminosity distance, and filter set are expected to agree between different objects.
        flux_obs : np.ndarray, (Nregions+1, Nfilters), float
            An array giving the observed fluxes in mJy: the first row should contain the unresolved "global" fluxes, and
            subsequent rows should contain the resolved fluxes for each region. Missing or NA entries should be set to NaN.
            As an example, suppose we have 2 regions with 4 resolved optical measurements and 2 unresolved IR measurements.
            Our flux array would be structured like
                              Opt1   Opt2   Opt3   Opt4  IR1  IR2
            unresolved     [[  NaN,   NaN,   NaN,   NaN, 0.2,   1],
            region1         [0.003, 0.002, 0.004, 0.002, NaN, NaN],
            region2         [  NaN, 0.001, 0.002, 0.001, NaN, NaN]]

            where we've supposed that region2 is undetected in band Opt1.
        flux_unc : np.ndarray, (Nregions+1, Nfilters), float
            An array giving the uncertainties on the fluxes in mJy. Here, missing and NA entries should be set to 0.0.
        model_unc : float
            Fractional model uncertainty to apply to all bands.
        Nregions : int
            Number of resolved regions. Only necessary if `lgh` is a single Lightning object.
        reg_names : list of str
            Names for individual regions, e.g. 'core', 'disk', 'clumpA', 'clumpB', etc. Defaults to 'reg1', 'reg2', etc.
            if not set.
        '''

        try:
            self.Nregions = len(lgh)
            self.lgh = lgh
            self.singlemodel = False
            self.Nfilters = lgh[0].Nfilters
            self.redshift = lgh[0].redshift
            self.DL = lgh[0].DL
        except TypeError as e:
            if (Nregions is None):
                print('If `lgh` is a single object, the `Nregions` keyword must be set!')
                raise(e)
            else:
                self.Nregions = int(Nregions)
                self.lgh = lgh
                self.singlemodel = True
                self.Nfilters = lgh.Nfilters
                self.redshift = lgh.redshift
                self.DL = lgh.DL

        if (reg_names is not None):
            assert (len(reg_names) == self.Nregions), 'Length of `reg_names` does not match `Nregions`.'
            self.reg_names = reg_names
        else:
            self.reg_names = ['reg%d' % (i+1) for i in np.arange(self.Nregions)]

        if self.singlemodel:
            self.Nparams = {self.reg_names[i]: lgh.Nparams for i in np.arange(self.Nregions)}
        else:
            self.Nparams = {self.reg_names[i]: lgh[i].Nparams for i in np.arange(self.Nregions)}

        # Convert fnu to lnu for likelihood computation
        fnu2lnu = 4 * np.pi * (self.DL * u.Mpc) ** 2
        self.lnu_obs = (fnu2lnu * (flux_obs * u.mJy)).to(u.Lsun / u.Hz).value
        self.lnu_unc = (fnu2lnu * (flux_unc * u.mJy)).to(u.Lsun / u.Hz).value

        self.unres_mask = ~np.isnan(self.lnu_obs[0,:])

        assert (self.lnu_obs.shape[1] == self.Nfilters), 'length of filter axis in flux array is not the same as `Nfilters`'
        assert (self.lnu_obs.shape == self.lnu_unc.shape), 'Shapes of flux array and uncertainty array do not match'

        if (model_unc is not None):
            self.model_unc = model_unc
        else:
            self.model_unc = 0.0

    def print_params(self, verbose=False):
        '''List parameters (or pretty-print, when verbose=True)'''

        if self.singlemodel:
            print('%d regions, model is the same for each:' % (self.Nregions))
            self.lgh.print_params(verbose=verbose)
        else:
            for i, l in enumerate(self.lgh):
                print('Region %d:' % (i+1))
                l.print_params(verbose=verbose)

    def get_model_log_prior(self, params, priors):
        '''Calculate prior probability of parameters.

        Parameters
        ----------
        params : dict
            A dictionary keyed on `reg_names`, where each entry is a numpy array with dimensions (Nmodels, Nparamsx),
            giving the parameters for each region.
            Note that in the most general case `Nparamsx` can be different for each region, while the first `Nmodels`
            axis is the vectorization axis, and should be the same for each region.
        priors : dict
            A dictionary keyed on `reg_names`, where each entry is a list of the `Nparamsx`-many prior functions for
            each region. Prior functions should be specified by lightning.priors objects.

        It's pretty self evident how you might have parameters *loosely* connected between regions (e.g. region2's tauV
        is only allowed to be so far away from region1's) but harder to figure out how to *exactly* connect them, since
        this is not so much a prior as a reduction in dimensionality.
        '''

        Nmodels = params[self.reg_names[0]].shape[0]

        prior_prob = 1 + np.zeros(Nmodels)

        for reg in self.reg_names:

            reg_params = params[reg]
            reg_priors = priors[reg]

            for j in np.arange(len(reg_priors)):
                if (reg_priors[j] is not None):
                    # The 'connection'-type priors require both parameters and work slightly different
                    # on the backend, so they're separated here.
                    if ('connection' not in reg_priors[j].model_name):
                        prior_prob *= reg_priors[j](reg_params[:,j])
                    else:
                        prior_prob *= reg_priors[j].evaluate(reg_params[:,j], params[reg_priors[j].target_name][:,reg_priors[j].param_idx])

        lnprior_prob = np.zeros_like(prior_prob)
        lnprior_prob[prior_prob == 0] = -1*np.inf
        lnprior_prob[prior_prob != 0] = np.log(prior_prob[prior_prob != 0])

        return lnprior_prob

    def get_model_log_like(self, params):
        '''Calculate model log-likelihood.

        Parameters
        ----------
        params : dict
            A dictionary keyed on `reg_names`, where each entry is a numpy array with dimensions (Nmodels, Nparamsx),
            giving the parameters for each region.
            Note that in the most general case `Nparamsx` can be different for each region, while the first `Nmodels`
            axis is the vectorization axis, and should be the same for each region.

        '''

        Nmodels = params[self.reg_names[0]].shape[0]
        Lmod_arr = np.zeros((self.Nregions + 1, Nmodels, self.Nfilters))

        Lmod_unres = np.zeros((Nmodels, self.Nfilters))
        chi2_total = np.zeros(Nmodels)

        for i,reg in enumerate(self.reg_names):
            if self.singlemodel:
                Lmod,_ = self.lgh.get_model_lnu(params[reg])
            else:
                Lmod,_ = self.lgh[i].get_model_lnu(params[reg])

            if Nmodels == 1: Lmod = Lmod.reshape(1,-1)

            lnu_unc_total = np.sqrt((self.lnu_unc[i+1,~self.unres_mask])[None,:]**2 + (self.model_unc * Lmod[:,~self.unres_mask])**2)
            chi2_reg = np.nansum((Lmod[:,~self.unres_mask] - (self.lnu_obs[i+1,~self.unres_mask])[None,:])**2 / lnu_unc_total**2, axis=1)
            chi2_total = chi2_total + chi2_reg

            # print('Lmod_%s:' % (reg))
            # print(Lmod)
            Lmod_unres = Lmod_unres + Lmod

        lnu_unc_total = np.sqrt((self.lnu_unc[0,self.unres_mask])[None,:]**2 + (self.model_unc * Lmod_unres[:, self.unres_mask])**2)

        chi2_unres = np.nansum((Lmod_unres[:,self.unres_mask] - (self.lnu_obs[0,self.unres_mask])[None,:])**2 / lnu_unc_total**2, axis=1)
        # print('Lmod_unres:')
        # print(Lmod_unres)
        # print('chi2_res = %.2f' % chi2_total[0])
        # print('chi2_unres = %.2f' % chi2_unres[0])
        chi2_total = chi2_total + chi2_unres

        return -0.5 * chi2_total


    def get_model_log_prob(self, params, priors, p_bound=np.inf):
        '''Calculate model log probability.

        Parameters
        ----------
        params : dict
            A dictionary keyed on `reg_names`, where each entry is a numpy array with dimensions (Nmodels, Nparamsx),
            giving the parameters for each region.
            Note that in the most general case `Nparamsx` can be different for each region, while the first `Nmodels`
            axis is the vectorization axis, and should be the same for each region.
        priors : dict
            A dictionary keyed on `reg_names`, where each entry is a list of the `Nparamsx`-many prior functions for
            each region. Prior functions should be specified by lightning.priors objects.
        p_bound : float
            The magnitude of the log-probability of zero (default: np.inf)
        '''

        Nmodels = params[self.reg_names[0]].shape[0]
        ob_mask = np.zeros(Nmodels, dtype='bool')

        for i,reg in enumerate(self.reg_names):
            if self.singlemodel:
                ob_mask = ob_mask | self.lgh._check_bounds(params[reg])
            else:
                ob_mask = ob_mask | self.lgh[i]._check_bounds(params[reg])

        ib_mask = ~ob_mask

        if np.count_nonzero(ib_mask) > 0:

            # construct dict for every region again, with the out-of-bounds models
            # removed
            params_upd = {reg: params[reg][ib_mask,:] for reg in self.reg_names}

            log_prior = np.zeros(Nmodels)
            log_prior[ob_mask] = -1 * p_bound
            log_prior[ib_mask] = self.get_model_log_prior(params_upd, priors)

            log_like = np.zeros(Nmodels)
            log_like[ib_mask] = self.get_model_log_like(params_upd)

            log_prob = log_prior + log_like

            # print('logP = %.2f' % (log_prob[0]))
            # print('logL = %.2f' % (log_like[0]))
            # print('logp = %.2f' % (log_prior[0]))

        else:
            log_prob = np.zeros(Nmodels) - p_bound


        return log_prob

    def _params_vec2dict(self, x, p0, priors):
        '''Take the squished parameter array that emcee works with and turn it back into a dict keyed on region names.

        Parameters
        ----------
        x : np.ndarray, (Nwalkers, Ndim)
            emcee's working array, with dimensions set by the number of walkers and the number of sampled
            dimensions in the model.
        p0 : dict
            A dictionary keyed on `reg_names`, where each entry is an array giving the initial parameters for
            each region. We need this here to figure out what values the contant parameters have.
        priors : dict
            A dictionary keyed on `reg_names`, where each entry is a list of the `Nparamsx`-many prior functions for
            each region. Prior functions should be specified by lightning.priors objects. Any prior == None is assumed
            to indicate a constant parameter. We need these here to know which parameters are constant and which are
            fixed to another region.

        Returns
        -------
        The parameters (including constant and fixed parameters) sorted into a dictionary, keyed on `reg_names`.

        '''

        const_dim = {reg: np.array([pr is None for pr in priors[reg]]) for reg in self.reg_names}
        fixed_dim = {reg: np.array([(pr is not None) and (pr.model_name == 'fixed-connection') for pr in priors[reg]]) for reg in self.reg_names}
        Nfixed = {reg: np.count_nonzero(fixed_dim[reg]) for reg in self.reg_names}
        const_vals = {reg: (p0[reg].flatten())[const_dim[reg]] for reg in self.reg_names}
        Nconst = {reg: np.count_nonzero(const_dim[reg]) for reg in self.reg_names}

        params = {}
        start = 0
        Nmodels = x.shape[0]
        for reg in self.reg_names:

            Nparams_reg = self.Nparams[reg]
            Nconst_reg = Nconst[reg]
            Nfixed_reg = Nfixed[reg]
            Nvar_reg = Nparams_reg - Nconst_reg - Nfixed_reg

            var_mask = (~const_dim[reg]) & (~fixed_dim[reg])

            arr = np.zeros((Nmodels, Nparams_reg))
            arr[:,var_mask] = x[:,start:start+Nvar_reg]
            arr[:,const_dim[reg]] = const_vals[reg]

            params[reg] = arr
            #print('%s:' % (reg), params[reg][0,:])

            start = start + Nvar_reg

        # Looping over the regions twice every likelihood call is obviously not ideal in the case
        # of larger numbers of regions, but we can only figure out what the fixed values
        # are supposed to be once we know what parameters correspond to which region.
        for reg in self.reg_names:

            Nfixed_reg = Nfixed[reg]

            if Nfixed_reg > 0:

                idcs = np.flatnonzero(fixed_dim[reg])
                prpr = np.array(priors[reg])[fixed_dim[reg]]

                for i,pr in zip(idcs, prpr):
                    params[reg][:,i] = params[pr.target_name][:,pr.param_idx]

        return params

    def fit(self, p0, priors, Nwalkers=64, Nsteps=30000, init_sigma=1e-3, progress=True, savefile=None):
        '''Fit the multi-region model with emcee.

        Parameters
        ----------
        priors : dict
            A dictionary keyed on `reg_names`, where each entry is a list of the `Nparamsx`-many prior functions for
            each region. Prior functions should be specified by lightning.priors objects. Any prior == None is assumed
            to indicate a constant parameter.
        Nwalkers : int
            Number of MCMC samplers for the emcee affine-invariant algorithm (default: 64)
        Nsteps : int
            Number of steps to run the MCMC (default: 30000)
        init_sigma : float
            Sigma for gaussian ball initialization (default: 1e-3)
        '''

        const_dim = {reg: np.array([pr is None for pr in priors[reg]]) for reg in self.reg_names}
        # Add distinction between constant and 'fixed' values -- which are attached to a model component in a different
        # region.
        fixed_dim = {reg: np.array([(pr is not None) and (pr.model_name == 'fixed-connection') for pr in priors[reg]]) for reg in self.reg_names}

        def _log_prob_func(x):
            '''
            Because of the way emcee works, our log prob function can only accept
            a single array `x` of parameter values that has shape (Nmodels, sum(Nparam - Nconst))
            where the first axis is the arbitrary vectorization axis and the second is all of the model
            parameters, minus constant parameters, flattened into one axis. We must reconstruct the dictionary.
            Somehow.
            '''

            params = self._params_vec2dict(x, p0, priors)

            return self.get_model_log_prob(params, priors)

        # Initialize
        rng = np.random.default_rng()
        x0_seed = []
        for reg in self.reg_names:
            var_mask = (~const_dim[reg]) & (~fixed_dim[reg])
            x0_seed = x0_seed + list((p0[reg].flatten())[var_mask])
        x0_seed = np.array(x0_seed)

        x0 = rng.normal(loc=0, scale=init_sigma, size=(Nwalkers, len(x0_seed))) + x0_seed[None,:]

        if savefile is not None:
            backend = emcee.backends.HDFBackend(savefile)
            backend.reset(Nwalkers, len(x0_seed))
        else:
            backend = None

        sampler = emcee.EnsembleSampler(Nwalkers, len(x0_seed), _log_prob_func, vectorize=True, backend=backend)
        state = sampler.run_mcmc(x0, Nsteps, progress=progress)

        return sampler
