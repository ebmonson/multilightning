#!/usr/bin/env python

import h5py
import sys
from pprint import pprint
import numpy as np
from astropy.table import Table

from lightning import Lightning
from lightning.priors import UniformPrior
from multi_lightning import MultiLightning
from multi_lightning.priors import FixedConnection, NormalConnection, UniformConnection

import matplotlib.pyplot as plt
plt.style.use('ebm-dejavu')

def main():

    phot = Table.read('NGC0628-photometry.fits')

    filter_labels = np.array([s.strip() for s in phot['FILTER_LABELS']])
    filter_labels[filter_labels == ['2MASS_K']] = '2MASS_Ks'

    Nfilters = len(filter_labels)

    #ascii.write(phot['FILTER_LABELS','WAVELENGTH'], format='fixed_width_two_line')

    fnu = phot['FNU_OBS'] * 1e3 # in mJy
    fnu_cent = phot['FNU_CENT_OBS'] * 1e3 # in mJy
    # The annulus photometry is the difference of the total and center regions.
    # In the future the input table may also contain fluxes for the annulus.
    fnu_out = fnu - fnu_cent

    cal_unc = np.array([0.15, 0.15, # GALEX
                        0.05, 0.05, 0.05, 0.05, 0.05, # SDSS
                        0.05, 0.05, 0.05, # 2MASS
                        0.05, 0.05, 0.05, 0.05, # IRAC
                        0.05, # MIPS 24 um
                        0.05, 0.05, 0.05, # PACS
                        0.15, 0.15, 0.15, # SPIRE
                        0.07, 0.07, 0.07, 0.07 # WISE
                        ])

    fnu_unc = cal_unc * fnu
    fnu_out_unc = cal_unc * fnu_out
    fnu_cent_unc = cal_unc * fnu_cent

    fnu_unc = cal_unc * fnu
    fnu_out_unc = cal_unc * fnu_out
    fnu_cent_unc = cal_unc * fnu_cent

    unresolved = (fnu != 0) & (fnu_cent == 0)
    resolved = ~unresolved

    fnu[fnu == 0] = np.nan
    fnu_out[fnu_cent == 0] = np.nan
    fnu_cent[fnu_cent == 0] = np.nan

    # The way the flux array is structured has changed to allow for more arbitrary
    # numbers and compositions of regions.
    fnu_arr = np.zeros((3, Nfilters))
    fnu_arr[1,:] = fnu_cent
    fnu_arr[2,:] = fnu_out
    fnu_arr[0,unresolved] = fnu[unresolved]
    fnu_arr[0,resolved] = np.nan

    fnu_unc_arr = np.zeros((3, Nfilters))
    fnu_unc_arr[1,:] = fnu_cent_unc
    fnu_unc_arr[2,:] = fnu_out_unc
    fnu_unc_arr[0,unresolved] = fnu_unc[unresolved]
    fnu_unc_arr[0,resolved] = 0.0

    reg_names = ['center', 'disk']

    redshift = 0.0
    DL = 7.3 # Mpc - collected from Moustakis+2010, originally measured by Sharina+(1996) from supergiant stars.

    # We need only a single Lightning object since
    # both regions use the same models.
    lgh = Lightning(filter_labels,
                    lum_dist=DL,
                    atten_type='Calzetti',
                    dust_emission=True)

    lgh.save_json('NGC0628_config.json')

    fig, ax = plt.subplots()

    ax.errorbar(lgh.wave_obs,
                lgh.nu_obs * fnu_arr[0,:],
                yerr=lgh.nu_obs * fnu_unc_arr[0,:],
                marker='D',
                linestyle='',
                capsize=2.0,
                label='unresolved')
    ax.errorbar(lgh.wave_obs,
                lgh.nu_obs * fnu_arr[1,:],
                yerr=lgh.nu_obs * fnu_unc_arr[1,:],
                marker='D',
                linestyle='',
                capsize=2.0,
                color='darkorange',
                label='center')
    ax.errorbar(lgh.wave_obs,
                lgh.nu_obs * fnu_arr[2,:],
                yerr=lgh.nu_obs * fnu_unc_arr[2,:],
                marker='D',
                linestyle='',
                capsize=2.0,
                color='dodgerblue',
                label='disk')

    ax.set_xlabel(r'Observed-Frame Wavelength $[\rm \mu m]$')
    ax.set_ylabel(r'$\nu f_{\nu}~[\rm mJy~Hz]$')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend(loc='lower right')

    fig.savefig('NGC0628_observed_sed.pdf')

    mlgh = MultiLightning(lgh,
                          fnu_arr,
                          fnu_unc_arr,
                          model_unc=0.10,
                          Nregions=2,
                          reg_names=reg_names)

    # This can also be the 'seed' for our initialization for the MCMC fit.
    # Note that since supplying the initial points for the MCMC would be a real
    # pain with this parameter construction, it's now handled inside the `fit` method.
    ptest = {'center': np.array([0.1,0.1,0.1,0.1,0.1,0.1,2,1,3e5,0.01,0.02]).reshape(1,-1),
             'disk': np.array([1,1,1,1,1,0.1,2,1,3e5,0.01,0.02]).reshape(1,-1)
            }

    ptest2 = {'center': np.array([1.52255314e-04, 2.06336822e-04, 7.45473387e-05, 2.29343426e-02,
                                  2.82105599e-05, 1.49994396e-01, 2.00000000e+00, 2.08370416e+00,
                                  3.00000000e+05, 1.36305518e-02, 4.44731813e-02]).reshape(1,-1),
               'disk': np.array([1.74203653e-01, 1.27957064e-01, 3.17438598e+00, 7.57288134e-02,
                                 5.45456659e-01, 4.75409991e-01, 2.00000000e+00, 2.08370416e+00,
                                 3.00000000e+05, 1.36305518e-02, 4.44731813e-02]).reshape(1,-1)
               }

    # We can make the prior process a little less obnoxious with
    # list arithmetic.
    # The prior `FixedConnection` takes two arguments: a region name and a parameter index. Here it fixes the
    # qPAH of the center region to the qPAH of the disk.
    # priors_cen = 5 * [UniformPrior([0,1])] + \
    #                  [UniformPrior([0,3])] + \
    #                  [None, UniformPrior([0.1, 25]), None, UniformPrior([0,1]), FixedConnection('disk', -1)]
    # The prior `NormalConnection` takes three arguments: a [mu, sigma] iterable, a region name, and a parameter index.
    # Here it loosely fixes the qPAH of the center region to the qPAH of the disk.
    priors_cen = 5 * [UniformPrior([0,1])] + \
                     [UniformPrior([0,3])] + \
                     [None, UniformPrior([0.1, 25]), None, UniformPrior([0,1]), NormalConnection([0,0.01], 'disk', -1)]
    priors_disk = 5 * [UniformPrior([0,10])] + \
                      [UniformPrior([0,3])] + \
                      [None, UniformPrior([0.1, 25]), None, UniformPrior([0,1]), UniformPrior([0.0047, 0.0458])]

    const_dim_cen = np.array([pr is None for pr in priors_cen])
    const_dim_disk = np.array([pr is None for pr in priors_disk])

    fixed_dim_cen = np.array([(pr is not None) and (pr.model_name == 'fixed-connection') for pr in priors_cen])
    fixed_dim_disk = np.array([(pr is not None) and (pr.model_name == 'fixed-connection') for pr in priors_disk])

    priors = {'center': priors_cen,
              'disk': priors_disk
             }

    lnlike = mlgh.get_model_log_like(ptest)
    lnprior = mlgh.get_model_log_prior(ptest, priors)
    lnprob = mlgh.get_model_log_prob(ptest, priors)

    lnlike2 = mlgh.get_model_log_like(ptest2)
    lnprior2 = mlgh.get_model_log_prior(ptest2, priors)
    lnprob2 = mlgh.get_model_log_prob(ptest2, priors)

    mcmc = mlgh.fit(ptest, priors,
                    progress=True)

    # Some of this should and will be put into MultiLightning, I'm just
    # sketching it out here.
    try:
        tau = mcmc.get_autocorr_time()
        print('tau = ', tau)
    except:
        print('Chains too short to properly estimate autocorrelation time.')
        print('Run a longer chain. Using tau = 500. Inspect products carefully...')
        tau = 500

    burnin = int(2.5 * np.max(tau))
    thin = int(0.5 * np.min(tau))

    samples = mcmc.get_chain(discard=burnin, flat=True, thin=thin)
    log_prob_samples = mcmc.get_log_prob(discard=burnin, flat=True, thin=thin)

    # Put the final chains (including constant parameters) in an hdf5 file.
    Nsamples_final = samples.shape[0]
    Nparams_cen = len(ptest['center'].flatten())
    Nparams_disk = len(ptest['disk'].flatten())

    params = mlgh._params_vec2dict(samples, ptest, priors)
    #samples_cen = np.zeros((Nsamples_final, Nparams_cen))
    #samples_disk = np.zeros((Nsamples_final, Nparams_disk))

    #var_dim_cen = (~const_dim_cen) & (~fixed_dim_cen)
    #var_dim_disk = (~const_dim_disk) & (~fixed_dim_disk)

    #samples_cen[:,var_dim_cen] = samples[:,0:8]
    #samples_cen[:,const_dim_cen] = (ptest['center'].flatten())[const_dim_cen][None,:]

    #samples_disk[:,var_dim_disk] = samples[:,8:]
    #samples_disk[:,const_dim_disk] = (ptest['disk'].flatten())[const_dim_disk][None,:]

    #samples_cen[:, fixed_dim_cen] = samples_disk[:,-1].reshape(-1,1)

    with h5py.File('NGC0628_res.h5', 'w') as f:
        f.create_dataset('phot/center/lnu', data=mlgh.lnu_obs[1,:])
        f.create_dataset('phot/center/lnu_unc', data=mlgh.lnu_unc[1,:])
        f.create_dataset('phot/disk/lnu', data=mlgh.lnu_obs[2,:])
        f.create_dataset('phot/disk/lnu_unc', data=mlgh.lnu_unc[2,:])
        f.create_dataset('phot/unresolved/lnu', data=mlgh.lnu_obs[0,:])
        f.create_dataset('phot/unresolved/lnu_unc', data=mlgh.lnu_unc[0,:])
        f.create_dataset('samples/center', data=params['center'])
        f.create_dataset('samples/disk', data=params['disk'])
        f.create_dataset('logprob_samples', data=log_prob_samples)



if __name__ == '__main__':
    main()
