#!/usr/bin/env python

'''
test_postproc_NGC0628.py

Construct diagnostic plots for the NGC0628 example.
'''

from astropy.table import Table
from lightning import Lightning
from lightning.plots import step_curve, sfh_plot
import numpy as np
import h5py
import corner
import matplotlib.pyplot as plt
plt.style.use('ebm-dejavu')

def main():

    f = h5py.File('NGC0628_res.h5', 'r')

    lgh = Lightning.from_json('NGC0628_config.json')

    param_labels = []
    for mod in lgh.model_components:
        if mod is not None:
            param_labels = param_labels + mod.param_names_fncy

    param_labels_log = ['log' + l if 'psi' in l else l for l in param_labels]
    param_labels_log = np.array(param_labels_log)
    param_labels = np.array(param_labels)

    samples_center = f['samples/center'][:,:]
    samples_disk = f['samples/disk'][:,:]
    const_dim = np.var(samples_center, axis=0) < 1e-10

    sfh_center = samples_center[:,:5]
    sfh_disk = samples_disk[:,:5]

    mstar_center = np.sum(lgh.stars.mstar[None,:] * sfh_center, axis=1)
    mstar_disk = np.sum(lgh.stars.mstar[None,:] * sfh_disk, axis=1)

    sfh_center_lo, sfh_center_med, sfh_center_hi = np.quantile(sfh_center, q=(0.16, 0.50, 0.84), axis=0)
    sfh_disk_lo, sfh_disk_med, sfh_disk_hi = np.quantile(sfh_disk, q=(0.16, 0.50, 0.84), axis=0)
    mstar_center_lo, mstar_center_med, mstar_center_hi = np.quantile(mstar_center, q=(0.16, 0.50, 0.84), axis=0)
    mstar_disk_lo, mstar_disk_med, mstar_disk_hi = np.quantile(mstar_disk, q=(0.16, 0.50, 0.84), axis=0)

    # Plot log SFH so that we can fit both
    # regions on the same corner plot
    tmp = f['samples/center'][:,~const_dim]
    tmp[:,:5] = np.log10(tmp[:,:5])
    fig1 = corner.corner(tmp,
                        labels=param_labels_log[~const_dim],
                        quantiles=[0.16, 0.50, 0.84],
                        #show_titles=True,
                        smooth=1,
                        color='darkorange')
    tmp = f['samples/disk'][:,~const_dim]
    tmp[:,:5] = np.log10(tmp[:,:5])
    fig1 = corner.corner(tmp,
                        labels=param_labels_log[~const_dim],
                        quantiles=[0.16, 0.50, 0.84],
                        #show_titles=True,
                        smooth=1,
                        fig=fig1,
                        color='dodgerblue')

    fig1.savefig('NGC0628_corner.pdf')

    # Chain plot
    fig3, axs = plt.subplots(9,2, figsize=(9,9))

    t = 1 + np.arange(samples_center.shape[0])
    for i, label in enumerate(param_labels[~const_dim]):
        axs[i,0].plot(t, samples_center[:,~const_dim][:,i], color='darkorange')
        axs[i,1].plot(t, samples_disk[:,~const_dim][:,i], color='dodgerblue')

        axs[i,0].set_ylabel(label)
        if i != 8:
            axs[i,0].set_xticklabels([])
            axs[i,1].set_xticklabels([])

        axs[i,0].set_xlim(t[0],t[-1])
        axs[i,1].set_xlim(t[0],t[-1])

    axs[0,0].set_title('Center', color='darkorange')
    axs[0,1].set_title('Disk', color='dodgerblue')
    axs[8,0].set_xlabel('Trial Number')
    axs[8,1].set_xlabel('Trial Number')

    fig3.savefig('NGC0628_chains.pdf')

    # Best-fit SED plot
    fig4 = plt.figure()
    ax41 = fig4.add_axes([0.1, 0.4, 0.8, 0.5])
    ax42 = fig4.add_axes([0.1, 0.1, 0.8, 0.3])

    bestfit = np.argmax(f['logprob_samples'])
    disk_lnu_best = lgh.get_model_components_lnu_hires(f['samples/disk'][bestfit,:])
    center_lnu_best = lgh.get_model_components_lnu_hires(f['samples/center'][bestfit,:])
    disk_lnu_best_total,_ = lgh.get_model_lnu_hires(f['samples/disk'][bestfit,:])
    center_lnu_best_total,_ = lgh.get_model_lnu_hires(f['samples/center'][bestfit,:])
    disk_lmod_best_total,_ = lgh.get_model_lnu(f['samples/disk'][bestfit,:])
    center_lmod_best_total,_ = lgh.get_model_lnu(f['samples/center'][bestfit,:])
    # total_total
    total_lmod_best_total = disk_lmod_best_total + center_lmod_best_total

    delchi_disk = (f['phot/disk/lnu'][:] - disk_lmod_best_total) / np.sqrt(f['phot/disk/lnu_unc'][:]**2 + (0.10 * disk_lmod_best_total)**2)
    delchi_center = (f['phot/center/lnu'][:] - center_lmod_best_total) / np.sqrt(f['phot/center/lnu_unc'][:]**2 + (0.10 * center_lmod_best_total)**2)
    delchi_unresolved = (f['phot/unresolved/lnu'][:] - total_lmod_best_total) / np.sqrt(f['phot/unresolved/lnu_unc'][:]**2 + (0.10 * center_lmod_best_total)**2)

    ax41.plot(lgh.wave_grid_obs,
              lgh.nu_grid_obs*disk_lnu_best_total,
              color='dodgerblue',
              label='Disk',
              alpha=0.5)
    ax41.plot(lgh.wave_grid_obs,
              lgh.nu_grid_obs*center_lnu_best_total,
              color='darkorange',
              label='Center',
              alpha=0.5)
    ax41.plot(lgh.wave_grid_obs,
              lgh.nu_grid_obs*(center_lnu_best_total + disk_lnu_best_total),
              color='black',
              label='Total',
              alpha=0.5)
    ax41.errorbar(lgh.wave_obs,
                  lgh.nu_obs * f['phot/disk/lnu'],
                  yerr=lgh.nu_obs * f['phot/disk/lnu_unc'],
                  marker='D',
                  linestyle='',
                  color='dodgerblue',
                  markerfacecolor='dodgerblue',
                  capsize=2)
    ax41.errorbar(lgh.wave_obs,
                  lgh.nu_obs * f['phot/center/lnu'],
                  yerr=lgh.nu_obs * f['phot/center/lnu_unc'],
                  marker='D',
                  linestyle='',
                  color='darkorange',
                  markerfacecolor='darkorange',
                  capsize=2)
    ax41.errorbar(lgh.wave_obs,
                  lgh.nu_obs * f['phot/unresolved/lnu'],
                  yerr=lgh.nu_obs * f['phot/unresolved/lnu_unc'],
                  marker='D',
                  linestyle='',
                  color='k',
                  markerfacecolor='k',
                  capsize=2)

    ax41.set_xscale('log')
    ax41.set_yscale('log')
    ax41.set_ylabel(r'$\nu L_{\nu}~[\rm L_{\odot}]$')
    ax41.legend(loc='lower left')

    ax42.axhline(-1, color='slategray', linestyle='--')
    ax42.axhline(0, color='slategray', linestyle='-')
    ax42.axhline(1, color='slategray', linestyle='--')

    ax42.errorbar(lgh.wave_obs,
                  delchi_disk,
                  yerr=np.ones_like(delchi_disk),
                  marker='D',
                  color='dodgerblue',
                  markerfacecolor='dodgerblue',
                  linestyle='',
                  capsize=2.0)
    ax42.errorbar(lgh.wave_obs,
                  delchi_center,
                  yerr=np.ones_like(delchi_center),
                  marker='D',
                  color='darkorange',
                  markerfacecolor='darkorange',
                  linestyle='',
                  capsize=2.0)
    ax42.errorbar(lgh.wave_obs,
                  delchi_unresolved,
                  yerr=np.ones_like(delchi_unresolved),
                  marker='D',
                  color='k',
                  markerfacecolor='k',
                  linestyle='',
                  capsize=2.0)
    ax42.set_xscale('log')
    ax42.set_xlim(ax41.get_xlim())
    ax42.set_xlabel(r'Observed-Frame Wavelength [$\rm \mu m$]')
    ax42.set_ylim(-5,5)
    ax42.set_ylabel(r'Residual [$\sigma$]')

    fig4.savefig('NGC0628_sed.pdf')

    # SFH plot -- Could also normalize the SFH somehow to better show the
    # differences in the shapes.
    fig5, ax5 = plt.subplots()
    fig5, ax5 = sfh_plot(lgh, samples_center,
                         shade_kwargs={'color':'darkorange', 'alpha':0.3, 'zorder':0},
                         line_kwargs={'color':'darkorange', 'zorder':1},
                         ax=ax5)
    fig5, ax5 = sfh_plot(lgh, samples_disk,
                         shade_kwargs={'color':'dodgerblue', 'alpha':0.3, 'zorder':0},
                         line_kwargs={'color':'dodgerblue', 'zorder':1},
                         ax=ax5)

    ax5.set_ylabel(r'$\psi (t)~[\rm M_{\odot}~yr^{-1}]$')
    ax5.set_xlabel(r'Stellar Age [yr]')
    ax5.set_xscale('log')
    ax5.set_yscale('log')
    ax5.set_xlim(1e6,13.6e9)
    ax5.set_ylim(1e-5, 10)

    fig5.savefig('NGC0628_sfh.pdf')



if __name__ == '__main__':
    main()
