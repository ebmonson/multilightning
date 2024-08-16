Background
==========

The original motivation for ``multilightning`` is fairly simple:

- The sub-arcsecond to arcsecond-scale resolution of UV-optical imaging of nearby (:math:`\lesssim 40` Mpc) galaxies contains
  detailed information about the spatial variation of star formation and dust obscuration.
- The comparatively worse :math:`\gtrsim 10-20` arcsecond-scale resolution of IR imaging loses much of that spatial information,
  but IR data provide a valuable handle on the overall normalization of the star formation history (SFH), and the inclusion of
  IR data in spectral energy distribution (SED) fits can allow us to use more complicated dust models.

We must therefore balance the need for spatial information with the need for an IR constraint on the SFH.

In the case of Lehmer et al. (2024), we required good constraints on the SFH over only the disks of some galaxies,
excluding the nuclear regions: the nuclear regions of star-forming galaxies suffer from crowding of X-ray sources,
such that individual sources can't be counted. However, these galaxies are not resolved in the available far-IR imaging.
Thus, instead of attempting to fit the SED and estimate the SFH using only the subset of fluxes that resolve the disk, we forward
model all the available data, using the IR data as a constraint on the SFH of both the disk and the nuclear region.

We have three sets of fluxes:

- Two sets of :math:`N_{\rm res}` fluxes, which resolve the two regions,
  :math:`\{F^{\rm disk}_{\nu,i}\}_{i=1}^{N_{\rm res}}` and :math:`\{F^{\rm center}_{\nu,i}\}_{i=1}^{N_{\rm res}}`
- One set of :math:`N_{\rm unres}` fluxes which cover the entire galaxy, :math:`\{F^{\rm total}_{\nu,j}\}_{j=1}^{N_{\rm unres}}`.

Each underlying region is modeled with a separate star formation history. We then have three :math:`\chi^2` terms:

.. math::

    \chi^2_{\rm disk} = \sum_{i=1}^{N_{\rm res}} \frac{(F_{\nu,i}^{\rm disk} - \hat F_{\nu,i}^{\rm disk})^2}{(\sigma F_{\nu,i}^{\rm disk})^2} \\
    \chi^2_{\rm center} = \sum_{i=1}^{N_{\rm res}} \frac{(F_{\nu,i}^{\rm center} - \hat F_{\nu,i}^{\rm center})^2}{(\sigma F_{\nu,i}^{\rm center})^2} \\
    \chi^2_{\rm total} = \sum_{j=1}^{N_{\rm unres}} \frac{(F_{\nu,j}^{\rm total} - \hat F_{\nu,j}^{\rm total})^2}{(\sigma F_{\nu,j}^{\rm total})^2}

where

.. math::

    \hat F_{\nu,j}^{\rm total} = \hat F_{\nu,j}^{\rm disk} + \hat F_{\nu,j}^{\rm center}.


We then have

.. math::

    \chi^2 = \chi^2_{\rm disk} + \chi^2_{\rm center} + \chi^2_{\rm total}

In `the included Jupyter notebook <examples/NGC628.ipynb>`_, we demonstrate this in practice for NGC 628.

.. This often results in building PSF-matched data cubes by convolving the available data to a common PSF, set by the worst resolution
.. of the data we use (c.f. Eufrasio et al. 2017, where multiwavelength data were convolved to a common :math:`25''` PSF to
.. estimate the SFH in a spatially-resolved way), smoothing out spatial variations on smaller scales than the PSF.
