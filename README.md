# MultiLightning
---

A package for jointly fitting multiple regions of a galaxy with the `Lightning` SED-fitting code.

For nearby galaxies, we often have a subsets of bands in which the galaxy is resolved and unresolved (e.g., the far-IR).
We can thus choose to either discard the spatial information of the resolved bands and fit the integrated light of the
galaxy, or fit without the unresolved data, potentially losing valuable constraints on the star formation rate.

Here, we define a joint likelihood for fitting multiple resolved regions, under the assumption that the unresolved
photometry is described by the sum of the models for each resolved region. By way of example, this method is applied to
the nearby galaxy NGC 628, jointly fitting the galactic disk and a small, crowded region at its center. This is a very
simple case, however, as the center region is small and contributes little to the overall IR SED.

This method could be extended to:

- Larger numbers of regions, to provide a computationally
  cheaper option to computing physical property maps.
- Separately modeling the nucleus and disk in nearby AGN.
- Determining the ages of individual star formation regions which are not separated in IR imaging?

## Requirements
- `Lightning`
- `Numpy/Scipy`
- `emcee`
- `astropy`

## Installing
Clone this repository to wherever you prefer, and then add

```
export PYTHONPATH="<install_dir>/multi_lightning/src/:$PYTHONPATH"
```

to your terminal startup file (e.g. `~/.bash_profile`). You will also need to have done this already for the main
`Lightning` repository.
