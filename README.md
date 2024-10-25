# MultiLightning
---

[![Documentation Status](https://readthedocs.org/projects/multilightning/badge/?version=latest)](https://multilightning.readthedocs.io/en/latest/?badge=latest)

A package for jointly fitting multiple regions of a galaxy with the python version of the `lightning` SED-fitting code.

For nearby galaxies, we often have a subsets of bands in which the galaxy is resolved and unresolved (e.g., the far-IR).
We can thus choose to either discard the spatial information of the resolved bands and fit the integrated light of the
galaxy, or fit without the unresolved data, potentially losing valuable constraints on the star formation rate.

Here, we define a joint likelihood for fitting multiple resolved regions, under the assumption that the unresolved
photometry is described by the sum of the models for each resolved region. By way of example, this method is applied to
the nearby galaxy NGC 628, jointly fitting the galactic disk and a small, crowded region at its center. This is a very
simple case, however, as the center region is small and contributes little to the overall IR SED.

This method could be extended to, hypothetically:

- Larger numbers of regions, to provide a computationally
  cheaper option to computing physical property maps.
- Separately modeling the nucleus and disk in nearby AGN.
- Determining the ages of individual star formation regions which are not separated in IR imaging.

## Requirements
Note that the requirements of `lightning` are implicit, and that `multilightning` has no required
dependencies beyond those of `lightning`. The below are the explicit dependencies:

- `lightning`
- `numpy`
- `scipy`
- `emcee`

## Installing

Before installing this package, you should install `lightning` by following the instructions in its documentation, optionally
creating a new virtual environment in e.g. `conda` to install it in.

Then run

```sh
git clone https://github.com/ebmonson/multilightning.git
cd multilightning

pip install .  --no-deps

```

to clone this repository and install the code locally.

> [!Note]
> I may not go to the extra step of uploading this package to pypi/conda-forge since it is a niche, hacky
> thing on top of `lightning`.

## Documentation
---
Online documentation for the package can be found on [its readthedocs page](https://multilightning.readthedocs.io/en/latest/),
and a compiled PDF version of the documentation is available [in this repository](https://github.com/ebmonson/multilightning/blob/main/docs/multilightning.pdf).
