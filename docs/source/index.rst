.. multilightning documentation master file, created by
   sphinx-quickstart on Thu Aug 15 14:04:15 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

multilightning
============================

``multilightning`` provides an object oriented interface to jointly fit spectral energy distributions of multiple
regions of a galaxy, using the python version of the ``Lightning`` SED-fitting code.

For nearby galaxies, we often have a subsets of bands in which the galaxy is resolved and unresolved (e.g., the far-IR).
We can thus choose to either discard the spatial information of the resolved bands and fit the integrated light of the
galaxy, or fit without the unresolved data, potentially losing valuable constraints on the star formation rate.

Here, we define a joint likelihood for fitting multiple resolved regions, under the assumption that the unresolved
photometry is described by the sum of the models for each resolved region. By way of example, this method is applied to
the nearby galaxy NGC 628, jointly fitting the galactic disk and a small, crowded region at its center. This is a very
simple case, however, as the center region is small and contributes little to the overall IR SED.

This method could be extended to:

- Larger numbers of regions, to provide a computationally cheaper option to computing physical property maps.
- Separately modeling the nucleus and disk in nearby AGN.
- Determining the ages of individual star formation regions which are not separated in IR imaging?


Installation
------------

Before installing this package, you should install `lightning` by following the instructions in its documentation, optionally
creating a new virtual environment in e.g. `conda` to install it in.

Then run::

    git clone https://github.com/ebmonson/multilightning.git
    cd multilightning

    pip install . --no-deps


to clone this repository and install the code locally.

.. note::

    I may not go to the extra step of uploading this package to pypi/conda-forge since it is a niche, hacky
    thing on top of ``lightning``.


Documentation
-------------

.. toctree::
   :maxdepth: 1

   recipe
   examples/NGC628.ipynb
   MultiLightning
   Priors


Attribution
-----------
``multilightning`` was originally used in `Lehmer et al. (2024, in press at ApJS) <https://ui.adsabs.harvard.edu/abs/2024arXiv241019901L/abstract>`_::

    @ARTICLE{2024arXiv241019901L,
        author = {{Lehmer}, Bret D. and {Monson}, Erik B. and {Eufrasio}, Rafael T. and {Amiri}, Amirnezam and {Doore}, Keith and {Basu-Zych}, Antara and {Garofali}, Kristen and {Oskinova}, Lidia and {Andrews}, Jeff J. and {Antoniou}, Vallia and {Geda}, Robel and {Greene}, Jenny E. and {Kovlakas}, Konstantinos and {Lazzarini}, Margaret and {Richardson}, Chris T.},
        title = "{An Empirical Framework Characterizing the Metallicity and Star-Formation History Dependence of X-ray Binary Population Formation and Emission in Galaxies}",
        journal = {arXiv e-prints},
        keywords = {Astrophysics - Astrophysics of Galaxies, Astrophysics - High Energy Astrophysical Phenomena},
        year = 2024,
        month = oct,
        eid = {arXiv:2410.19901},
        pages = {arXiv:2410.19901},
        archivePrefix = {arXiv},
        eprint = {2410.19901},
        primaryClass = {astro-ph.GA},
        adsurl = {https://ui.adsabs.harvard.edu/abs/2024arXiv241019901L},
        adsnote = {Provided by the SAO/NASA Astrophysics Data System}
    }

and can be cited using its ASCL identifier::

    @software{multilightning
              TKTKTKTKTKTK
             }


License
-------
``multilightning`` is released under the terms of the MIT license.


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
