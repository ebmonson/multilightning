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


.. toctree::
   :maxdepth: 1
   :caption: Documentation

   recipe
   examples/NGC628.ipynb
   MultiLightning
   Priors

Attribution
-----------
``multilightning`` was originally used in Lehmer et al. 2024::

    @article{lehmer2024
             TKTKTKTKTKTK
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
