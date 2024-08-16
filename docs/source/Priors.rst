Priors
------

``multilightning`` is compatible with the priors defined in ``lightning.priors``. However, the addition of multi-region
fitting introduces the need for priors which *connect parameters together* between different models. ``multilightning``
implements three connections: Normal and Uniform connections, which assume the distance and absolute distance in the
parameter values are normally and uniformly distributed, respectively; and a fixed connection, which pins the parameter
to a reference value.

These are slightly hacked together, which you can tell by the weird way you initialize them. Each prior takes as arguments:

- Its own parameters: :math:`\mu` and :math:`\sigma` for the normal connection, :math:`a` and :math:`b` for the uniform,
  nothing for the fixed connection.
- The name of the region to target, as a string.
- The index of the parameter to target.

Therefore ``NormalConnection([0.0, 1.0], 'disk', 0)`` assumes that the distance between the given parameter and the first
parameter (i.e., the parameter at index 0) of the ``'disk'`` component is drawn from the standard normal distribution.

.. note::

    This construction actually means you can use these priors to link parameters within the same region: the above
    construction could be used to link the second piecewise-constant SFH coefficient in the disk component to the first,
    for example.

.. autoclass:: multi_lightning.priors.NormalConnection
    :show-inheritance:
    :members: evaluate, quantile, sample

.. autoclass:: multi_lightning.priors.UniformConnection
    :show-inheritance:
    :members: evaluate, quantile, sample

.. autoclass:: multi_lightning.priors.FixedConnection
    :members: evaluate
