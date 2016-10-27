.. _volume:

.. currentmodule:: openpathsampling

Volume Functions
================

:class:`openpathsampling.Volume`

    >>> import openpathsampling as paths
    >>> volume = paths.Vnsemble()

basic volumes
-------------
.. autosummary::
    :toctree: api/generated/

    Volume
    FullVolume
    EmptyVolume

set-based volume combinations
-----------------------------
.. autosummary::
    :toctree: api/generated/

    IntersectionVolume
    UnionVolume
    SymmetricDifferenceVolume
    RelativeComplementVolume
    VolumeCombination

collective variable-based volumes
---------------------------------
.. autosummary::
    :toctree: api/generated/

    CVDefinedVolume
    PeriodicCVDefinedVolume

voronoi cell volumes
--------------------
.. autosummary::
    :toctree: api/generated/

    VoronoiVolume

volume factory
--------------
.. autosummary::
    :toctree: api/generated/

    VolumeFactory