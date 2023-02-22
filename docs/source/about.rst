About
^^^^^

==========
Background
==========

**meltPT** is a Python package for whole-rock major-element themormobarometric
analyses of basaltic (mafic) rocks. It contains modules for:

*  Correcting sample compositions for effects of olivine crystallisation
*  Estimating pressures and temperatures at which samples were last in
   equilibrium with the mantle
*  Estimating melt fractions and potential temperatures on an individual basis
   or for a suite of samples

You can find more background information and some example analyses in our paper
in *Volcanica* `here <https://doi.org/10.30909/vol.06.01.6376>`_.

=============
Citing meltPT
=============

If you use our code, please cite us!

* M\ :sup:`c`\ Nab, F. and Ball, P. W. (2023), ``meltPT``: A ``Python`` package
  for basaltic whole-rock thermobarometric analysis with application to Hawai'i,
  *Volcanica*, 6(1), p. 63--76, `doi: 10.30909/vol.06.01.6376 <https://doi.org/10.30909/vol.06.01.6376>`_.

You should also refer to the specific release of the code you used. For example,
the most recent **meltPT** release is archived in our Zenodo repository:

*  M\ :sup:`c`\ Nab, F. and Ball, P. W. (2023), meltPT, version 1.2.0,
   *Zenodo*, `doi: 10.5281/zenodo.6948030 <https://doi.org/10.5281/zenodo.6948030>`_.

We also urge you to cite the original literature on which our code is based.
The sample backtracking method is based on that of 
`Lee et al. (2009, EPSL) <https://doi.org/10.1016/j.epsl.2008.12.020>`_. For
a list of themormobarometric schemes available in **meltPT**, and links to the
original papers, see the :doc:`thermobarometers` section.

If you use our melt-path fitting routines you should also cite:

*  pyMelt: `Matthews et al. (2022, Volcanica) <https://doi.org/10.30909/vol.05.02.469475>`_
*  melting model for our examples:
   `Katz et al. (2003, G-cubed) <https://doi.org/10.1029/2002GC000433>`_