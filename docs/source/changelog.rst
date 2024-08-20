Changelog
^^^^^^^^^

======
v1.2.0
======

No major changes. Creating neat new release to coincide with appearance of
M\ :sup:`c`\ Nab and Ball (2023, *Volcanica*, `doi: 10.30909/vol.06.01.6376
<https://doi.org/10.30909/vol.06.01.6376>`_), with citation information, etc.

======
v1.1.0
======

Bug fixes; class implementation for backtracking; new tutorials/docs.

* Lots of additions/updates to the Tutorials.

* Some basic instructions for installing on Windows.

* Some changes made to allow smooth running of the code with Binder.

* Previously, correcting sample compositions for fractional crystallisation of
  olivine was achieved using a series of stand-alone functions and was
  hard-coded into the Suite class. Now, we have collected these functions into
  a single class, "BacktrackOlivineFractionation", which is instantiated then
  passed to Suite. This change will facilitate future incorporation of
  alternative backtracking schemes.
  
* Some minor changes to the way bounding potential temperatures for a suite of
  pressure/temperature data are found.

======
v1.0.1
======

* The melting parameterisations used in **meltPT**, via **pyMelt**, use a
  quadratic function for the solidus, so that for large potential temperatures
  no melting path can be calculated. To avoid optimisation algorithms 
  attempting to test such high potential temperatures, we included an
  arbritrary upper bound on potential temperature. In this new version, rather
  than being set abritrarily, the upper bound is chosen based on the properties
  of the mantle object passed to the optimisation routine.
  
* A bug has been fixed which caused sample-composition backtracking with a
  constant value for the partition coefficient to fail in verbose mode.
  
* A discrepency between the script Examples/Tutorials/Tutorial_1.py and the
  version rendered in the documentation has been fixed.  

======
v1.0.0
======

Initial release.