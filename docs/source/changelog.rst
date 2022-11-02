Changelog
^^^^^^^^^

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