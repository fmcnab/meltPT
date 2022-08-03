Thermobarometers
^^^^^

This list documents the thermobarometers and thermometers currently implemented within **meltPT**. To use each option, the abbreviated name at the beginning of each option can be typed within the method flag of compute_pressure_temperature(), compute_temperature() and/or compute_pressure(). The associated equation names in the list below refer to the equations in the original papers. See the code documentation for more details.

=============
Currently Available Options
=============

* P08       - `Putirka (2008, Revs. in Min. and Geo.) <https://doi.org/10.2138/rmg.2008.69.3>`_
              Simultaneous calculation of pressure and temperature using Equations 42 and 22, respectively.
* L09       - `Lee et al. (2009, EPSL) <https://doi.org/10.1016/j.epsl.2008.12.020>`_
              T and P(T) calculated using Equations 3 and 2, respectively. 
* TGK12_PLG - `Till et al. (2012, JGR: Solid Earth) <https://doi.org/10.1029/2011JB009044>`_
              P and T(P) calculated using Lines 1 and 2 of Table 5, respectively. Plagioclase must be stable mantle phase.
* TGK12_SPL - `Till et al. (2012, JGR: Solid Earth) <https://doi.org/10.1029/2011JB009044>`_
              P and T(P) calculated using Lines 3 and 4 of Table 5, respectively. Spinel must be stable mantle phase.
* PF16      - `Plank and Forsyth (2016, G-cubed) <https://doi.org/10.1002/2015GC006205>`_
              T and P(T) calculated using Equations 1 and 2, respectively. Corrected for CO2 using Equation 3.
* SD20      - `Sun and Dasgupta (2020, EPSL) <https://doi.org/10.1016/j.epsl.2020.116549>`_
              Simultaneous calculation of pressure and temperature using Equations 4 and 6, respectively.
* BK21      - `Brown Krien et al. (2021, JGR: Solid Earth) <https://doi.org/10.1029/2020JB021292>`_
              Selects and applies the best thermobarometer for a given sample composition from BK21_PLG, BK21_SPL, BK21_GNT.
* BK21_PLG  - `Brown Krien et al. (2021, JGR: Solid Earth) <https://doi.org/10.1029/2020JB021292>`_
              P and T(P) calculated using Equations 1P-P and 1P-T, respectively. Assumes plagioclase is stable mantle phase.
* BK21_SPL  - `Brown Krien et al. (2021, JGR: Solid Earth) <https://doi.org/10.1029/2020JB021292>`_
              P and T(P) calculated using Equations 1S-P and 1S-T, respectively. Assumes plagioclase is stable mantle phase.
* BK21_GNT  - `Brown Krien et al. (2021, JGR: Solid Earth) <https://doi.org/10.1029/2020JB021292>`_
              P and T(P) calculated using Equations 1G-P and 1G-T, respectively. Assumes plagioclase is stable mantle phase.

**Thermometer Only**

* B93       - `Beattie, (1993, Contrib. to Min. and Pet.) <https://doi.org/10.1007/BF00712982>`_
              Calculates temperate of olivine liquidus using Equation 10.
* P07_2     - `Putirka et al. (2007, Chemical Geology) <https://doi.org/10.1016/j.chemgeo.2007.01.014>`_
              Calculates temperate using Equation 2.
* P07_4     - `Putirka et al. (2007, Chemical Geology) <https://doi.org/10.1016/j.chemgeo.2007.01.0145>`_
              Calculates temperate using Equation 4.
* HA15      - `Herzberg and Asimow (2015, G-cubed) <https://doi.org/10.1002/2014GC005631>`_
              Calculate temperate using B93 converted to desired pressure using Equation 12.

**Unfinished Thermobarometers**

* G13       - Grove et al. (2013, Contrib. Min. and Pet.)  <https://doi.org/10.1007/s00410-013-0899-9>`_