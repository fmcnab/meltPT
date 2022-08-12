==============
The Input File
==============

All analyses in **meltPT** require an input file containing the compositions of
samples to be analaysed in the form of a simple csv. In the source code you
will find many examples -- here we briefly describe what the file must contain,
what it may contain, and what will be done with different types of information.

------------------------
Major elements: Required
------------------------

An important step in **meltPT** analyses is backtracking fractional crystallisation
by olivine to estimate sample primary compositions (i.e., melt compositions
when they were last in equilibrium with the mantle). We assume that
crystallising olivine consists of SiO\ :sub:`2`\ , MgO and FeO. We need to
estimate the composition of olivine that would be in equilibrium with our melt,
which is a function of Mg and Fe content. For the these reasons, we require
input files to include columns with the headings: "SiO2", "MgO" and "FeO"; if
these columns are not present, the program will crash.

------------------------
Major elements: Optional
------------------------

Each of the thermobarometers we inlude in **meltPT** parameterise pressure and
temperature as functions of different major elements. As such, we will look
for the following input columns:
  - Al2O3
  - Fe2O3 & FeO_tot (i.e., combined FeO and Fe2O3)
  - CaO
  - Na2O
  - K2O
  - TiO2
  - MnO
  - Cr2O3
  - P2O5
  - NiO
  - CoO
  - H2O
  - CO2
Though none of these are required.

If any of these columns are not present, they will generally be assigned zero
values. This procedure is necessary for the smooth running of the code.
However, it will lead to inaccurate results if you apply a thermobarometer that
relies on an element for which you do not provide a concentration. We therefore
urge users to familiarise themselves with the assumptions behind any
thermobarometer they use, and carefully consider whether it is appropriate to
apply a specific thermobarometer to their data set.

-------------------------------------------------
Trace elements and major-element parameterisation
-------------------------------------------------

Trace-element concentrations are not generally required for analyses in
**meltPT** but will be used in some specific circumstances.

Concentrations of volatile phases such as H\ :sub:`2`\ O and CO\ :sub:`2` are
notoriously difficult to measure and often not available. An alternative is to
use a calibrated proxy to estimate their concentrations. In meltPT, if 
H\ :sub:`2`\ O concentrations are set to zero, we search for a Ce column, 
expected in ppm, and estimate water content using a specified
H\ :sub:`2`\ O\ /Ce value. For CO\ :sub:`2` we optionally use a parameterisation
based on SiO\ :sub:`2` content 
`(Sun and Dasgupta, 2020, EPSL) <https://doi.org/10.1016/j.epsl.2020.116549>`_.

Finally, since NiO, CoO, and Cr\ :sub:`2`\ O\ :sub:`3` are not always treated
as major elements, if they are not present, we will search for columns with
headings Ni, Co and Cr. Assuming these columns contain concentrations in ppm,
they are then used to calculate oxide concentrations in percent.

--------------------
Primary compositions
--------------------

We acknowledge that alternative methods exist to correct sample
compositions for the effects of fractionaly crystallisation of olivine (and
other phases). We therefore include an option to read the contents of the 
input file as primary compositions. Column headings should be as above, but
they will be treated by meltPT as primary compositions. See
:doc:`Tutorial_4` for examples.

-----------------
Other information
-----------------

Any other information (e.g. sample metadata, other elemental or isotopic
concentrations), will be read in and stored, maintaining the given column 
headings. These data will not be used by meltPT's functions but will be
readily available to facilitate easy comparision with thermobarometric results.