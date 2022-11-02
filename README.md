[![DOI](https://zenodo.org/badge/430704582.svg)](https://zenodo.org/badge/latestdoi/430704582)

# meltPT

**meltPT** is a Python package for whole-rock major-element themormobarometric analyses of basaltic (mafic) rocks. It contains modules for:
- Correcting sample compositions for effects of olivine crystallisation
- Estimating pressures and temperatures at which samples were last in equilibrium with the mantle
- Estimating melt fractions and potential temperatures on an individual basis or for a suite of samples

## Installation

### Linux/Mac

To avoid incompatibility issues with other packages you may have installed, 
we recommend running **meltPT** in a virtual environment. Start a terminal and
enter:

```
$ # Create virtual environment.
$ # Replace /path/to/virtual/environment with your desired path.
$ python3 -m venv /path/to/virtual/environment

$ # Activate new virtual environment.
$ source /path/to/virtual/environment/bin/activate
```

Note that, when using a virtual environment, some users may experience issues
trying to display images created with matplotlib (e.g., in our tutorials): the
default matplotlib backend is sometimes "agg", which cannot show figures. There
does not seem to be a universal solution to this issue, so please see
[here](https://matplotlib.org/3.1.3/faq/virtualenv_faq.html) for some options
on how to proceed. 

What you do next depends on whether you just want to use meltPT as is or if
you want to edit the source code.

#### Basic Usage

If you just want to use **meltPT** as it is you can now simply type:

```
(meltpt) $ pip install meltPT
```
  
This will take a minute or two and print lots of stuff. If it completes
successfully, you are now ready to use **meltPT**!

#### Development usage

If you want to edit **meltPT**'s source code, you first need to download it.
Navigate to your desired location and type:

```
(meltpt) $ git clone git@github.com:fmcnab/meltPT.git
```

You should now have a directory called "meltPT". To install the package, enter
this directory and run the setup script. Using the -e flag means that the code
will be installed in "editable" mode, and changes you make locally will be
incorporated without the need for a fresh install.

```
(meltpt) $ cd meltPT
(meltpt) $
(meltpt) $ # Standard install
(meltpt) $ pip install -e .
(meltpt) $
(meltpt) $ # Including extra dependencies - 
(meltpt) $ # required for those who wish to compile documentation locally.
(meltpt) $ pip install .[docs]
(meltpt) $ pip install -e .[docs] # (editable)
```

This will take a minute or two and print lots of stuff. If it completes
successfully, you are now ready to use **meltPT**!

### Windows

Unfortunately, we don't know much about using Python on a Windows system.
But we plan to find out!

## Using meltPT

Under ./Examples/Tutorials, you will find some exaple scripts to help you get
started.

For more information, the **meltPT** [ReadtheDocs](https://meltpt.readthedocs.io)
page has full documentation for the
[codebase](https://meltpt.readthedocs.io/en/latest/codedoc.html), 
a series of informative
[tutorials](https://meltpt.readthedocs.io/en/latest/tutorials.html), information
about [contributing](https://meltpt.readthedocs.io/en/latest/contributing.html)
and our [liscence](https://meltpt.readthedocs.io/en/latest/license.html).

## Citing meltPT

If you use our code, please cite us. Currently each **meltPT** release is
archived in our Zenodo repository:

*  M<sup>c</sup>Nab, F. and Ball, P. W. (2022), meltPT, version 1.0.1,
   *Zenodo*, [doi: 10.5281/zenodo.6948030](https://doi.org/10.5281/zenodo.6948030).

Please be sure to include the version number of the code you used, so that
others can reproduce your results. We are working on an accompanying
publication which we hope will be available soon.

We also urge you to cite the original literature on which our code is based.
The sample backtracking method is based on that of 
[Lee et al. (2009, EPSL)](https://doi.org/10.1016/j.epsl.2008.12.020). For
a list of themormobarometric schemes available in **meltPT**, and links to the
original papers, see the\
[Thermobarometers](https://meltpt.readthedocs.io/en/latest/thermobarometers.html)
section of our documentation.

If you use our melt-path fitting routines you should also cite:

*  pyMelt: [Matthews et al. (in rev., Volcanica)](https://doi.org/10.31223/X5JP7X)
*  melting model for our examples:
   [Katz et al. (2003, G-cubed)](https://doi.org/10.1029/2002GC000433)