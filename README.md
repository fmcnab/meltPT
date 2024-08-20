[![DOI](https://zenodo.org/badge/430704582.svg)](https://zenodo.org/badge/latestdoi/430704582)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/fmcnab/meltPT/master)

# meltPT

**meltPT** is a Python package for whole-rock major-element themormobarometric 
analyses of basaltic (mafic) rocks. It contains modules for:
- Correcting sample compositions for effects of olivine crystallisation
- Estimating pressures and temperatures at which samples were last in 
  equilibrium with the mantle
- Estimating melt fractions and potential temperatures on an individual basis
  or for a suite of samples
  
Full documentation can be found [here](https://meltpt.readthedocs.io/en/latest/).

You can find more background information and some example analyses in our paper
in *Volcanica* [here](https://doi.org/10.30909/vol.06.01.6376).

## Installation

### Quickstart

To try out **meltPT** without having to install anything, you can run our
Tutorials in your browser via our 
[Binder](https://mybinder.org/v2/gh/fmcnab/meltPT/master) page. You will find 
our jupyter notebook tutorials under "./Tutorials/jupyter".
Alternatively, you can use the terminal function to run the tutorial scripts.
For example:

```
$ # Navigate to script directory.
$ cd Tutorials/scripts

$ # Run the first tutorial.
$ python3 Tutorial_1.py
```

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

We do not have much experience running Python on a Windows system, but outline
here some basic steps that we have tested and believe should be accessible
to most users.

First, you need to install a Python package manager, if you don't have one
already. We tested the Miniforge3 installer from 
[Miniforge](https://github.com/conda-forge/miniforge>), which is free to
anyone regardless of affiliation. If you are already running a different Conda 
distribution, don't worry, the following steps should still work.

Once you have a package manager installed, open the program. You should see
a command prompt. Create a new environment in which to install **meltPT**,
called, for example, 'meltpt':

```
(base) > conda create meltpt
```

Next, activate the environment:

```
(base) > conda activate meltpt
(meltpt) >
```
  
If you want to install an Interactive Development Environment (IDE), allowing
you to edit scripts, use a Python interface etc., now is the time to install
it. For example, [Spyder](https://www.spyder-ide.org/>):

```
(meltpt) > conda install spyder
```
  
Finally, install **meltPT** using pip:

```
(meltpt) > pip install meltPT
```

If you wish to edit the source code, see the Linux/Mac instructions above for 
alternative pip commands; these should also work in Conda. You should now be 
ready to use **meltPT**!

## Using meltPT

Under ./Tutorials, you will find some exaple scripts to help you get
started.

For more information, the **meltPT** [ReadtheDocs](https://meltpt.readthedocs.io)
page has full documentation for the
[codebase](https://meltpt.readthedocs.io/en/latest/codedoc.html), 
a series of informative
[tutorials](https://meltpt.readthedocs.io/en/latest/tutorials.html), information
about [contributing](https://meltpt.readthedocs.io/en/latest/contributing.html)
and our [liscence](https://meltpt.readthedocs.io/en/latest/license.html).

## Citing meltPT

If you use our code, please cite us!

* M<sup>c</sup>Nab, F. and Ball, P. W. (2023), ``meltPT``: A ``Python`` package
  for basaltic whole-rock thermobarometric analysis with application to Hawai'i,
  *Volcanica*, 6(1), p. 63&ndash;76, [doi: 10.30909/vol.06.01.6376](https://doi.org/10.30909/vol.06.01.6376).

You should also refer to the specific release of the code you used. For example,
the most recent **meltPT** release is archived in our Zenodo repository:

*  M<sup>c</sup>Nab, F. and Ball, P. W. (2023), meltPT, version 1.2.0,
   *Zenodo*, [doi: 10.5281/zenodo.6948030](https://doi.org/10.5281/zenodo.6948030).

We also urge you to cite the original literature on which our code is based.
The sample backtracking method is based on that of 
[Lee et al. (2009, EPSL)](https://doi.org/10.1016/j.epsl.2008.12.020). For
a list of themormobarometric schemes available in **meltPT**, and links to the
original papers, see the\
[Thermobarometers](https://meltpt.readthedocs.io/en/latest/thermobarometers.html)
section of our documentation.

If you use our melt-path fitting routines you should also cite:

*  pyMelt: [Matthews et al. (2022, Volcanica)](https://doi.org/10.30909/vol.05.02.469475)
*  melting model for our examples:
   [Katz et al. (2003, G-cubed)](https://doi.org/10.1029/2002GC000433)