[![DOI](https://zenodo.org/badge/430704582.svg)](https://zenodo.org/badge/latestdoi/430704582)

# meltPT

**meltPT** is a Python package for whole-rock major-element themormobarometric analyses of basaltic (mafic) rocks. It contains modules for:
- Correcting sample compositions for effects of olivine crystallisation
- Estimating pressures and temperatures at which samples were last in equilibrium with the mantle
- Estimating melt fractions and potential temperatures on an individual basis or for a suite of samples

## Getting started

### Installation

We recommend running **meltPT** in a virtual environment. Start a terminal and enter:

```
# Create virtual environment.
# Replace /path/to/virtual/environment with your desired path.
python3 -m venv /path/to/virtual/environment

# Activate new virtual environment.
source /path/to/virtual/environment/bin/activate
```

Next, download the code. Navigate to your desired location and type:

```
git clone git@github.com:fmcnab/meltPT.git
```

You should now have a directory called "meltPT". To install the package, enter this directory and run the setup script:

```
cd meltPT
pip install .
```

This will take a minute or two and print lots of stuff. If it completes successfully, you are now ready to use **meltPT**!

For more information and a list of dependencies visit our [ReadtheDocs](https://meltpt.readthedocs.io/en/latest/installation.html) page.

### Running meltPT

We provide a command-line tool that will allow you start running the code without needing to write your own scripts. All you need to do is edit a parameter file that contains paths to your chosen input and output files, values for certain variables, and the analyses you would like to run.

As an example, we recreate a result provided by [Plank and Forsyth (2016, *G<sup>3</sup>*)](https://doi.org/10.1002/2015GC006205)  . From the base directory:

```
# Navigate to the Plank and Forsyth directory.
cd ./Examples/Plank_and_Forsyth_2016

# Run meltPT with the example parameter file provided.
meltPT parameters.yaml
```

The program should produce a file called "PF16_S7_out.csv", containing the sample's backtracked composition, estimated equilibration pressure and temperature, and potential temperature. See also an example parameter file in the "Examples/Hawaii" directory.

Elsewhere in the "Examples" directory you will find lots of scripts that utilise all the different functionalities of meltPT, and will hopefully help you write your own.

## Using meltPT

For more information, the **meltPT** [ReadtheDocs](https://meltpt.readthedocs.io) page has full documentation for the [codebase](https://meltpt.readthedocs.io/en/latest/codedoc.html), a series of informative [tutorials](https://meltpt.readthedocs.io/en/latest/tutorials.html), information about [contributing](https://meltpt.readthedocs.io/en/latest/contributing.html) and our [liscence](https://meltpt.readthedocs.io/en/latest/license.html).

## Citing meltPT

If you use our code, please cite us. We are currently preparing an accompnaying manuscript for publication. In the meantime, you can cite the package with its Zenodo doi.

We also urge you to cite the original literature on which our code is based:
- Sample backtracking and thermobarometry: [Lee et al. (2009, *EPSL*)](https://doi.org/10.1016/j.epsl.2008.12.020)
- Thermobarometry: please cite the appropriate thermobarometer [(see list)](https://meltpt.readthedocs.io/en/latest/thermobarometers.html).

If you use our melt-path fitting routines you should also cite:
- pyMelt: [Matthews et al. (in rev., *Volcanica*)](https://doi.org/10.31223/X5JP7X)
- melting model for our examples: [Katz et al. (2003, *G<sup>3</sup>*)](https://doi.org/10.1029/2002GC000433)
