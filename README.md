# meltPT

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
python3 setup.py install
```

You are now ready to use **meltPT**!

### Running meltPT

We provide a command-line tool that will allow you start running the code without needing to write your own scripts. All you need to do is edit a parameter file that contains paths to your chosen input and output files, values for certain variables, and the analyses you would like to run.

As an example, we recreate a result provided by Plank and Forsyth (2016, G<sup>3</sup>). From the base directory:

```
# Navigate to the Plank and Forsyth directory.
cd ./Examples/Plank_and_Forsyth_2016

# Run meltPT with the example parameter file provided.
meltPT parameters.yaml
```

The program should produce a file called "PF16_S7_out.csv", containing the sample's backtracked composition, estimated equilibration pressure and temperature, and potential temperature.

Elsewhere in the "Examples" directory you will find lots of scripts that utilise all the different functionalities of meltPT, and will hopefully help you write your own.