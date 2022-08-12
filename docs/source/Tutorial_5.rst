=================================
Tutorial 5: The Command-Line Tool
=================================

We provide a command-line tool that will allow you start running the code
without needing to write your own scripts. All you need to do is edit a
parameter file that contains paths to your chosen input and output files, 
values for certain variables, and the analyses you would like to run.

As an example, we repeat the analysis shown in Tutorial 1
`Plank and Forsyth (2016-cubed) <https://doi.org/10.1002/2015GC006205>`_.
Under "./Examples/Tutorials" you should find an example
parameter file called "Tutorial_5.yaml" which contains:

.. code-block:: python

  # ---------------------------------------------------------------------------- #
  # -------------------- PARAMETER FILE FOR USE WITH meltPT -------------------- #
  #                                                                              #
  # To run, at the command line type:                                            #
  #                                                                              #
  # $ meltPT <parameter_file_name>.yaml                                          #
  #                                                                              #
  # ---------------------------------------------------------------------------- #


  # ---- Input / output filenames and options ---------------------------------- #

  # Path to input data file.
  data_file: "Data/PF16_UT09DV04.csv"

  # Path to output data file.
  output_file: "Data/PF16_UT09DV04_out.csv"

  # How to read input CSV.
  read_as_primary: False
  read_PT: False

  # ---- Parameters for backtracking compositions ------------------------------ #

  backtrack:
    
    apply: True
    
    # Ration of Ce to H2O in the melt.
    # Only used if H2O not provided.
    Ce_to_H2O: 200.
    
    # Whether to apply Sun & Dasgupta (2022) SiO2 --> CO2 parameterisation.
    param_co2: False
    
    # Ratio of Fe3+ to total Fe in the melt. 
    src_FeIII_totFe: 0.17

    # Minimum MgO for a sample to be processed in wt%. Should be above ~8.
    min_MgO: 8.

    # Target Forsterite number for samples during backtracking.
    target_Fo: 0.9
    
    # Partition coefficient for Mg and Fe exchange between olivine and melt.
    # Can be a float (e.g. ~0.3) or if False will be calculated as a function
    # of olivine Forsterite number.
    Kd: 0.3
    
    # Fraction of olivine to add at each iteration.
    dm: 0.0005
    
    # Maximum fraction of olivine to add before giving up on backtracking.
    max_olivine_addition: 0.3
    
    # Save final primary compositions.
    save: True


  # ---- Options for computing pressure and temperature of melting ------------- #

  PT:
    
    apply: True
    save: True


  # ---- Options for fitting melting paths ------------------------------------- #

  # Fitting a melting path to the entire suite.
  # Each option is True / False.
  suite_Tp:
    apply: False
    plot: False
    save: False
    
  # Fitting a melting path to each individual sample.
  # Note: with a lot of samples will be very slow!
  # Each option is True / False.
  individual_Tp:
    apply: True
    plot: True
    save: True

To run the program, simply type:

.. code-block:: console

  $ meltPT Tutorial_5.yaml

The program should produce a file called "PF16_UT09DV04_out.csv", containing
the sample's backtracked composition, estimated equilibration pressure and
temperature, and potential temperature.

To perform your own analyses, simply edit the various parameters in the
parameters.yaml file as you wish. The meltPT command should work from
anywhere in your system. (If you install meltPT in a virtual environment, 
just make sure the environment activated.)