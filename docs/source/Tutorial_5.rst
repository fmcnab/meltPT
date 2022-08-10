=================================
Tutorial 5: The Command-Line Tool
=================================

We provide a command-line tool that will allow you start running the code
without needing to write your own scripts. All you need to do is edit a
parameter file that contains paths to your chosen input and output files, 
values for certain variables, and the analyses you would like to run.

As an example, we repeat the analysis shown in Tutorial 1
`Plank and Forsyth (2016-cubed) <https://doi.org/10.1002/2015GC006205>`_. 
From the base directory:

.. code-block:: console

  $ # Navigate to Tutorial 5's directory.
  $ cd ./Examples/Tutorials/Tutorial_5

  $ # Run meltPT with the example parameter file provided.
  $ meltPT parameters.yaml

The program should produce a file called "PF16_UT09DV04_out.csv", containing
the sample's backtracked composition, estimated equilibration pressure and
temperature, and potential temperature.

To perform your own analyses, simply edit the various parameters in the
parameters.yaml file as you wish. The meltPT command should work from
anywhere in your system. (If you install meltPT in a virtual environment, 
just make sure the environment activated.)