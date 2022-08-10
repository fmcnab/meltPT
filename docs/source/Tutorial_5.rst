==================================
Tutorial 5 - The Command-Line Tool
==================================

We provide a command-line tool that will allow you start running the code
without needing to write your own scripts. All you need to do is edit a
parameter file that contains paths to your chosen input and output files, 
values for certain variables, and the analyses you would like to run.

As an example, we recreate a result provided by 
`Plank and Forsyth (2016-cubed) <https://doi.org/10.1002/2015GC006205>`_. 
From the base directory:

.. code-block:: console

  $ # Navigate to the Plank and Forsyth directory.
  $ cd ./Examples/Plank_and_Forsyth_2016

  $ # Run meltPT with the example parameter file provided.
  $ meltPT parameters.yaml

The program should produce a file called "PF16_S7_out.csv", containing the
sample's backtracked composition, estimated equilibration pressure and
temperature, and potential temperature. See also an example parameter file
in the "Examples/Hawaii" directory.
