Getting Started
^^^^^^^^^^^^^^^

============
Installation
============

---------
Linux/Mac
---------

We recommend running **meltPT** in a virtual environment. Start a terminal and enter:

.. code-block:: console

  $ # Create virtual environment.
  $ # Replace /path/to/virtual/environment with your desired path.
  $ python3 -m venv /path/to/virtual/environment

  $ # Activate new virtual environment.
  $ source /path/to/virtual/environment/bin/activate

Next, download the code. Navigate to your desired location and type:

.. code-block:: console

  $ git clone git@github.com:fmcnab/meltPT.git

You should now have a directory called "meltPT". To install the package, enter this directory and run the setup script:

.. code-block:: console

  $ cd meltPT
  $ python3 setup.py install

This will take a minute or two and print lots of stuff. If it completes successfully, you are now ready to use **meltPT**!

-------
Windows
-------

Unfortunately, we don't know much about using Python on a Windows system. But we will find out!

==============
Running meltPT
==============

We provide a command-line tool that will allow you start running the code without needing to write your own scripts. All you need to do is edit a parameter file that contains paths to your chosen input and output files, values for certain variables, and the analyses you would like to run.

As an example, we recreate a result provided by `Plank and Forsyth (2016-cubed) <https://doi.org/10.1002/2015GC006205>`_. From the base directory:

.. code-block:: console

  $ # Navigate to the Plank and Forsyth directory.
  $ cd ./Examples/Plank_and_Forsyth_2016

  $ # Run meltPT with the example parameter file provided.
  $ meltPT parameters.yaml

The program should produce a file called "PF16_S7_out.csv", containing the sample's backtracked composition, estimated equilibration pressure and temperature, and potential temperature. See also an example parameter file in the "Examples/Hawaii" directory.

Elsewhere in the "Examples" directory you will find lots of scripts that utilise all the different functionalities of **meltPT**, and will hopefully help you write your own.
