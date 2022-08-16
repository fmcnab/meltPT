Installation
^^^^^^^^^^^^

=========
Linux/Mac
=========

To avoid incompatibility issues with other packages you may have installed, 
we recommend running **meltPT** in a virtual environment. To create an
environment called "meltpt", start a terminal and enter:

.. code-block:: console

  $ # Create virtual environment.
  $ # Replace /path/to/virtual/environment with your desired path.
  $ python3 -m venv /path/to/virtual/environment/meltpt

  $ # Activate new virtual environment.
  $ source /path/to/virtual/environment/meltpt/bin/activate

What you do next depends on whether you just want to use meltPT as is or if
you want to edit the source code.

-----------
Basic usage
-----------

If you just want to use **meltPT** as it is you can now simply type:

.. code-block:: console

  (meltpt) $ pip install meltPT
  
This will take a minute or two and print lots of stuff. If it completes
successfully, you are now ready to use **meltPT**!

-----------------
Development usage
-----------------

If you want to edit **meltPT**'s source code, you first need to download it.
Navigate to your desired location and type:

.. code-block:: console

  $ git clone git@github.com:fmcnab/meltPT.git

You should now have a directory called "meltPT". To install the package, enter
this directory and run the setup script. Using the -e flag means that the code
will be installed in "editable" mode, and changes you make locally will be
incorporated without the need for a fresh install.

.. code-block:: console

  $ cd meltPT
  $
  $ # Standard install
  $ pip install -e .
  $
  $ # Including extra dependencies - 
  $ # required for those who wish to compile documentation locally.
  $ pip install .[docs]
  $ pip install -e .[docs] # (editable)

As above, this will take a minute or two and print lots of stuff. If it
completes successfully, you should now be ready to go.

=======
Windows
=======

Unfortunately, we don't know much about using Python on a Windows system.
But we plan to find out!

============
Dependencies
============

**meltPT** requires various other packages to be installed in order to work
correctly. We have tested meltPT with the following versions. They are:

* pandas, v1.4
* numpy, v1.23
* matplotlib, v3.5
* scipy, v1.8
* shapely, v1.8
* pyMelt, v1.96
* pyyaml, v6.0
* sympy, v1.10

The additional packages required to compile documentation locally are:

* sphinx, v5.1
* nbsphinx, v0.8
* sphinx-rtd-theme, v1.0
* ipykernel, v6.15

If you follow the steps above, these dependencies will automatically be
installed alongside **meltPT**, with approximately these versions. Note that,
if you don't use a virtual environment as described above, this might mean
that versions of common packages you have already installed may change, and
other packages you have installed may no longer work; this is why we
recommend virtual environments.