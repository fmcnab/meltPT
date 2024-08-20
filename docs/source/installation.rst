Installation
^^^^^^^^^^^^

==========
Quickstart
==========

To try out **meltPT** without having to install anything, you can run our
Tutorials in your browser via our 
`Binder <https://mybinder.org/v2/gh/fmcnab/meltPT/master>`_ page. You will find 
our jupyter notebook tutorials under "./Examples/Tutorials/jupyter".
Alternatively, you can use the terminal function to run the tutorial scripts.
For example:

.. code-block:: console

  $ # Navigate to script directory.
  $ cd Examples/Tutorials/scripts

  $ # Run the first tutorial.
  $ python3 Tutorial_1.py

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

Note that, when using a virtual environment, some users may experience issues
trying to display images created with matplotlib (e.g., in our tutorials): the
default matplotlib backend is sometimes "agg", which cannot show figures. There
does not seem to be a universal solution to this issue, so please see
`here <https://matplotlib.org/3.1.3/faq/virtualenv_faq.html>`_ for some options
on how to proceed. 

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

  (meltpt) $ git clone git@github.com:fmcnab/meltPT.git

You should now have a directory called "meltPT". To install the package, enter
this directory and run the setup script. Using the -e flag means that the code
will be installed in "editable" mode, and changes you make locally will be
incorporated without the need for a fresh install.

.. code-block:: console

  (meltpt) $ cd meltPT
  (meltpt) $
  (meltpt) $ # Standard install
  (meltpt) $ pip install -e .
  (meltpt) $
  (meltpt) $ # Including extra dependencies - 
  (meltpt) $ # required for those who wish to compile documentation locally.
  (meltpt) $ pip install .[docs]
  (meltpt) $ pip install -e .[docs] # (editable)

As above, this will take a minute or two and print lots of stuff. If it
completes successfully, you should now be ready to go.

=======
Windows
=======

We do not have much experience running Python on a Windows system, but outline
here some basic steps that we have tested and believe should be accessible
to most users.

First, you need to install a Python package manager, if you don't have one
already. We tested the Miniforge3 installer from 
`Miniforge <https://github.com/conda-forge/miniforge>`_, which is free to
anyone regardless of affiliation. If you are already running a different Conda 
distribution, don't worry, the following steps should still work.

Once you have a package manager installed, open the program. You should see
a command prompt. Create a new environment in which to install **meltPT**,
called, for example, 'meltpt':

.. code-block:: console

  (base) > conda create meltpt
  
Next, activate the environment:

.. code-block:: console

  (base) > conda activate meltpt
  (meltpt) >
  
If you want to install an Interactive Development Environment (IDE), allowing
you to edit scripts, use a Python interface etc., now is the time to install
it. For example, `Spyder <https://www.spyder-ide.org/>`_:

.. code-block:: console

  (meltpt) > conda install spyder
  
Finally, install **meltPT** using pip:

.. code-block:: console

  (meltpt) > pip install meltPT

If you wish to edit the source code, see the Linux/Mac instructions above for 
alternative pip commands; these should also work in Conda. You should now be 
ready to use **meltPT**!

============
Dependencies
============

**meltPT** requires various other packages to be installed in order to work
correctly. We have tested meltPT with the following versions:

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

If you follow the steps above, these dependencies will be automatically
installed alongside **meltPT**, with approximately these versions. Note that,
if you don't use a virtual environment as described above, this might mean
that versions of common packages you have already installed may change, and
other packages you have installed may no longer work; this is why we
recommend virtual environments.