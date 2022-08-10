Installation
^^^^^^^^^^^^

=========
Linux/Mac
=========

We recommend running **meltPT** in a virtual environment. Start a terminal and
enter:

.. code-block:: console

  $ # Create virtual environment.
  $ # Replace /path/to/virtual/environment with your desired path.
  $ python3 -m venv /path/to/virtual/environment

  $ # Activate new virtual environment.
  $ source /path/to/virtual/environment/bin/activate

Next, download the code. Navigate to your desired location and type:

.. code-block:: console

  $ git clone git@github.com:fmcnab/meltPT.git

You should now have a directory called "meltPT". To install the package, enter
this directory and run the setup script:

.. code-block:: console

  $ cd meltPT
  $ python3 setup.py install

This will take a minute or two and print lots of stuff. If it completes
successfully, you are now ready to use **meltPT**!

=======
Windows
=======

Unfortunately, we don't know much about using Python on a Windows system.
But we plan to find out!

============
Dependencies
============

**meltPT** requires various other packages to be installed in order to work
correctly. They are:

* pandas
* numpy
* scipy
* matplotlib
* shapely
* pyMelt
* pyyaml
* sympy

If you follow the steps above, these dependencies will automatically be
installed alongside **meltPT**. The additional packages required to compile
documentation locally are:

* sphinx, version 4.4.0
* nbsphinx
* sphinx-rtd-theme
* ipykernel