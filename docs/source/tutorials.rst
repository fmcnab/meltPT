Tutorials
^^^^^^^^^

In the following tutorials we explore **meltPT**'s range of functionality. All
analyses in meltPT require an input file containing compositions of samples to
be analysed. We first describe what this file must/can contain, and it how it
will be dealt with by **meltPT**.

In Tutorials 1-4 we go through example Python scripts that utilise meltPT
various capabilities. If you are used to working in Python, we hope neltPT's
various classes and functions will be straightforward to get used to and you
will be writing your own scripts in no time.

In Tutorial 5, we show how to use **meltPT**'s built in command-line tool. This
mode of working is less flexible but does not require you to write your own
Python code.

In Tutorial 6, we outline how to implement your own thermobarometer in the
context of **meltPT**.

All scripts and data corresponding to these Tutorials can be found on our
`GitHub page <https://github.com/fmcnab/meltPT/tree/master/Examples/Tutorials>`_
or in the `Zenodo archive <https://doi.org/10.5281/zenodo.6948030>`_. Moreover, 
the tutorials can also be downloaded as Jupiter Notebook files `here
<https://github.com/fmcnab/meltPT/tree/master/docs/source>`_.

Note that, if you installed **meltPT** in a virtual environment, you might have
some problems displaying images created with matplotlib: the default matplotlib
backend is sometimes "agg", which cannot show figures. There does not seem to
be a universal solution to this issue, so please see
`here <https://matplotlib.org/3.1.3/faq/virtualenv_faq.html>`_ for some options
on how to proceed. 

.. toctree::
   :maxdepth: 2
   :caption: Contents:
   
   input_file
   Tutorial_1
   Tutorial_2
   Tutorial_3
   Tutorial_4
   Tutorial_5
   Tutorial_6