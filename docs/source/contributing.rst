Contributing
^^^^^^^^^^^^

================
Getting involved
================

Contributions to **meltPT** are encouraged! The library has be designed to
support additional thermobarometers, melting models and fractional
crystallisation schemes. Please do not hesitate to contact Fergus McNab
(mcnab@gfz-potsdam.de) or Patrick Ball (paddywball@gmail.com) at any stage
in the editing process.

Any edits pushed to the GitHub page will be confirmed by the lead editors
(F. McNab and P. Ball) before they are integrated into the system. Additional
thermobarometers, melting models and fractional crystallisation schemes will
not be included without an accompanying benchmarking test.

We are always open to discussions about adding additional functionality to
**meltPT** -- just get in touch!


======================
Creating a new release
======================

Here are some notes for developers outlining the various steps in creating
a new **meltPT** release.

#.  *Update version numbers*

    At the moment version nnumbers are hard coded into the following files;
    remember to update them!
    
    * setup.py
    * README.md (citation)
    * docs/source/conf.py
    * docs/source/about.rst (citation)