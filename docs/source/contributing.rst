.. highlight:: console

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

    At the moment version nnumbers are hard coded into the following files:
    
    * setup.py
    * README.md (citation)
    * docs/source/conf.py
    * docs/source/about.rst (citation)
    
    Please remember to update them!

#.  *Create GitHub release*

    Create a new release on GitHub and check it worked properly on
    `Zenodo <https://doi.org/10.5281/zenodo.6948030>`_.

#.  *Build*

    Make sure build is up to date. Then run build.

    From the base directory::

        $ python3 -m pip install --upgrade build
        $ python3 -m build

#.  *Upload to TestPyPI*

    You will need an account on TestPyPI for this. Upgrade twine, then,
    being sure to include the right version number, upload::

        $ python3 -m pip install --upgrade twine
        $ python3 -m twine upload --repository testpypi dist/meltpt-VERSION_NUMBER*

#.  *Try out TestPyPI release*

    Create a new virtual environment, try installing::

        $ python3 -m venv test_env
        $ source test/bin/activate
        (test) $ python3 -m pip install \
                    --index-url https://test.pypi.org/simple/ \
                    --extra-index-url https://pypi.org/simple \
                    meltPT

    Make sure the expected version was installed and everything works
    (e.g., run the tutorials).

#.  *Upload to the real PyPI*

    As before::
    
        $ python3 -m twine upload dist/meltpt-VERSION_NUMBER*

#.  *Test out PyPI release*

    As before::

        $ python3 -m venv test_env
        $ source test/bin/activate
        (test) $ python3 -m pip install meltPT
        
    Make sure the expected version was installed and everything works
    (e.g., run the tutorials).