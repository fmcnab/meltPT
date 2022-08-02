from setuptools import setup

setup(
    name='meltPT',
    version='v0.0.2-alpha',
    author = 'Fergus McNab, Paddy Ball',
    packages=['meltPT'],
    description=("Calculate pressures and temperatures of melting for basaltic rocks and fit with melting models."),
    install_requires=[
        'pandas',
        'numpy',
        'matplotlib',
        'scipy',
        'shapely',
        'pyMelt',
        'sphinx==4.4.0',
        'pyyaml',
        'sympy',
        'nbsphinx',
        'sphinx-rtd-theme'
    ],
    scripts=['meltPT/meltPT'],
)