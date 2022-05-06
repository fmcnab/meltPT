from setuptools import setup

setup(
    name='meltPT',
    version='0.0.1',
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
        'PyQt5',
        'sphinx',
        'pyyaml',
    ],
    scripts=['meltPT/meltPT'],
)