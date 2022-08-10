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
        'pyyaml',
        'sympy'
    ],
    extras_require={'docs': [
            'sphinx==4.4.0',
            'nbsphinx',
            'sphinx-rtd-theme',
            'ipykernel'
        ]
    },
    scripts=['meltPT/meltPT'],
)