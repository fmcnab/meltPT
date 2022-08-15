from setuptools import setup

setup(
    name='meltPT',
    version='v0.0.7-alpha',
    author = 'Fergus McNab, Paddy Ball',
    packages=['meltPT'],
    description=("Calculate pressures and temperatures of melting for basaltic rocks and fit with melting models."),
    install_requires=[
        'pandas~=1.4',
        'numpy~=1.23',
        'matplotlib~=3.5',
        'scipy~=1.8',
        'shapely~=1.8',
        'pyMelt~=1.96',
        'pyyaml~=6.0',
        'sympy~=1.10'
    ],
    extras_require={'docs': [
            'sphinx~=5.1',
            'nbsphinx~=0.8',
            'sphinx-rtd-theme~=1.0',
            'ipykernel~=6.15'
        ]
    },
    scripts=['meltPT/meltPT'],
)