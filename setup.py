from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='meltPT',
    version='v1.0.1',
    author = 'F. McNab, P. W. Ball',
    author_email="mcnab@gfz-potsdam.de",
    url="https://github.com/fmcnab/meltPT",
    long_description=long_description,
    long_description_content_type="text/markdown",
    description=("Calculate pressures and temperatures of melting for basaltic rocks and fit with melting models."),
    packages=['meltPT'],
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
    python_requires='>=3.8'
)