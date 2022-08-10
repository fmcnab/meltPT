"""
=====
parse
=====

Read data from a csv file.

"""

import warnings

import numpy as np
import pandas as pd

def parse_csv(infile, Ce_to_H2O=0., src_FeIII_totFe=0., min_MgO=0., param_co2=False):
    """
    Read a csv and return a dataframe after some processing.
    
    Processing steps are:
    - check SiO2, MgO, and FexOx are specified; if not will crash
    - check other major elements are specified; if not will be set to zero
    - try to set some values via trace elements
    - redistribute Fe according to src_FeIII_totFe
    - if desired, estimate CO2 from SiO2
    - normalise major elements to 100%
    - reject samples with MgO less than some threshold.
      
    Parameters
    ----------
    infile : str
        Path to a csv containing data to be read.
    Ce_to_H2O : float
        Ratio of Ce to H2O in mantle source.
    src_FeIII_totFe : float
        Ratio of Fe3+ to total Fe in the mantle source.
    min_MgO : float, optional
        Minimum amount of MgO in sample to be accepted.
    read_as_primary : bool
        If true, data from input_csv is assumed to be primary
        and backtracking is avoided.
    param_co2 : bool
        If true, CO2 is calculated from SiO2 concentration.
        
    Returns
    -------
    df : pandas dataframe
        Dataframe containing processed data.
    """

    # Read in file
    df = pd.read_csv(infile, delimiter=",")
    
    # Check for compulsory columns
    compulsory_cols = ['SiO2', 'MgO']
    for col in compulsory_cols:
        if col not in df.columns:
            raise Exception("Input csv must contain a %s column." % col)
    if ("FeO" not in df.columns and 
        "Fe2O3" not in df.columns and 
        "FeO_tot" not in df.columns):
        raise Exception("Input csv must contain one of FeO, Fe2O3, FeO_tot columns.")

    # Replace empties and NaNs with zeros
    df = df.replace(r'^\s*$', 0., regex=True)
    df = df.replace(np.nan, 0., regex=True)

    # If major element columns do not exist in the input make them and
    # give them zeros for every row.
    major_columns = ['SiO2', 'TiO2', 'Al2O3', 'FeO', 'Fe2O3', 'MgO', 'CaO', 
                     'Na2O', 'K2O', 'MnO', 'Cr2O3', 'P2O5', 'NiO', 'CoO',
                     'H2O', 'CO2']
    for col in major_columns + ["FeO_tot"]:
        if col not in df.columns:
            message = "Input csv does not contain a %s column: we will try to fill it for you, or set it to zero." % col
            warnings.warn(message)
            df[col] = 0.
            
    # If not given compute some major oxides from trace element concentrations
    # df.loc[df['P2O5'] == 0, 'P2O5'] = df['P'] * 141.942524 / 2. / 10000.
    try:
        df.loc[df['NiO']==0., 'NiO'] = df['Ni'] * 1.2725 / 10000.
    except KeyError:
        pass
    try:
        df.loc[df['CoO']==0., 'CoO'] = df['Co'] * 1.2715 / 10000.
    except KeyError:
        pass
    try:
        df.loc[df['Cr2O3']==0., 'Cr2O3'] = df['Cr'] * 1.4616 / 10000.
    except KeyError:
        pass
        
    # Calculate H2O value if H2O value is zero
    # parameterizes water by converting Ce to H20
    if (df['H2O'] == 0.).any():
        try:
            df.loc[df['H2O']==0., 'H2O'] = df['Ce'] * Ce_to_H2O / 10000.
        except KeyError:
            message = "Some sample's H2O still zero after parameterization with Ce."
            warnings.warn(message)

    # Add chosen FeIII_totFe value if none are given
    if 'src_FeIII_totFe' not in df.columns:
        df['src_FeIII_totFe'] = src_FeIII_totFe
    df.loc[df['src_FeIII_totFe'] == 0, 'src_FeIII_totFe'] = src_FeIII_totFe
    df.loc[df['FeO_tot']==0., 'FeO_tot'] = df['FeO'] + (df['Fe2O3'] * (74.84 * 2.) / 159.69)
    
    # Calculate FeO and Fe2O3 based on FeIII_totFe    
    # Only applicable if FeO and Fe2O3 not already stated.
    df.loc[(df['FeO']==0.) | (df['Fe2O3']==0.), 'FeO'] = df['FeO_tot'] * (1. - df['src_FeIII_totFe'])
    df.loc[(df['FeO']==0.) | (df['Fe2O3']==0.), 'Fe2O3'] = df['FeO_tot'] * df['src_FeIII_totFe'] * 159.69 / 71.84 / 2.      

    # Calculate CO2 value if desired and CO2 value is zero
    # parameterize CO2 using Equation 8 of Sun & Dasgupta, 2020    
    if param_co2:
        df.loc[df['CO2'] == 0, 'CO2'] = 43.77 - 0.9 * df['SiO2']

    # Calculate total & normalise
    df['Total'] = df[major_columns].sum(axis=1)
    for col in major_columns:
        df[col] = df[col] / df['Total'] * 100.
    df['Total'] = df[major_columns].sum(axis=1)

    # Filter dataset to only include usable samples
    # Remove lines with MgO < 8.5 wt%
    df = df.loc[(df['MgO'] > min_MgO)]

    return df