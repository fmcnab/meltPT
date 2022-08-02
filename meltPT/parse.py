"""
=====
parse
=====

Read data from a csv file.

"""

import warnings

import numpy as np
import pandas as pd

def parse_csv(infile, Ce_to_H2O=0., src_FeIII_totFe=0., min_MgO=0., read_as_primary=False, param_co2=False):
    """
    Parse csv.
    """


    # Read in file
    df = pd.read_csv(infile, delimiter=",")

    # Replace empties and NaNs with zeros
    df = df.replace(r'^\s*$', 0., regex=True)
    df = df.replace(np.nan, 0., regex=True)

    # If these columns do not exist in the input make them and
    # give them zeros for every row.
    check_cols=[
        'Cr2O3', 'MnO', 'H2O', 'CO2', 'FeO', 'Fe2O3', 'src_FeIII_totFe', 'Ce', 
        'P2O5', 'NiO', 'CoO', 'Ni', 'Co', 'Cr']
    df = df.reindex(df.columns.union(check_cols, sort=False), axis=1, fill_value=0)

    # If not given compute some major oxides from trace element concentrations
    # df.loc[df['P2O5'] == 0, 'P2O5'] = df['P'] * 141.942524 / 2. / 10000.
    df.loc[df['NiO'] == 0, 'NiO'] = df['Ni'] * 1.2725 / 10000.
    df.loc[df['CoO'] == 0, 'CoO'] = df['Co'] * 1.2715 / 10000.
    df.loc[df['Cr2O3'] == 0, 'Cr2O3'] = df['Cr'] * 1.4616 / 10000.

    # Calculate H2O value if H2O value is zero
    # parameterizes water by converting Ce to H20
    df.loc[df['H2O'] == 0, 'H2O'] = df['Ce'] * Ce_to_H2O / 10000.

    # Add chosen FeIII_totFe value if none are given
    df.loc[df['src_FeIII_totFe'] == 0, 'src_FeIII_totFe'] = src_FeIII_totFe
    df['FeO_tot'] = df['FeO'] + (df['Fe2O3'] * (74.84 * 2.) / 159.69)
    
    # Calculate FeO and Fe2O3 based on FeIII_totFe    
    # Only applicable if FeO and Fe2O3 not already stated.
    df.loc[(df['FeO']==0.) | (df['Fe2O3']==0.), 'FeO'] = df['FeO_tot'] * (1. - df['src_FeIII_totFe'])
    df.loc[(df['FeO']==0.) | (df['Fe2O3']==0.), 'Fe2O3'] = df['FeO_tot'] * df['src_FeIII_totFe'] * 159.69 / 71.84 / 2.      
    
    if param_co2:
        
        # Calculate CO2 value if CO2 value is zero
        # parameterize CO2 using Equation 8 of Sun & Dasgupta, 2020
        df.loc[df['CO2'] == 0, 'CO2'] = 43.77 - 0.9 * df['SiO2']

    # Calculate total
    major_cols = ['SiO2','Al2O3','FeO','Fe2O3','MgO','CaO','Na2O','K2O','TiO2','MnO','Cr2O3','H2O','CO2']
    df['Total'] = df[major_cols].sum(axis=1)

    # Normalise major elements to sum to 100%
    df['SiO2'] = df['SiO2'] / df['Total'] * 100.
    df['Al2O3'] = df['Al2O3'] / df['Total'] * 100.
    df['FeO'] = df['FeO'] / df['Total'] * 100.
    df['Fe2O3'] = df['Fe2O3'] / df['Total'] * 100.
    df['MgO'] = df['MgO'] / df['Total'] * 100.
    df['CaO'] = df['CaO'] / df['Total'] * 100.
    df['Na2O'] = df['Na2O'] / df['Total'] * 100.
    df['K2O'] = df['K2O'] / df['Total'] * 100.
    df['TiO2'] = df['TiO2'] / df['Total'] * 100.
    df['MnO'] = df['MnO'] / df['Total'] * 100.
    df['Cr2O3'] = df['Cr2O3'] / df['Total'] * 100.
    df['H2O'] = df['H2O'] / df['Total'] * 100.
    df['CO2'] = df['CO2'] / df['Total'] * 100.    
    df['Total'] = df[major_cols].sum(axis=1)

    # Filter dataset to only include usable samples
    # Remove lines with MgO < 8.5 wt%
    df = df.loc[(df['MgO'] > min_MgO)]
    
    # Create warning!
    # Remove lines with no water values
    #df = df.loc[(df['H2O'] > 0.)]

    return df