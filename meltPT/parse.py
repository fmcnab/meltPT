import warnings

import numpy as np
import pandas as pd

def parse_csv(infile, Ce_to_H2O=200., src_FeIII_totFe=0.2, min_SiO2=0., min_MgO=0., read_as_primary=False):

    # Read in file
    df = pd.read_csv(infile, delimiter=",")

    # Replace empties and NaNs with zeros
    # df = df.replace(r'^\s*$', 0., regex=True)
    # df = df.replace(np.nan, 0., regex=True)

    # If these columns do not exist in the input make them and
    # give them zeros for every row.
    check_cols=[
        'Cr2O3', 'MnO', 'H2O', 'FeO', 'Fe2O3', 'src_FeIII_totFe', 'Ce', 
        'P2O5', 'NiO', 'CoO', 'Ni', 'Co']
    df = df.reindex(df.columns.union(check_cols, sort=False), axis=1, fill_value=0)

    # If not given compute some major oxides from trace element concentrations
    # df.loc[df['P2O5'] == 0, 'P2O5'] = df['P'] * 141.942524 / 2. / 10000.
    df.loc[df['NiO'] == 0, 'NiO'] = df['Ni'] * 74.69239999999999 / 10000.
    df.loc[df['CoO'] == 0, 'CoO'] = df['Co'] * 74.932195 / 10000.

    if not read_as_primary:

        # Calculate H2O value if H2O value is zero
        # parameterizes water by converting Ce to H20
        df.loc[df['H2O'] == 0, 'H2O'] = df['Ce'] * Ce_to_H2O / 10000.
        
        # Add chosen FeIII_totFe value if none are given
        df.loc[df['src_FeIII_totFe'] == 0, 'src_FeIII_totFe'] = src_FeIII_totFe

        # Calculate FeO and Fe2O3 based on FeIII_totFe
        df['FeO_tot'] = df['FeO'] + (df['Fe2O3'] * (74.84 * 2.) / 159.69)
        df['FeO'] = df['FeO_tot'] * (1. - df['src_FeIII_totFe'])
        df['Fe2O3'] = df['FeO_tot'] * df['src_FeIII_totFe'] * 159.69 / 71.84 / 2.

    # Calculate total
    major_cols = ['SiO2','Al2O3','FeO','Fe2O3','MgO','CaO','Na2O','K2O','TiO2','MnO','Cr2O3','H2O']
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
    df['Total'] = df[major_cols].sum(axis=1)

    # Filter dataset to only include usable samples
    # Remove lines with SiO2 < 40 wt%
    df = df.loc[(df['SiO2'] > min_SiO2)]
    # Remove lines with MgO < 8.5 wt%
    df = df.loc[(df['MgO'] > min_MgO)]
    # Remove lines with no water values
    df = df.loc[(df['H2O'] > 0.)]

    # df_wt_hydrous = df[major_cols]
    # to_drop = np.hstack((major_cols, "Sample", "Total"))
    # df_other_info = df.drop(to_drop, axis=1)
    #
    # samples = []
    # for i in range(len(df)):
    #     sample = Sample(df["Sample"][i], df_wt_hydrous.to_dict("records")[i], df_other_info.to_dict("records")[i])
    #     samples.append(sample)

    return df