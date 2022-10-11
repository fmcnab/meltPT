"""
======================
backtrack_compositions
======================

Backtrack sample compositions to 'primary' compositions.

"""

import warnings

import numpy as np
import pandas as pd


MAJOR_OXIDES = [
    'SiO2','Al2O3','FeO','Fe2O3','MgO','CaO','Na2O','K2O','TiO2','MnO',
    'Cr2O3', 'P2O5', 'NiO', 'CoO', 'H2O', 'CO2']

def normalise(in_comp):
    """
    Normalise dictionary so that values sum to 100%.
    
    Parameters
    ---------
    in_comp : dict
        The values to be normalised.
        
    Returns
    -------
    out_comp : dict
        The normalised values.
    """
    out_comp = {}
    total = sum(in_comp.values())
    for phase in in_comp:
        out_comp[phase] = in_comp[phase] / total * 100.
    return out_comp

def oxide_wt_to_cation_mole(in_comp):
    """
    Convert oxide concentrations in weight percent to cation concentrations.
    
    Parameters
    ----------
    in_comp : dict
        The oxide concentrations to be converted.
        
    Returns
    -------
    out_comp : dict
        The cation concentrations.
    """

    weights = {
        # Combined atomic weights, e.g., SiO2 = Si + O + O = 60.08 (g/mol).
        # Atomic weights given per cation, e.g., Al2O3 is divided by 2.
        'SiO2': 60.08,
        'Al2O3': 101.96 / 2.,
        'FeO': 71.84,
        'Fe2O3': 159.69 / 2.,
        'MgO': 40.3,
        'CaO': 56.08,
        'Na2O': 61.98 / 2.,
        'K2O': 94.2 / 2.,
        'TiO2': 79.86,
        'MnO': 70.94,
        'Cr2O3': 151.99 / 2.,
        'P2O5': 141.942524 / 2.,
        'NiO': 74.69239999999999,
        'CoO': 74.932195,
        'H2O': 18.014680000000002 / 2.,
        'CO2': 44.009
        }
    keys = {
        'SiO2': 'Si',
        'Al2O3': 'Al',
        'FeO': 'FeII',
        'Fe2O3': 'FeIII',
        'MgO': 'Mg',
        'CaO': 'Ca',
        'Na2O': 'Na',
        'K2O': 'K',
        'TiO2': 'Ti',
        'MnO': 'Mn',
        'Cr2O3': 'Cr',
        'P2O5': 'P',
        'NiO': 'Ni',
        'CoO': 'Co',
        'H2O': 'H',
        'CO2': 'C'
    }
    out_comp = {}
    for phase in in_comp:
        out_comp[keys[phase]] = in_comp[phase] / weights[phase]
    return normalise(out_comp)

def oxide_wt_to_oxide_mole(in_comp):
    """
    Convert oxide concentrations in weight percent to oxide concentrations
    in mole percent.
    
    Parameters
    ----------
    in_comp : dict
        The oxide concentrations in weight percent to be converted.
        
    Returns
    -------
    out_comp : dict
        The oxide concentrations in mole percent.
    """

    weights = {
        # Combined atomic weights. E.g., Si, O and O = 60.08 (g/mol).
        'SiO2': 60.08,
        'Al2O3': 101.96,
        'FeO': 71.84,
        'Fe2O3': 159.69,
        'MgO': 40.3,
        'CaO': 56.08,
        'Na2O': 61.98,
        'K2O': 94.2,
        'TiO2': 79.86,
        'MnO': 70.94,
        'Cr2O3': 151.99,
        'P2O5': 141.942524,
        'NiO': 74.69239999999999,
        'CoO': 74.932195,
        'H2O': 18.014680000000002,
        'CO2': 44.009
        }
    out_comp = {}
    for phase in in_comp:
        out_comp[phase] = in_comp[phase] / weights[phase]
    return normalise(out_comp)

def oxide_mole_to_mole_species(in_comp):
    """
    Convert oxide concentrations in mole percent to mole species concentrations.
    
    Parameters
    ----------
    in_comp : dict
        The oxide concentrations to be converted.
        
    Returns
    -------
    out_comp : dict
        The mole species concentrations.
    """

    out_comp = {
        'Si4O8':
            ((in_comp['SiO2'] -
            0.5*(in_comp['FeO']+in_comp['MgO']+in_comp['CaO']+in_comp['MnO']) -
            in_comp['Na2O'] -
            in_comp['K2O']) * 0.25),
        'Al16/3O8': (in_comp['Al2O3'] - in_comp['Na2O']) * (3./8.),
        'Fe4Si2O8': in_comp['FeO'] * 0.25,
        'Fe16/3O8': in_comp['Fe2O3'] * (3./8.),
        'Mg4Si2O8': in_comp['MgO'] * 0.25,
        'Ca4Si2O8': in_comp['CaO'] * 0.25,
        'Na2Al2Si2O8': in_comp['Na2O'],
        'K2Al2Si2O8': in_comp['K2O'],
        'Ti4O8': in_comp['TiO2'] * 0.25,
        'Mn4Si2O8': in_comp['MnO'] * 0.25,
        'Cr16/3O8': in_comp['Cr2O3']  * (3./8.),
        'H16O8': in_comp['H2O'] * 0.125
    }
    return normalise(out_comp)

def compute_partition_coefficient(in_comp):
    """
    Compute the partition coefficient.
    
    Expression from Tamura et al. (2000, J. Pet), via Lee et al. (2009)
    spreadsheet. Spreadsheet has more significant figures than paper...
    """
    Kd = 0.25324 + 0.0033663*(in_comp['Mg'] + 0.33*in_comp['FeII'])
    return Kd

def compute_forsterite_number(in_comp, Kd=None):
    """
    Compute forsterite number from oxide weight compositions.
    
    Parameters
    ----------
    in_comp : dict
        The oxide weight compositions.
    Kd : float or NoneType, optional
        Partition coefficient to be used.
        If None, calculated from the sample composition provided.
        
    Returns
    -------
    Fo : float
        The calculated forsterite number.
    """
    cation = oxide_wt_to_cation_mole(in_comp)
    if not Kd:
        Kd = compute_partition_coefficient(cation)
    Fo = 1. / (1. + (Kd * (cation['FeII'] / cation['Mg'])))
    return Fo

def add_olivine(in_comp, Kd=None, dm=0.0005):
    """
    Add olivine in equilibrium with given melt composition.
    
    Parameters
    ----------
    in_comp : dict
        The initial composition.
        Should be hydrous oxide weight concentrations.
    Kd : float or NoneType, optional
        Partition coefficient to be used.
        If None, calculated from the sample composition provided.
    dm : float, optional
        The fraction of olivine to be added.
        
    Returns
    -------
    out_comp : dict
        The updated concentrations.
    """

    # Compute melt forsterite number
    Fo = compute_forsterite_number(in_comp, Kd=Kd)

    # Compute olivine composition in equilibrium with melt
    oxide_wt_olivine = {
        'FeO': 2. * (1 - Fo) * 71.85,
        'MgO': 2. * Fo * 40.3,
        'SiO2': 60.08
        }
    oxide_wt_olivine = normalise(oxide_wt_olivine)

    # loop over phases adding olivine
    out_comp = {}
    for phase in in_comp:
        if phase == 'SiO2' or phase == 'FeO' or phase == 'MgO':
            out_comp[phase] = (in_comp[phase] + dm*oxide_wt_olivine[phase]) / (1. + dm)
        else:
            out_comp[phase] = in_comp[phase] / (1. + dm)

    return normalise(out_comp)

def fill_dict_with_nans(in_dict):
    """
    Fill a dictionary with numpy.nan values.
    
    Parameters
    ----------
    in_dict : dict
        The dictionary to be filled with nans.
        
    Returns
    out_dict : dict
        The dictionary with all values replaces with nan.
    """
    out_dict = {}
    for key in in_dict:
        out_dict[key] = np.nan
    return out_dict

def backtrack_sample_composition(
    df, target_Fo=0.9, Kd=None, dm=0.0005, verbose=False, 
    max_olivine_addition=0.3, return_all=False):
    """
    Backtrack composition to desired mantle forsterite number.
    
    Iteratively adds olivine in equilibrium with melt until desired composition
    is reached.
    
    Parameters
    ----------
    df : pandas dataframe
        Dataframe containing the initial composition to be backtracked.
        Should contain only one row. To use with a multi-row dataframe use
        df.apply().
    target_Fo : float, optional
        The forsterite number of the mantle source.
        We add iteratively add olivine to the sample composition until this
        value is reached.
    Kd : float or NoneType, optional
        Partition coefficient to be used.
        If None, calculated from the sample composition provided.
    dm : float, optional
        The fraction of olivine to be added at each iteration.
    verbose : bool, optional
        If True, will print messages with information at each iteration.
    max_olivine_addition : float, optional
        Maximum fraction of olivine to add before backtracking is abandoned.
    return_all : bool
        Return intermediate backtracking compositions.
    
    Returns
    -------
    primary_oxide : df
        The backtracked compositions.
    composition_through_addition : list
        If return_all is True, returns list of dictionaries containing
        intermediate backtracking compositions.
    """

    # Get major oxides from data frame
    oxide_wt_hydrous = {}
    for ox in MAJOR_OXIDES:
        oxide_wt_hydrous[ox] = df[ox]
    oxide_wt_hydrous = normalise(oxide_wt_hydrous)
    
    # Testing Lee et al. spreadsheet
    # Seems to not update Kd during iteration
    # if not Kd:
    #     cation = oxide_wt_to_cation_mole(oxide_wt_hydrous)
    #     Kd = compute_partition_coefficient(cation)

    # Check Fo is below mantle Fo
    Fo = compute_forsterite_number(oxide_wt_hydrous)
    if target_Fo-Fo < 0.001:
        oxide_wt_hydrous = fill_dict_with_nans(oxide_wt_hydrous)
        dm_tot = np.nan
        message = df.Sample + ": backtracking failed! Starting Fo above mantle Fo."
        warnings.warn(message)
    # Otherwise add olvine until primary Fo is reached
    else:

        if verbose:

            print("Backtracking sample %s to primary composition:" % df.Sample)

        dm_tot = 0.
        composition_through_addition = []
        # while abs(target_Fo - Fo) > 1.e-15:
        while target_Fo - Fo > 0.0002:
            
            oxide_wt_hydrous = add_olivine(oxide_wt_hydrous, Kd=Kd, dm=dm)
            composition_through_addition.append(oxide_wt_hydrous)
            dm_tot += dm
            Fo = compute_forsterite_number(oxide_wt_hydrous, Kd=Kd)
            if verbose:
                cation = oxide_wt_to_cation_mole(oxide_wt_hydrous)
                if not Kd:
                    Kd_i = compute_partition_coefficient(cation)
                else:
                    Kd_i = Kd
                print(
                    "    - %.2f%% olivine added, melt Fo = %.4f, Kd = %.4f." %
                    (dm_tot/(1.+dm_tot)*100., Fo, Kd_i)
                    )
            if dm_tot/(1. + dm_tot) > max_olivine_addition:
                oxide_wt_hydrous = fill_dict_with_nans(oxide_wt_hydrous)
                message = (
                    df.Sample + ": backtracking failed! Olivine addition exceeding %d%%" 
                    % (max_olivine_addition*100.))
                warnings.warn(message)
                break


    # # Anhydrous concentrations: remove water and renormalise
    # oxide_wt_anhydrous = {phase:oxide_wt_hydrous[phase] for phase in oxide_wt_hydrous.keys() if phase != "H2O"}
    # oxide_wt_anhydrous = normalise(oxide_wt_anhydrous)
    # 
    # # Compute other concentrations
    # oxide_mole_hydrous = oxide_wt_to_oxide_mole(oxide_wt_hydrous)
    # mole_species_hydrous = oxide_mole_to_mole_species(oxide_mole_hydrous)
    # mole_species_anhydrous = {phase:mole_species_hydrous[phase] for phase in mole_species_hydrous.keys() if phase != "H16O8"}
    # mole_species_anhydrous = normalise(mole_species_anhydrous)
    # 
    # # Package up
    primary_oxide = {}
    for phase in oxide_wt_hydrous:
        primary_oxide[phase + "_primary_wt"] = oxide_wt_hydrous[phase]
    # for phase in oxide_wt_anhydrous:
    #     primary_oxide[phase + "_primary_wt_dry"] = oxide_wt_anhydrous[phase]
    # for phase in oxide_mole_hydrous:
    #     primary_oxide[phase + "_primary_mol"] = oxide_mole_hydrous[phase]
    # for phase in mole_species_hydrous:
    #     primary_oxide[phase] = mole_species_hydrous[phase]
    # for phase in mole_species_anhydrous:
    #     primary_oxide[phase + "_dry"] = mole_species_anhydrous[phase]
        
    # Add amount olivine added
    primary_oxide['ol_added'] = dm_tot / (1. + dm_tot)

    if not return_all:
        return primary_oxide
    else:
        return primary_oxide, composition_through_addition
