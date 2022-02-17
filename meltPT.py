import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from copy import deepcopy
import sys
import pyMelt as m
from scipy.optimize import minimize, minimize_scalar
import shapely.geometry as shp

def parse_csv(infile, Ce_to_H2O=200., src_FeIII_totFe=0.2, min_SiO2=0., min_MgO=0.):

    # Read in file
    df = pd.read_csv(infile, delimiter=",")

    # Replace empties and NaNs with zeros
    df = df.replace(r'^\s*$', 0., regex=True)
    df = df.replace(np.nan, 0., regex=True)

    # If these columns do not exist in the input make them and
    # give them zeros for every row.
    check_cols=['Cr2O3', 'MnO', 'H2O', 'FeO', 'Fe2O3', 'src_FeIII_totFe', 'Ce']
    df = df.reindex(df.columns.union(check_cols, sort=False), axis=1, fill_value=0)

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

MAJOR_OXIDES = ['SiO2','Al2O3','FeO','Fe2O3','MgO','CaO','Na2O','K2O','TiO2','MnO','Cr2O3','H2O']

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
        'H2O': 18.014680000000002 / 2.
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
        'H2O': 'H'
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
        'H2O': 18.014680000000002
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
        'Cr16/3O8': in_comp['Cr2O3']  * (3./8.)
        # 'H16O8': in_comp['H2O'] * 0.125
    }
    return normalise(out_comp)

def compute_partition_coefficient(in_comp):
    """
    Compute the partition coefficient.
    """
    Kd = 0.25324 + 0.003363*(in_comp['Mg'] + 0.33*in_comp['FeII'])
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
    max_olivine_addition=0.3):
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
    
    Returns
    -------
    primary_oxide : df
        The backtracked compositions.
    """

    # Get major oxides from data frame
    oxide_wt_hydrous = {}
    for ox in MAJOR_OXIDES:
        oxide_wt_hydrous[ox] = df[ox]

    # Check Fo is below mantle Fo
    Fo = compute_forsterite_number(oxide_wt_hydrous)
    if target_Fo-Fo < 0.005:
        oxide_wt_hydrous = fill_dict_with_nans(oxide_wt_hydrous)
        dm_tot = np.nan
        print(df.Sample + ": backtracking failed! Starting Fo above mantle Fo.")
    # Otherwise add olvine until primary Fo is reached
    else:

        if verbose:
            print("Backtracking sample %s to primary composition:" % df.Sample)

        dm_tot = 0.
        # while abs(target_Fo - Fo) > 1.e-15:
        while abs(target_Fo - Fo) > dm:
            oxide_wt_hydrous = add_olivine(oxide_wt_hydrous, Kd=Kd, dm=dm)
            dm_tot += dm
            Fo = compute_forsterite_number(oxide_wt_hydrous, Kd=Kd)
            if verbose:
                print(
                    "    - iteration %d: %.2f%% olivine added, melt Fo = %.4f." %
                    (i, dm_tot/(1.+dm_tot)*100., Fo)
                    )
            if dm_tot/(1. + dm_tot) > max_olivine_addition:
                oxide_wt_hydrous = fill_dict_with_nans(oxide_wt_hydrous)
                print(
                    df.Sample + ": backtracking failed! Olivine addition exceeding %d%%" 
                    % (max_olivine_addition*100.))
                break


    # Anhydrous concentrations: remove water and renormalise
    oxide_wt_anhydrous = {phase:oxide_wt_hydrous[phase] for phase in oxide_wt_hydrous.keys() if phase != "H2O"}
    oxide_wt_anhydrous = normalise(oxide_wt_anhydrous)

    # Compute other concentrations
    oxide_mole_hydrous = oxide_wt_to_oxide_mole(oxide_wt_hydrous)
    mole_species = oxide_mole_to_mole_species(oxide_mole_hydrous)

    # Package up
    primary_oxide = {}
    for phase in oxide_wt_hydrous:
        primary_oxide[phase + "_primary_wt"] = oxide_wt_hydrous[phase]
    for phase in oxide_wt_anhydrous:
        primary_oxide[phase + "_primary_wt_dry"] = oxide_wt_anhydrous[phase]
    for phase in oxide_mole_hydrous:
        primary_oxide[phase + "_primary_mol"] = oxide_mole_hydrous[phase]
    for phase in mole_species:
        primary_oxide[phase] = mole_species[phase]

    # Add amount olivine added
    primary_oxide['ol_added'] = dm_tot / (1. + dm_tot)

    return primary_oxide

def compute_water_correction(df):
    """
    Calculate temperature correction from water content.

    Equation (1), Plank and Forsyth (2016, G-cubed).
    
    Parameters
    ----------
    df : pandas dataframe
        The sample compositions.
        
    Returns
    -------
    correction : pandas dataframe
        The temperature corrections.
    """
    
    correction = 40.4*df['H2O_primary_wt'] - 2.97*df['H2O_primary_wt']**2. + 0.0761*df['H2O_primary_wt']**3.
    return correction

def compute_CO2_correction(df):
    """
    Calculate temperature correction for samples with equlibration pressures
    above 2 GPa (assumed to be influenced by presence of CO2).

    Equation (3), Plank and Forsyth (2016, G-cubed).
    
    Parameters
    ----------
    df : pandas dataframe
        The sample compositions.
        
    Returns
    -------
    correction : pandas dataframe
        The temperature corrections.
    """
    correction = (df['SiO2_primary_wt_dry'] - 50.3) / -0.12804
    return correction

def compute_temperature(df):
    """
    Calculate equlibration temperature.

    Equation (1), Plank and Forsyth (2016, G-cubed).
    
    Parameters
    ----------
    df : pandas dataframe
        The sample compositions.
        
    Returns
    -------
    temperature : pandas dataframe
        The calculated temperature(s).
    """
    temperature = (
        1264.5 +
        7.85*df['Mg4Si2O8'] +
        8545./df['Si4O8'] -
        5.96*df['Al16/3O8']
        )
    return temperature

def compute_pressure(df, T):
    """
    Calculate equlibration pressure.

    Equation (2), Plank and Forsyth (2016, G-cubed).

    Parameters
    ----------
    df : pandas dataframe
        The sample compositions.
        
    Returns
    -------
    pressure : pandas dataframe
        The calculated pressure(s).
    """
    pressure = (
        np.log(df['Si4O8']) -
        4.045 +
        0.0114*df['Fe4Si2O8'] +
        0.00052*df['Ca4Si2O8']**2. +
        0.0024*df['Mg4Si2O8']) / (-336.3/T - 0.0007*np.sqrt(T)
        )
    return pressure

def compute_sample_pressure_temperature(df):
    """
    Compute equilibration pressure and temperature for a given sample.

    Parameters
    ----------
    df : pandas dataframe
        Dataframe containing the sample primary composition.
        Should contain only one row. To use with a multi-row dataframe use
        df.apply().
        
    Returns
    -------
    out : dict
        The equilibration pressure and temperature result.
        Temperature is in degrees celcius.
        Pressure is in GPa.
    """
    water_correction = compute_water_correction(df)
    T = compute_temperature(df) - water_correction
    P = compute_pressure(df, T)
    if P > 2.:
        T -= compute_CO2_correction(df)
        P = compute_pressure(df, T)
    out = {'P': P, 'T': T - 273.15}
    return out

# ---- Melt Path Fitting

def melt_fraction_to_pressure_temperature(F, path):
    """
    Find pressure and temperature corresponding to specified melt fraction
    for a specified melting path.
    
    Parameters
    ----------
    F : float
        The melt fraction.
    path : instance of pyMelt.meltingcolumn_classes.meltingColumn
        The melt path.
        Can be instance of any class containing arrays of melt fraction (F),
        pressure (P) and temperature (T).
        
    Returns
    -------
    P : float
        Pressure corresponding to specified melt fraction.
    T : float 
        Temperature corresponding to specified melt fraction.
    """
    P = np.interp(F, path.F, path.P)
    T = np.interp(F, path.F, path.T)
    return P, T
    
def compute_sample_melt_fraction_misfit(F, df, path, P_err=0.24, T_err=39.):
    """
    Compute misfit between sample and point along given melting path,
    specified by melt path.
    
    Parameters
    ----------
    F : float
        The melt fraction with which to compare the pressure/temperature
        point.
    df : pandas dataframe
        Dataframe containing sample pressure and temperature estimates.
    path : instance of pyMelt.meltingcolumn_classes.meltingColumn
        The melt path.
        Can be instance of any class containing arrays of melt fraction (F),
        pressure (P) and temperature (T).
    P_err : float, optional
        Uncertainty in pressure observation(s).
    T_err : float, optional
        Uncertainty in temperature observation(s).

    Returns
    -------
    misfit : float
        Misfit between observed pressure and temperature and point along
        melt path. Defined as distance in normalised pressure-temperature
        space. Pressure and temperature are normalised by their respective
        uncertainties.
    """
    model_P, model_T = melt_fraction_to_pressure_temperature(F, path)
    misfit = np.sqrt( 
        ((df['P']-model_P)/P_err)**2. + ((df['T']-model_T)/T_err)**2.
        )
    return misfit
    
def find_sample_melt_fraction(df, path):
    """
    Find best-fitting melt fraction for sample along given melting path.
    
    Uses scipy's minimize_scalar function to find melt fraction that minimizes 
    misfit between estimated pressure and temperature and the melt path.
    Options used are:
        method: "bounded"
        bounds:  min --> 0
                 max --> Maximum melt fraction in "path"
        bracket: min --> 0
                 max --> Maximum melt fraction in "path"
        
    Parameters
    ----------
    df : pandas dataframe
        Dataframe containing sample pressure and temperature estimates.
        Should contain only one row. To use with a multi-row dataframe use
        df.apply().
    path : instance of pyMelt.meltingcolumn_classes.meltingColumn
        The melting path.
        Can be instance of any class containing arrays of melt fraction (F),
        pressure (P) and temperature (T).
        
    Returns
    -------
    out : dict
        Various properties of the result:
            F : float 
                The best-fitting melt fraction.
            P : float
                Pressure on path corresponding to best-fitting melt fraction.
            T : float
                Temperature on path corresponding to best-fitting melt fraction.
            misfit : float
                Distance between sample and closest point on melting path.
        If sample pressure or temperature are nan, all returned values are nan.
    """
    if np.isnan(df['P']) or np.isnan(df['T']):
        out = {'F': np.nan, 'P': np.nan, 'T': np.nan, 'misfit': np.nan}
    else:
        fit = minimize_scalar(
            compute_sample_melt_fraction_misfit, 
            bounds=(0.,max(path.F)), 
            bracket=(0.,max(path.F)), 
            args=(df, path), 
            method="bounded")
        P, T = melt_fraction_to_pressure_temperature(fit.x, path)
        out = {'F': fit.x, 'P': P, 'T': T, 'misfit': fit.fun}
    return out

def compute_sample_potential_temperature_misfit(Tp, df, mantle):
    """
    Compute a melting path for a given potential temperature then find misfit
    between it and a sample pressure-temperature estimate.
    
    First uses mantle.adiabaticMelt() to compute melt path for specified
    potential temperature. Then finds closest melt fraction 
    
    Parameters
    ----------
    Tp : float
        The potential temperature to be used.
    df : pandas dataframe
        Dataframe containing sample pressure and temperature estimates.
        Should contain only one row. To use with a multi-row dataframe use
        df.apply().
    mantle : instance of pyMelt.mantle_class.mantle
        The mantle object to be used to calculate the melting path.
    
    Returns
    -------
    misfit : float
        Distance between sample pressure-temperature estimate and its nearest
        point on the calculated melting path.
    """
    path = mantle.adiabaticMelt(
        Tp, 
        Pstart=max(mantle.solidusIntersection(Tp))+0.01,
        dP=-0.01)
    misfit = find_sample_melt_fraction(df, path)['misfit']
    return misfit

def find_sample_potential_temperature(df, mantle):
    """
    Find best-fitting potential temperature for a sample pressure-temperature
    estimate.
    
    Uses scipy's minimize_scalar function to find potential temperature that
    minimizes misfit between estimated pressure and temperature and the
    corresponding  melt path. Options used are:
        method: "bounded"
        bounds:  min --> 
                 max --> 1600 oC
        bracket: min --> intersection of solidus with surface
                 max --> 1600 oC
                 
    Parameters
    ----------
    df : pandas dataframe
        Dataframe containing sample pressure and temperature estimates.
        Should contain only one row. To use with a multi-row dataframe use
        df.apply().
    mantle : instance of pyMelt.mantle_class.mantle
        The mantle object to be used to calculate the melting paths.
        
    Returns
    -------
    out : dict
        Various properties of the result:
            F : float 
                The melt fraction along the best-fitting melting path.
            P : float
                The pressure along the best-fitting melting path.
            T : float
                The temperature along the best-fitting melting path.
            misfit : float
                Distance between sample and closest point on melting path.
            Tp : float
                The best-fitting potential temperature.
            path : instance of pyMelt.meltingcolumn_classes.meltingColumn
                The best-fitting melting path.
        If sample pressure or temperature are nan, all returned values are nan.
    """
    if np.isnan(df['P']) or np.isnan(df['T']):
        out = {
            'F': np.nan, 'P': np.nan, 'T': np.nan, 
            'misfit': np.nan, 'Tp': np.nan, 'path': np.nan}
    else:    
        Tp_fit = minimize_scalar(
            compute_sample_potential_temperature_misfit, 
            bracket=(min([lith.TSolidus(0.) for lith in mantle.lithologies]),1600.), 
            bounds=(min([lith.TSolidus(0.) for lith in mantle.lithologies]),1600.), 
            args=(df,mantle), 
            method="bounded")
        path = mantle.adiabaticMelt(
            Tp_fit.x, 
            Pstart=max(mantle.solidusIntersection(Tp_fit.x))+0.01, 
            dP=-0.01)
        out = find_sample_melt_fraction(df, path)
        out['Tp'] = Tp_fit.x 
        out['path'] = path
    return out
        
def compute_suite_potential_temperature_misfit(Tp, df, mantle):
    """
    Compute a melt path for a specified potential temperature then calculate
    the misfit between it and one or more pressure-temperature estimates.
    
    Parameters
    ----------
    Tp : float
        The potential temperature to be used.
    df : pandas dataframe
        Dataframe containing sample pressure and temperature estimates.
        Can contain any number of samples.
    mantle : instance of pyMelt.mantle_class.mantle
        The mantle object to be used to calculate the melting paths.
    
    Returns
    -------
    misfit : float
        The average distance between sample pressure-temperature estimates
        and their nearest points on the calculated melting path.
    """
    path = mantle.adiabaticMelt(
        Tp,
        Pstart=max(mantle.solidusIntersection(Tp))+0.01, 
        dP=-0.01)
    melt_fraction_fits = df.apply(
        find_sample_melt_fraction, 
        axis=1, 
        result_type="expand", 
        args=(path,))
    misfit = np.nanmean(melt_fraction_fits['misfit'])
    return misfit

def find_bounding_potential_temperature(points, starting_temperature, mantle, lower=False, threshold=(2./3.)):
    """
    Find either upper or lower bound on best-fitting potential temperature for
    suite of pressure-temperature estimates.
    
    Works by computing melting paths progressively further away from the best-
    fitting path, until a threshold number of points lie between the two paths.
    
    Parameters
    ----------
    points : list of shapely.geometry.point.Point objects
        The points to be bounded.
    starting_temperature : float
        The best-fitting potential temperature for the suite.
    mantle : instance of pyMelt.mantle_class.mantle
        The mantle object to be used to calculate the melting paths.
    lower : bool, optional
        Specify whether an upper or lower bound is to be found.
    threshold : float
        The threshold fraction of points to be lie between the best-fitting
        and bounding temperature melting paths.
    
    Returns
    -------
    bounding_temperature : float
        The estimated bounding temperature.
    bounding_path : instance of pyMelt.meltingcolumn_classes.meltingColumn
        The melting path corresonding to the bounding potential temperature.
    """
    
    # First check whether an upper or lower bound is to be found and define
    # initial guess and the adjustment direction accordingly.
    # If searching for an upper bound, initial guess given by rounding up
    # best-fitting temperature to nearest whole degree and increment is
    # positive. Otherwise, we round down and the increment is negative.
    if not lower:
        bounding_temperature = np.ceil(starting_temperature)
        adjustment = 1.
    else:
        bounding_temperature = np.floor(starting_temperature)
        adjustment = -1.
    
    # Initialise array specifying which points are inside (1) vs. outside (0)
    # the bounding melting paths.
    inside = np.zeros(len(points))

    # Compute melting path corresponding to best-fitting potential temperature
    # Use largest of solidus pressure or maximum pressure in points as the
    # starting pressure.
    max_P = max([p.coords.xy[1][0] for p in points])
    main_path = mantle.adiabaticMelt(
        starting_temperature,
        Pstart=max(max(mantle.solidusIntersection(starting_temperature))+0.01, max_P),
        dP=-0.01)

    # Start incrementally expanding bounds.
    while True:
        
        # Compute melting bath corresponding to bounding potential temperature.
        bounding_path = mantle.adiabaticMelt(
            bounding_temperature,
            Pstart=max(max(mantle.solidusIntersection(bounding_temperature))+0.01, max_P), 
            dP=-0.01)
        
        # Create bounding polygon.
        bounds = np.vstack(( 
            np.column_stack(( main_path.T, main_path.P )),
            np.column_stack(( bounding_path.T[::-1], bounding_path.P[::-1] ))
            ))
        poly = shp.polygon.Polygon(bounds)
        
        # Check which points lie inside polygon.
        for i,p in enumerate(points):
            if poly.contains(p):
                inside[i] = 1.
        
        # Stop if we have reached threshold, otherwise increment and continue.
        if sum(inside) / len(points) > threshold:
            break
        else:
            bounding_temperature += adjustment
            
    return bounding_temperature, bounding_path


def combine(df):
    return {'fit': df.to_numpy().all()}


# ---- Suite class

class Suite:
    """
    Store and process compositions of basaltic rocks.
    
    Includes methods to find primary compositions (i.e. correcting for the
    crystallisation of olivine; Lee et al., 2009, EPSL), compute equilibration
    pressures and temperatures (Plank & Forsyth, 2016, G-cubed), and fit
    melting paths to those pressure-temperature estimates.
    
    Parameters
    ----------
    input_csv : str
        Path to a csv containing data to be read.
    src_FeIII_totFe : float
        Ratio of Fe3+ to total Fe in the mantle source.
    min_SiO2 : float, optional
        Minimum amount of SiO2 in sample to be accepted.
    min_MgO : float, optional
        Minimum amound of MgO in sample to be accepted.
        
    Properties
    ----------
    data : pandas dataframe
        The raw data from the provided csv.
    primary : pandas dataframe or NoneType
        If backtrack_compositions has been run, will contain the estimated
        primary compositions in various forms.
    PT : pandas dataframe or NoneType
        If compute_pressure_temperature has been run, will contain the
        esimated equilibration pressures and temperatures.
    PT_to_fit : pandas dataframe or NoneType
        If any melt-path fitting has been run, will contain pressure and
        temperature estimates of samples that have been selected for fitting.
    individual_melt_fractions : pandas dataframe or NoneType
        If find_individual_melt_fractions, will contain results of fitting
        suite of pressure-temperature estimates to a specified melt path.
    individual_potential_temperatures : pandas dataframe or NoneType
        If find_individual_potential_temperatures has been run, will contain
        results of fitting melting paths to each individual pressure-
        temperature estimate.
    suite_melt_fractions : pandas dataframe or NoneType
        If find_suite_potential_temperature has been run, will contain closest
        points on best-fitting melting path for each pressure-temperature
        estimate.
    potential_temperature : float or NoneType
        If find_suite_potential_temperature has been run, will be the best-
        fitting potential temperature for the suite of pressure-temperature
        estimates.
    upper_potential_temperature : float or NoneType
        If find_suite_potential_temperature has been run, with find_bounds,
        will be the upper bound on potential temperature for the suite of
        pressure-temperature estimates.
    lower_potential_temperature : float or NoneType
        If find_suite_potential_temperature has been run, with find_bounds,
        will be the lower bound on potential temperature for the suite of
        pressure-temperature estimates.
    path : instance of pyMelt.meltingcolumn_classes.meltingColumn
        If find_suite_potential_temperature has been run, will be the best-
        fitting melting path.
    upper_path : instance of pyMelt.meltingcolumn_classes.meltingColumn
        If find_suite_potential_temperature has been run, with find_bounds,
        will be the upper-bounding melting path.
    lower_path : instance of pyMelt.meltingcolumn_classes.meltingColumn
        If find_suite_potential_temperature has been run, with find_bounds,
        will be the lower-bounding melting path.
    
    Methods
    -------
    backtrack_compositions :
        Backtrack compositions for entire dataframe.
    compute_pressure_temperature :
        Compute equilibration pressures and temperatures for entire suite.
    check_samples_for_fitting :
        Determine whether a sample should be fit to or not.
    find_individual_melt_fractions :
        Find best-fitting melt fractions for each sample relative to given melt
        path.
    find_individual_potential_temperatures :
        Find best-fitting potential temperatures and corresponding melt
        fractions for each sample.
    find_suite_potential_temperature :
        Find best-fitting potential temperature for entire suite.
    write_to_csv :
        Write results to csv.
    """

    def __init__(self, input_csv, src_FeIII_totFe=0.19, min_SiO2=0., min_MgO=0.):
        self.data = parse_csv(
            input_csv,
            src_FeIII_totFe=src_FeIII_totFe,
            min_SiO2=min_SiO2,
            min_MgO=min_MgO)
        self.primary = None
        self.PT = None
        self.PT_to_fit = None
        self.individual_melt_fractions = None
        self.individual_potential_temperatures = None
        self.suite_melt_fractions = None
        self.potential_temperature = None
        self.upper_potential_temperature = None
        self.lower_potential_temperature = None
        self.path = None
        self.upper_path = None
        self.lower_path = None

    def backtrack_compositions(self, target_Fo=0.9, Kd=False, dm=0.0005, verbose=False):
        """
        Backtrack compositions for entire dataframe.
        """
        self.primary = self.data.apply(
            backtrack_sample_composition,
            axis=1,
            result_type="expand",
            args=(target_Fo,Kd,dm,verbose)
            )

    def compute_pressure_temperature(self):
        """
        Compute equilibration pressures and temperatures for entire suite.
        """
        self.PT = self.primary.apply(compute_sample_pressure_temperature, axis=1, result_type="expand")
        
    def check_samples_for_fitting(self, mantle, filters=(None,), args=((None,))):
        """
        Determine whether sample should be fit to or not.
        
        Checks if within error of the solidus for given mantle composition.

        Also applies any other filters provided. Filters should be in form of
        function that returns True or False.
        """
                
        def above_solidus(df, mantle):
            return {'fit': max([l.TSolidus(df['P']) for l in mantle.lithologies]) - df['T'] < 39.}

        to_fit = self.PT.apply(above_solidus, axis=1, result_type="expand", args=(mantle,))
        if filters[0]:
            for f,a in zip(filters,args):
                to_fit = pd.concat([to_fit, self.PT.apply(f, axis=1, args=a)], axis=1)
            to_fit = to_fit.apply(combine, axis=1, result_type="expand")

        self.PT_to_fit = self.PT.copy()
        for i,fit in enumerate(to_fit.to_numpy()):
            if not fit:
                self.PT_to_fit.iloc[i]['P'] = np.nan
                self.PT_to_fit.iloc[i]['T'] = np.nan
    
    def find_individual_melt_fractions(self, mantle, path, filters=(None,), filter_args=(None,)):
        """
        Find best-fitting melt fractions for each sample relative to given melt
        path. 
        """
        self.check_samples_for_fitting(mantle, filters, filter_args)
        self.individual_melt_fractions = self.PT_to_fit.apply(
            find_sample_melt_fraction, 
            axis=1, 
            result_type="expand", 
            args=(path,))

    def find_individual_potential_temperatures(self, mantle, filters=(None,), filter_args=(None,)):
        """
        Find best-fitting potential temperatures and corresponding melt
        fractions for each sample.
        """
        self.check_samples_for_fitting(mantle, filters, filter_args)
        self.individual_potential_temperatures = self.PT_to_fit.apply(
            find_sample_potential_temperature, 
            axis=1, 
            result_type="expand", 
            args=(mantle,))

    def find_suite_potential_temperature(self, mantle, find_bounds=False, filters=(None,), filter_args=(None,)):
        """
        Find best-fitting potential temperature for entire suite.
        """
        self.check_samples_for_fitting(mantle, filters, filter_args)
        Tp_fit = minimize_scalar(
            compute_suite_potential_temperature_misfit, 
            bracket=(min([lith.TSolidus(0.) for lith in mantle.lithologies]),1600.), 
            bounds=(min([lith.TSolidus(0.) for lith in mantle.lithologies]),1600.),
            args=(self.PT_to_fit, mantle),
            method="bounded")
        self.path = mantle.adiabaticMelt(Tp_fit.x, Pstart=max(mantle.solidusIntersection(Tp_fit.x))+0.01, dP=-0.01)
        self.suite_melt_fractions = self.PT_to_fit.apply(find_sample_melt_fraction, axis=1, result_type="expand", args=(self.path,))
        self.potential_temperature = Tp_fit.x
        
        if find_bounds:
            upper_points = []
            lower_points = []
            for i in range(len(self.PT_to_fit)):
                if self.PT_to_fit['T'].iloc[i] - self.suite_melt_fractions['T'].iloc[i] > 0:
                    upper_points.append(shp.Point(self.PT_to_fit['T'].iloc[i], self.PT_to_fit['P'].iloc[i]))
                elif self.PT_to_fit['T'].iloc[i] - self.suite_melt_fractions['T'].iloc[i] < 0.:
                    lower_points.append(shp.Point(self.PT_to_fit['T'].iloc[i], self.PT_to_fit['P'].iloc[i]))

        self.upper_potential_temperature, self.upper_path = find_bounding_potential_temperature(upper_points, self.potential_temperature, mantle)
        self.lower_potential_temperature, self.lower_path = find_bounding_potential_temperature(lower_points, self.potential_temperature, mantle, lower=True)

    def write_to_csv(self, outfile, write_primary=True, write_PT=True):
        """
        Write results to csv.
        """
        output_df = self.data.copy()
        if write_primary and self.primary is not None:
            output_df = pd.concat([output_df, self.primary], axis=1)
            output_df = output_df.drop([
                'SiO2_primary_wt_dry','Al2O3_primary_wt_dry',
                'FeO_primary_wt_dry','Fe2O3_primary_wt_dry',
                'MgO_primary_wt_dry','CaO_primary_wt_dry',
                'Na2O_primary_wt_dry','K2O_primary_wt_dry',
                'TiO2_primary_wt_dry','MnO_primary_wt_dry',
                'Cr2O3_primary_wt_dry','SiO2_primary_mol',
                'Al2O3_primary_mol','FeO_primary_mol',
                'Fe2O3_primary_mol','MgO_primary_mol',
                'CaO_primary_mol','Na2O_primary_mol',
                'K2O_primary_mol','TiO2_primary_mol',
                'MnO_primary_mol','Cr2O3_primary_mol',
                'H2O_primary_mol',
                'Si4O8','Al16/3O8', 'Fe4Si2O8','Fe16/3O8','Mg4Si2O8',
                'Ca4Si2O8','Na2Al2Si2O8','K2Al2Si2O8','Ti4O8','Mn4Si2O8',
                'Cr16/3O8'], axis=1)
        if write_PT and self.PT is not None:
            output_df = pd.concat([output_df, self.PT], axis=1)
        output_df.to_csv(outfile, index=False)