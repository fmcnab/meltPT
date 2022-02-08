import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from copy import deepcopy
import sys
import pyMelt as m
from scipy.optimize import minimize, minimize_scalar
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon

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
    """

    out_comp = {}
    total = sum(in_comp.values())
    for phase in in_comp:
        out_comp[phase] = in_comp[phase] / total * 100.
    return out_comp

def oxide_wt_to_cation_mole(in_comp):
    """
    Convert dictionary with oxide concentrations in weight percent to
    cation concentrations.
    """

    weights = {
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
    Convert dictionary with oxide concentrations in weight percent to
    oxide concentrations in mole percent.
    """

    weights = {
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
    Convert dictionary with oxide concentrations in mole percent to
    mole species concentrations.
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
    return 0.25324 + 0.003363*(in_comp['Mg'] + 0.33*in_comp['FeII'])

def forsterite_number(in_comp, Kd=None):
    """
    Compute forsterite number for given sample dictionary.
    """
    cation = oxide_wt_to_cation_mole(in_comp)
    if not Kd:
        Kd = compute_partition_coefficient(cation)
    return 1. / (1. + (Kd * (cation['FeII'] / cation['Mg'])))

def add_olivine(in_comp, Kd=None, dm=0.0005):
    """
    Add olivine in equilibrium with given melt composition.
    """

    # Compute melt forsterite number
    Fo = forsterite_number(in_comp, Kd=Kd)

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
    out_dict = {}
    for key in in_dict:
        out_dict[key] = np.nan
    return out_dict

def backtrack_sample_composition(df, target_Fo=0.9, Kd=None, dm=0.0005, verbose=False):
    """
    Backtrack composition to desired mantle forsterite number. Iteratively adds
    olivine in equilibrium with melt until desired composition is reached.
    """

    # Get major oxides from data frame
    oxide_wt_hydrous = {}
    for ox in MAJOR_OXIDES:
        oxide_wt_hydrous[ox] = df[ox]

    # Check Fo is below mantle Fo
    Fo = forsterite_number(oxide_wt_hydrous)
    if target_Fo-Fo < 0.005:
        oxide_wt_hydrous = fill_dict_with_nans(oxide_wt_hydrous)
        print(df.Sample + ": backtracking failed! Starting Fo above mantle Fo.")
    # Otherwise add olvine until primary Fo is reached
    else:

        if verbose:
            print("Backtracking sample %s to primary composition:" % df.Sample)

        i = 0
        dm_tot = 0.
        # while abs(target_Fo - Fo) > 1.e-15:
        while abs(target_Fo - Fo) > dm:
            oxide_wt_hydrous = add_olivine(oxide_wt_hydrous, Kd=Kd, dm=dm)
            dm_tot += dm
            Fo = forsterite_number(oxide_wt_hydrous, Kd=Kd)
            if verbose:
                print(
                    "    - iteration %d: %.2f%% olivine added, melt Fo = %.4f." %
                    (i, dm_tot/(1.+dm_tot)*100., Fo)
                    )
            i += 1
            if i>10000:
                oxide_wt_hydrous = fill_dict_with_nans(oxide_wt_hydrous)
                print(df.Sample + ": backtracking failed! Iterations: %d" % i)
                break
            # if target_Fo - Fo < 2.*dm:
            #     dm /= 1.5

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

    return primary_oxide

def backtrack_compositions(df, target_Fo=0.9, Kd=False, dm=0.0005, verbose=False):
    """
    Backtrack compositions for entire dataframe and append results.
    """
    primary = df.apply(backtrack_sample_composition, axis=1, result_type="expand", args=(target_Fo,Kd,dm,verbose))
    return pd.concat([df, primary], axis=1)

def compute_water_correction(df):
    """
    Calculate temperature correction from water content.

    Equation (1), Plank and Forsyth (2016, G-cubed).
    """
    return (
        40.4*df.H2O_primary_wt -
        2.97*df.H2O_primary_wt**2. +
        0.0761*df.H2O_primary_wt**3.
        )

def compute_CO2_correction(df):
    """
    Calculate temperature correction for samples with equlibration pressures
    above 2 GPa (assumed to be influenced by presence of CO2).

    Equation (3), Plank and Forsyth (2016, G-cubed).
    """
    return (df['SiO2_primary_wt_dry'] - 50.3) / -0.12804

def compute_temperature(df):
    """
    Calculate equlibration temperature.

    Equation (1), Plank and Forsyth (2016, G-cubed).
    """
    return (
        1264.5 +
        7.85*df['Mg4Si2O8'] +
        8545./df['Si4O8'] -
        5.96*df['Al16/3O8']
        )

def compute_pressure(df, T):
    """
    Calculate equlibration pressure.

    Equation (2), Plank and Forsyth (2016, G-cubed).
    """
    return (
        np.log(df['Si4O8']) -
        4.045 +
        0.0114*df['Fe4Si2O8'] +
        0.00052*df['Ca4Si2O8']**2. +
        0.0024*df['Mg4Si2O8']) / (-336.3/T - 0.0007*np.sqrt(T)
        )

def compute_sample_pressure_temperature(df):
    """
    Compute equilibration pressure and temperature for a given sample.
    """
    water_correction = compute_water_correction(df)
    T = compute_temperature(df) - water_correction
    P = compute_pressure(df, T)
    if P > 2.:
        T -= compute_CO2_correction(df)
        P = compute_pressure(df, T)
    return {'P': P, 'T': T - 273.15}

def compute_pressure_temperature(df):
    """
    Compute equilibration pressures and temperatures for entire dataframe
    and append results.
    """
    PT = df.apply(compute_sample_pressure_temperature, axis=1, result_type="expand")
    return pd.concat([df, PT], axis=1)

# ---- Melt Path Fitting

def melt_fraction_to_pressure_temperature(F, path):
    """
    Find pressure and temperature of given melt fraction along given melting
    path.
    """
    P = np.interp(F, path.F_total, path.P)
    T = np.interp(F, path.F_total, path.T)
    return P, T
    
def compute_sample_melt_fraction_misfit(F, df, path, full_output=True):
    """
    Compute misfit between sample and point along given melting path.
    Point on melt path is specified by melt fraction.
    """
    model_P, model_T = melt_fraction_to_pressure_temperature(F, path)
    misfit = np.sqrt( (abs(df['P']-model_P)/0.24)**2. + (abs(df['T']-model_T)/39.)**2. )
    return misfit
    
def find_sample_melt_fraction(df, path, full_output=True):
    """
    Find best-fitting melt fraction between sample and given melting path.
    
    If "full_output" is True, returns dictionary with best-fitting melt
    fraction, corresponding pressure and temperature, and misfit at minimum.
    If False, only returns misfit (for use in fitting melting paths to entire
    suite.)
    """
    if np.isnan(df['P']) or np.isnan(df['T']):
        return {'F_path_ind': np.nan, 'P_path_ind': np.nan, 'T_path_ind': np.nan, 'misfit': np.nan}
    else:
        fit = minimize_scalar(
            compute_sample_melt_fraction_misfit, 
            bounds=(0.,max(path.F_total)), 
            bracket=(0.,max(path.F_total)), 
            args=(df, path,False), 
            method="bounded")
        P, T = melt_fraction_to_pressure_temperature(fit.x, path)
        return {'F_path_ind': fit.x, 'P_path_ind': P, 'T_path_ind': T, 'misfit': fit.fun}    

def compute_sample_potential_temperature_misfit(Tp, df, mantle):
    path = mantle.AdiabaticMelt_1D(Tp, Pstart=max(mantle.solidus_intersection(Tp))+0.01, steps=101)
    fit = find_sample_melt_fraction(df, path, full_output=False)
    return fit

def find_sample_potential_temperature(df, mantle):
    if np.isnan(df['P']) or np.isnan(df['T']):
        return {'F_path_ind': np.nan, 'P_path_ind': np.nan, 'T_path_ind': np.nan, 'misfit': np.nan, 'Tp_ind': np.nan, 'path_ind': np.nan}
    else:    
        Tp_fit = minimize_scalar(
            compute_sample_potential_temperature_misfit, 
            bracket=(min([lith.TSolidus(0.) for lith in mantle.lithologies]),1600.), 
            bounds=(min([lith.TSolidus(0.) for lith in mantle.lithologies]),1600.), 
            args=(df,mantle), 
            method="bounded")
        path = mantle.AdiabaticMelt_1D(Tp_fit.x, Pstart=max(mantle.solidus_intersection(Tp_fit.x))+0.01, steps=101)
        out_dict = find_sample_melt_fraction(df, path)
        out_dict['Tp'] = Tp_fit.x 
        out_dict['path'] = path
        return out_dict
        
def compute_suite_potential_temperature_misfit(Tp, df, mantle):
    path = mantle.AdiabaticMelt_1D(Tp, Pstart=max(mantle.solidus_intersection(Tp))+0.01, steps=101)
    melt_fraction_fits = df.apply(find_sample_melt_fraction, axis=1, result_type="expand", args=(path,))
    return np.nanmean(melt_fraction_fits['misfit'])

def find_bound(points, starting_temperature, mantle, lower=False):
    
    if not lower:
        bounding_temperature = np.ceil(starting_temperature)
        adjustment = 1.
    else:
        bounding_temperature = np.floor(starting_temperature)
        adjustment = -1.
    
    inside = np.zeros(len(points))
    
    max_P = max([p.coords.xy[1][0] for p in points])
    main_path = mantle.AdiabaticMelt_1D(starting_temperature, Pstart=max(max(mantle.solidus_intersection(starting_temperature))+0.01, max_P), steps=101)

    while True:
        bounding_path = mantle.AdiabaticMelt_1D(bounding_temperature, Pstart=max(max(mantle.solidus_intersection(bounding_temperature))+0.01, max_P), steps=101)
        bounds = np.vstack(( 
            np.column_stack(( main_path.T, main_path.P )),
            np.column_stack(( bounding_path.T[::-1], bounding_path.P[::-1] ))
            ))
        poly = Polygon(bounds)
        for i,p in enumerate(points):
            if poly.contains(p):
                inside[i] = 1.
        if sum(inside) / len(points) > 2./3.:
            break
        else:
            bounding_temperature += adjustment
            
    return bounding_temperature, bounding_path


def combine(df):
    return {'fit': df.to_numpy().all()}


# ---- Suite class

class Suite:

    def __init__(self, input_csv, src_FeIII_totFe=0.19, min_SiO2=0., min_MgO=0.):
        self.data = parse_csv(input_csv, src_FeIII_totFe=src_FeIII_totFe, min_SiO2=min_SiO2, min_MgO=min_MgO)
        self.primary = None
        self.PT = None

    def backtrack_compositions(self, target_Fo=0.9, Kd=False, dm=0.0005, verbose=False):
        """
        Backtrack compositions for entire dataframe and append results.
        """
        self.primary = self.data.apply(
            backtrack_sample_composition,
            axis=1,
            result_type="expand",
            args=(target_Fo,Kd,dm,verbose)
            )

    def compute_pressure_temperature(self):
        """
        Compute equilibration pressures and temperaturesfor entire dataframe
        and append results.
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
        self.path = mantle.AdiabaticMelt_1D(Tp_fit.x, Pstart=max(mantle.solidus_intersection(Tp_fit.x))+0.01, steps=101)
        self.suite_melt_fractions = self.PT_to_fit.apply(find_sample_melt_fraction, axis=1, result_type="expand", args=(self.path,True))
        self.potential_temperature = Tp_fit.x
        
        if find_bounds:
            upper_points = []
            lower_points = []
            for i in range(len(self.PT_to_fit)):
                if self.PT_to_fit['T'].iloc[i] - self.suite_melt_fractions['T_path_ind'].iloc[i] > 0:
                    upper_points.append(Point(self.PT_to_fit['T'].iloc[i], self.PT_to_fit['P'].iloc[i]))
                elif self.PT_to_fit['T'].iloc[i] - self.suite_melt_fractions['T_path_ind'].iloc[i] < 0.:
                    lower_points.append(Point(self.PT_to_fit['T'].iloc[i], self.PT_to_fit['P'].iloc[i]))

        self.upper_potential_temperature, self.upper_path = find_bound(upper_points, self.potential_temperature, mantle)
        self.lower_potential_temperature, self.lower_path = find_bound(lower_points, self.potential_temperature, mantle, lower=True)

