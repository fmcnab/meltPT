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

class BacktrackOlivineFractionation:
    
    def __init__(
        self, Kd=None, dm=0.0005, verbose=False, 
        max_olivine_addition=0.3):
    
        self.fixed_Kd = Kd
        self.dm = dm
        self.verbose = verbose
        self.max_olivine_addition = max_olivine_addition

    @property
    def cation_mole(self):
        """
        Convert oxide concentrations in weight percent to cation concentrations.
            
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
        for phase in self.oxide_wt_hydrous:
            out_comp[keys[phase]] = self.oxide_wt_hydrous[phase] / weights[phase]
        return normalise(out_comp)
    
    @property
    def Kd(self):
        """
        Compute the partition coefficient.

        If a fixed value has been specified, it is returned. Otherwise,
        calculate as function of Mg & FeII content, using expression from 
        Tamura et al. (2000, J. Pet), via Lee et al. (2009) spreadsheet.
        """

        if self.fixed_Kd:
            return self.fixed_Kd
        else:
            cation = self.cation_mole
            return 0.25324 + 0.0033663*(cation['Mg'] + 0.33*cation['FeII'])

    @property
    def Fo(self):
        """
        Compute forsterite number from oxide weight compositions.

        Returns
        -------
        Fo : float
            The calculated forsterite number.
        """

        cation = self.cation_mole
        Fo = 1. / (1. + (self.Kd * cation['FeII'] / cation['Mg']))
        return Fo
        
    def add_olivine(self):
        """
        Add olivine in equilibrium with given melt composition.
        
        Returns
        -------
        out_comp : dict
            The updated concentrations.
        """

        Fo = self.Fo
        
        # Compute olivine composition in equilibrium with melt
        oxide_wt_olivine = {
            'FeO': 2. * (1 - Fo) * 71.85,
            'MgO': 2. * Fo * 40.3,
            'SiO2': 60.08
            }
        oxide_wt_olivine = normalise(oxide_wt_olivine)

        # loop over phases adding olivine
        out_comp = {}
        for phase in self.oxide_wt_hydrous:
            if phase == 'SiO2' or phase == 'FeO' or phase == 'MgO':
                out_comp[phase] = (self.oxide_wt_hydrous[phase] + self.dm*oxide_wt_olivine[phase]) / (1. + self.dm)
            else:
                out_comp[phase] = self.oxide_wt_hydrous[phase] / (1. + self.dm)

        return normalise(out_comp)


    def backtrack_sample_composition(self, df, return_all=False):
        """
        Backtrack composition to desired mantle forsterite number.
    
        Iteratively adds olivine in equilibrium with melt until desired composition
        is reached.
    
        Parameters
        ----------
        df : pandas dataframe
            Dataframe containing the initial composition to be backtracked and the
            target forsterite number. Should contain only one row. To use with a 
            multi-row dataframe use df.apply().
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
        self.oxide_wt_hydrous = {}
        for ox in MAJOR_OXIDES:
            self.oxide_wt_hydrous[ox] = df[ox]
        self.oxide_wt_hydrous = normalise(self.oxide_wt_hydrous)
        # Fo = self.compute_forsterite_number()

        # Check Fo is below mantle Fo
        if df['src_Fo']-self.Fo < 0.001:
            self.oxide_wt_hydrous = fill_dict_with_nans(self.oxide_wt_hydrous)
            dm_tot = np.nan
            message = df.Sample + ": backtracking failed! Starting Fo above mantle Fo."
            warnings.warn(message)
        # Otherwise add olvine until primary Fo is reached
        else:

            if self.verbose:

                print("Backtracking sample %s to primary composition:" % df.Sample)

            dm_tot = 0.
            composition_through_addition = []
            while df['src_Fo'] - self.Fo > 0.0002:
                
                self.oxide_wt_hydrous = self.add_olivine()
                composition_through_addition.append(self.oxide_wt_hydrous)
                dm_tot += self.dm
                if self.verbose:
                    print(
                        "    - %.2f%% olivine added, melt Fo = %.4f, Kd = %.4f." %
                        (dm_tot/(1.+dm_tot)*100., self.Fo, self.Kd)
                        )
                if dm_tot/(1. + dm_tot) > self.max_olivine_addition:
                    self.oxide_wt_hydrous = fill_dict_with_nans(self.oxide_wt_hydrous)
                    message = (
                        df.Sample + ": backtracking failed! Olivine addition exceeding %d%%" 
                        % (self.max_olivine_addition*100.))
                    warnings.warn(message)
                    break

        # Package up
        primary_oxide = {}
        for phase in self.oxide_wt_hydrous:
            primary_oxide[phase + "_primary_wt"] = self.oxide_wt_hydrous[phase]
            
        # Add amount olivine added
        primary_oxide['ol_added'] = dm_tot / (1. + dm_tot)

        if not return_all:
            return primary_oxide
        else:
            return primary_oxide, composition_through_addition

