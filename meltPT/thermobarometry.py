"""
===============
thermobarometry
===============

Perform thermobarometric calculations.

"""


import warnings
import numpy as np
import pandas as pd
import sympy as sym
import inspect

def Mg_num(df):
    """
    Calculate Mg#.
        
    Parameters
    ----------
    df : pandas dataframe 
       sample compositions in anhydrous mol%.
                
    Returns
    -------
    Mg_num : float
        The Mg#.
    """
    Mg_num = df['MgO_primary_mol_dry'] / (df['MgO_primary_mol_dry'] + 
                                            df['FeO_primary_mol_dry'])
    return Mg_num

def NaK_num(df):
    """
    Calculate NaK#.
    
    Must previously run compute_components_species(),
    compute_components_cation(), compute_components_compound()
    or compute components_per_ox().
    
    Parameters
    ----------
    df : pandas dataframe with sample compositions in anhydrous mole%
            
    Returns
    -------
    NaK_num : float
        The NaK#.        
    """            
    NaK_num = (
        (df['Na2O_primary_wt_dry'] + df['K2O_primary_wt_dry']) / 
        (df['Na2O_primary_wt_dry'] + df['K2O_primary_wt_dry'] + 
            df['CaO_primary_wt_dry'])
        )
    return NaK_num

def normalize_v2(df, list_parameters, input_suffix, output_suffix, end_sum):
    """
    Normalise sample compositions so they sum to a given value.
    
    Concentrations in list_parameters are normalised so that they sum to 
    value given by end_sum. New values added to df with list_parameter
    name appended by output_suffix.
    
    Parameters
    ----------
    df : pandas dataframe
        Should contain only one row. To use with a multi-row dataframe use
        df.apply().
    list_parameters : list 
        Dataframe columns to be normalised.
    input_suffix : str
        Original suffix for columns to use.
    output_suffix : str
        Suffix to be added to create new column name.
    end_sum : float
        Value that columns should total after normalisation.

    Returns
    -------
    df : pandas datframe
        The updated dataframe containing normalised values.
    """
    SUM = 0.
    for i in range(0, len(list_parameters), 1):
        parameter = list_parameters[i]
        SUM = df[parameter + input_suffix] + SUM

    for i in range(0, len(list_parameters), 1):
        parameter = list_parameters[i]        
        df[parameter + output_suffix] = end_sum * df[parameter + input_suffix] / SUM
        
    return df

def compute_DMg(df):
    """
    Calculate DMg.

    Follows method outlined in Beattie, (1993), Contrib. Min. and Pet.
    Must previously run compute_components_cation().
    
    Parameters
    ----------
    df : pandas dataframe
        Dataframe containing sample compositions in anhydrous mole%.
            
    Returns
    -------
    DMg : float
        Calculated DMg.
    """ 
    
    Sum_A = (df["FeO_primary_mol_dry"]*0.279)
    Sum_A += (df["MnO_primary_mol_dry"]*0.259)
    Sum_A += (df["MgO_primary_mol_dry"]*1.)
    Sum_A += (df["CaO_primary_mol_dry"]*0.0056)
    Sum_A += (df["NiO_primary_mol_dry"]*3.346)
    #####
    Sum_B = (df["FeO_primary_mol_dry"]*0.031)
    Sum_B -= (df["MnO_primary_mol_dry"]*0.049)
    Sum_B += (df["CaO_primary_mol_dry"]*0.0135)
    Sum_B += (df["NiO_primary_mol_dry"]*-3.665)
    Sum_B += 0.0001*(df["TiO2_primary_mol_dry"] + 
                     df["Al2O3_primary_mol_dry"] +
                     df["Cr2O3_primary_mol_dry"] +
                     df["Fe2O3_primary_mol_dry"] + 
                     df["Na2O_primary_mol_dry"] + 
                     df["K2O_primary_mol_dry"] + 
                     df["P2O5_primary_mol_dry"])
    #####
    DMg = (2./3. - Sum_B) / Sum_A  

    return DMg

def compute_minerals(df):
    """
    Calculate proportions of mineral components.
    
    Using equations listed in Grove (1982, Contrib. Min. Pet., 80:160-182), 
    compute stable mineral assemblage.
        
    Parameters
    ----------
    df : pandas dataframe 
        Datframe containing sample compositions. Expects columns with suffix 
        '_primary_mol_dry'; should be mol% anhydrous compositions.
        Should contain only one row. To use with a multi-row dataframe use
        df.apply().
                
    Returns
    -------
    MINS: numpy array
        Predicted mineral proportions.
        0: Olivine, 1: cpx, 2: plg, 3: qz, 4: opx, 5: spl, 6: ap.
    """      
    Sum = df["SiO2_primary_mol_dry"] + df["TiO2_primary_mol_dry"] + (df["Al2O3_primary_mol_dry"]*2.) 
    Sum += (df["Cr2O3_primary_mol_dry"]*2) + df["FeO_primary_mol_dry"]
    Sum += df["MnO_primary_mol_dry"] + df["MgO_primary_mol_dry"] 
    Sum += df["CaO_primary_mol_dry"] + (df["Na2O_primary_mol_dry"]*2) + (df["K2O_primary_mol_dry"]*2)
    Sum += df["NiO_primary_mol_dry"] + (df["P2O5_primary_mol_dry"]*2)
    ####    
    sio2 = 100. * df["SiO2_primary_mol_dry"] / Sum
    tio2 = 100. * df["TiO2_primary_mol_dry"] / Sum
    al2o3 = 200. * df["Al2O3_primary_mol_dry"] / Sum
    cr2o3 = 200. * df["Cr2O3_primary_mol_dry"] / Sum
    feo = 100. * df["FeO_primary_mol_dry"] / Sum
    mgo = 100. * df["MgO_primary_mol_dry"] / Sum
    cao = 100. * df["CaO_primary_mol_dry"] / Sum
    na2o = 200. * df["Na2O_primary_mol_dry"] / Sum
    k2o = 200. * df["K2O_primary_mol_dry"] / Sum
    p2o5 = 200. * df["P2O5_primary_mol_dry"] / Sum
    ####     
    Total = sio2 - cao - 2*(k2o+na2o) + cr2o3 + tio2 + (2*p2o5)
    ####    
    ol = (0.5*(feo+mgo+0.5*(al2o3-k2o-na2o)-cao-tio2-cr2o3+5/3*p2o5))/Total*2
    ####
    qz = (sio2-0.5*(feo+mgo)-1.5*cao-0.25*al2o3-2.75*(k2o+na2o)+0.5*cr2o3+0.5*tio2+2.5*p2o5)/Total
    ####  
    cpx = (cao-0.5*al2o3+0.5*(k2o+na2o)-5/3*p2o5)/Total*3
    ####    
    plg = (0.5*(al2o3+na2o-k2o))/Total*4
    ####    
    opx = k2o/Total*4
    #### 
    spl = (cr2o3+tio2)/Total*1.5
    ####   
    ap = p2o5/3/Total*6
    ####     
    Min_Tot = ol + qz + plg + cpx + opx + spl + ap
    ####       
    MINS = [ol/Min_Tot, cpx/Min_Tot, plg/Min_Tot, qz/Min_Tot, opx/Min_Tot, spl/Min_Tot, ap/Min_Tot]

    return MINS

def compute_parameter_BK21(df, OL, P, regression_params):
    """
    Calculate parameter values using equation
    of Brown Krein et al., (2021)
        
    Parameters
    ----------
    df : pandas dataframe
        sample compositions in anhydrous wt% need suffix '_primary_wt_dry'
        sample compositions in anhydrous mol% need suffix '_primary_mol_dry'.
        Should contain only one row. To use with a multi-row dataframe use
        df.apply().
    OL : flota
        Olivine concentration.
    P: float
        Pressure in GPa.
    regression_params: array
        Array containing eight calibrated parameter values.
                
    Returns
    -------
    result: float
        either T or predicted mineral percentages.
    """ 
    
    result = regression_params[0]
    result += regression_params[1] * OL 
    result += regression_params[2] * P
    result += regression_params[3] * Mg_num(df)
    result += regression_params[4] * NaK_num(df)
    result += regression_params[5] * df['TiO2_primary_wt_dry'] 
    result += regression_params[6] * df['K2O_primary_wt_dry']
    result += regression_params[7] * (df['CaO_primary_wt_dry']/df['Al2O3_primary_wt_dry'])        
    return result

def compute_parameter_TGK12(df, MINS, P, regression_params):
    """
    Calculate parameter values using equation
    of Till et al., (2012)
        
    Parameters
    ----------
    df : pandas dataframe
        sample compositions in anhydrous wt% need suffix '_primary_wt_dry'
        sample compositions in anhydrous mol% need suffix '_primary_mol_dry'.
        Should contain only one row. To use with a multi-row dataframe use
        df.apply().
    MINS : array
        Concentrations of each mineral.
    P: float
        Pressure in GPa.
    regression_params: array
        one row array of 8 calibrated prameter values.
                
    Returns
    -------
    result: float
        either T or predicted mineral percentages
    """ 
    
    result = regression_params[0]
    result += regression_params[1] * (MINS[0]/(MINS[0]+MINS[1]+MINS[2]+MINS[3]))
    result += regression_params[2] * P
    result += regression_params[3] * (1. - Mg_num(df))
    result += regression_params[4] * NaK_num(df)
    result += regression_params[5] * df['TiO2_primary_wt_dry'] 
    result += regression_params[6] * df['K2O_primary_wt_dry']    
    return result

def compute_components_species(df):
    """
    Calculate mole species proportions.
    
    Use equations listed in Appendix A of Lee et al (2009, G-cubed) to compute
    proportions of each mole species.
        
    Parameters
    ----------
    df : pandas dataframe
        Dataframe containing sample compositions in wt%.
        Should contain only one row. To use with a multi-row dataframe use
        df.apply().
                
    Returns
    -------
    df : pandas dataframe
        Original dataframe with additional columns:
          - anhydrous wt% with suffix '_primary_wt_dry'
          - oxide mol% with suffix '_primary_mol'
          - species proportions in mol%, e.g., 'Si4O8'
          - anhydrous species proportions in mol%, e.g., 'Si4O8_dry'.
    """   
    
    MAJOR_OXIDES = ['SiO2','Al2O3','FeO','Fe2O3','MgO','CaO','Na2O','K2O',
        'TiO2','MnO','Cr2O3', 'P2O5', 'H2O']
    OXIDES_DRY = ['SiO2','Al2O3','FeO','Fe2O3','MgO','CaO','Na2O','K2O',
        'TiO2','MnO','Cr2O3', 'P2O5']
    OXIDE_WEIGHT = [60.08,101.96,71.84,159.69,40.3,56.08,61.98,
        94.2,79.86,70.94,151.99,141.942524,18.014680000000002]
    SPECIES = ['Si4O8','Al16/3O8','Fe4Si2O8','Fe16/3O8','Mg4Si2O8',
        'Ca4Si2O8','Na2Al2Si2O8','K2Al2Si2O8','Ti4O8','Mn4Si2O8','Cr16/3O8',
        'P16/5O8']     
    SPECIES_DRY = ['Si4O8','Al16/3O8','Fe4Si2O8','Fe16/3O8','Mg4Si2O8',
        'Ca4Si2O8','Na2Al2Si2O8','K2Al2Si2O8','Ti4O8','Mn4Si2O8','Cr16/3O8',
        'P16/5O8']         

    normalize_v2(df, OXIDES_DRY, '_primary_wt', '_primary_wt_dry', 100.)

    for i in range(0, len(MAJOR_OXIDES), 1):
        OXIDE = MAJOR_OXIDES[i]
        df[OXIDE + '_primary_mol'] = df[OXIDE + '_primary_wt'] / OXIDE_WEIGHT[i]

    normalize_v2(df, MAJOR_OXIDES, '_primary_mol', '_primary_mol', 100.)
        
    df['Si4O8'] = (df['SiO2_primary_mol'] - 0.5*(df['FeO_primary_mol']
        +df['MgO_primary_mol']+df['CaO_primary_mol']+df['MnO_primary_mol'])
        - df['Na2O_primary_mol'] - df['K2O_primary_mol']) * 0.25
    df['Al16/3O8'] = (df['Al2O3_primary_mol'] - df['Na2O_primary_mol']) * (3./8.)
    df['Fe4Si2O8'] = df['FeO_primary_mol'] * 0.25
    df['Fe16/3O8'] = df['Fe2O3_primary_mol'] * (3./8.)
    df['Mg4Si2O8'] = df['MgO_primary_mol'] * 0.25
    df['Ca4Si2O8'] = df['CaO_primary_mol'] * 0.25
    df['Na2Al2Si2O8'] = df['Na2O_primary_mol']
    df['K2Al2Si2O8'] = df['K2O_primary_mol']
    df['Ti4O8'] = df['TiO2_primary_mol'] * 0.25
    df['Mn4Si2O8'] = df['MnO_primary_mol'] * 0.25
    df['Cr16/3O8'] = df['Cr2O3_primary_mol']  * (3./8.)
    df['P16/5O8'] = df['P2O5_primary_mol'] * 0.625
    df['H16O8'] = df['H2O_primary_mol'] * 0.125

    normalize_v2(df, SPECIES, '', '', 100.)
    normalize_v2(df, SPECIES_DRY, '', '_dry', 100.)

    return

def compute_components_compound(df):
    """
    Calculate proportions of oxides in mol% assuming all Fe is in FeO
    and the sample is anhydrous.
        
    Parameters
    ----------
    df : pandas dataframe
        The sample compositions in wt%.
        Should contain only one row. To use with a multi-row dataframe use
        df.apply().
                
    Returns
    -------
    df : pandas dataframe
        Original dataframe with additional columns:
          - anhydrous wt% with suffix '_primary_wt_dry'
          - anhydrous oxide mol% with suffix '_primary_mol_dry'.
    """ 
        
    MAJOR_OXIDES = ['SiO2','Al2O3','FeO','MgO','CaO','Na2O','K2O','TiO2',
        'MnO','Cr2O3','P2O5', 'NiO']
    OXIDE_WEIGHT = [60.08,101.96,71.84,40.3,56.08,61.98,
        94.2,79.86,70.94,151.99,141.945,74.71]

    df['FeO_primary_wt'] = df['FeO_primary_wt'] + 0.8998 * df['Fe2O3_primary_wt']
    df['Fe2O3_primary_wt'] = 0.

    normalize_v2(df, MAJOR_OXIDES, '_primary_wt', '_primary_wt_dry', 100.)

    for i in range(0, len(MAJOR_OXIDES), 1):
        OXIDE = MAJOR_OXIDES[i]
        df[OXIDE + '_primary_mol_dry'] = df[OXIDE + '_primary_wt_dry'] / OXIDE_WEIGHT[i]
    
    normalize_v2(df, MAJOR_OXIDES, '_primary_mol_dry', '_primary_mol_dry', 100.)      
    
    return

def compute_components_cation(df):
    """
    Calculate proportions of cations in mol% assuming the sample is anhydrous.
        
    Parameters
    ----------
    df : pandas dataframe
        Dataframe containing sample compositions in wt%.
        Should contain only one row. To use with a multi-row dataframe use
        df.apply().
                
    Returns
    -------
    df : pandas dataframe
        Original dataframe with additional columns:
          - anhydrous wt% with suffix '_primary_wt_dry'
          - cation mol% with suffix '_primary_mol'
          - anhydrous cation mol% with suffix '_primary_mol_dry'.
    """
        
    MAJOR_OXIDES = ['SiO2','Al2O3','FeO','Fe2O3','MgO','CaO','Na2O','K2O',
        'TiO2','MnO','Cr2O3', 'P2O5', 'NiO', 'H2O']
    CAT_WEIGHT = [60.08,50.98,71.85,159.96/2.,40.3,56.08,30.99,
        47.1,79.9,70.94,151.99/2.,141.945/2.,74.71, 18.014680000000002]    
    MAJOR_OXIDES_DRY = ['SiO2','Al2O3','FeO','Fe2O3','MgO','CaO','Na2O','K2O',
        'TiO2','MnO','Cr2O3', 'P2O5', 'NiO']
    CAT_WEIGHT_DRY = [60.08,50.98,71.85,159.96/2.,40.3,56.08,30.99,
        47.1,79.9,70.94,151.99/2.,141.945/2.,74.71] 

    normalize_v2(df, MAJOR_OXIDES_DRY, '_primary_wt', '_primary_wt_dry', 100.)    
    normalize_v2(df, MAJOR_OXIDES, '_primary_wt', '_primary_wt', 100.)

    for i in range(0, len(MAJOR_OXIDES_DRY), 1):
        OXIDE = MAJOR_OXIDES_DRY[i]
        df[OXIDE + '_primary_mol_dry'] = df[OXIDE + '_primary_wt_dry'] / CAT_WEIGHT_DRY[i]

    for i in range(0, len(MAJOR_OXIDES), 1):
        OXIDE = MAJOR_OXIDES[i]
        df[OXIDE + '_primary_mol'] = df[OXIDE + '_primary_wt'] / CAT_WEIGHT[i]

    normalize_v2(df, MAJOR_OXIDES, '_primary_mol', '_primary_mol', 1.)     
    normalize_v2(df, MAJOR_OXIDES_DRY, '_primary_mol_dry', '_primary_mol_dry', 1.)  
    
    return

def compute_components_per_oxygen(df):
    """
    Calculate proportions of oxides in mol% assuming all Fe is in FeO
    and the sample is anhydrous.
    Normalise to each to one oxygen per compound.
        
    Parameters
    ----------
    df : pandas dataframe
        Dataframe containing sample compositions in wt%.
        Should contain only one row. To use with a multi-row dataframe use
        df.apply().
                
    Returns
    -------
    df : pandas dataframe
        Original dataframe with additional columns:
          - anhydrous wt% with suffix '_primary_wt_dry'
          - volatile free oxide mol content per unit oxygen with suffix '_primary_mol_dry'.
    """ 
    MAJOR_OXIDES = ['SiO2','Al2O3','FeO','MgO','CaO','Na2O','K2O',
        'TiO2','MnO','Cr2O3', 'P2O5', 'NiO', 'H2O', 'CO2'] 
    MAJOR_OXIDES_DRY = ['SiO2','Al2O3','FeO','MgO','CaO','Na2O','K2O','TiO2',
        'MnO','Cr2O3','P2O5', 'NiO']
    OXIDE_WEIGHT = [60.08,101.96,71.84,40.3,56.08,61.98,
        94.2,79.86,70.94,151.99,141.945,74.71]    
    #OXIDE_WEIGHT = [60.843, 101.9612, 71.8444, 40.3039, 56.0774, 61.979,
    #                94.196, 79.8658, 70.9374, 151.9904, 141.945, 74.6928]
    ANO_N = [2., 3., 1., 1., 1., 1., 1., 2., 1., 3., 5., 1.]
    CAT_N = [1., 2., 1., 1., 1., 2., 2., 1., 1., 2., 2., 1.]
    
    df['FeO_primary_wt'] = df['FeO_primary_wt'] + 0.8998 * df['Fe2O3_primary_wt']
    df['Fe2O3_primary_wt'] = 0.

    normalize_v2(df, MAJOR_OXIDES_DRY, '_primary_wt', '_primary_wt_dry', 100.)    
    normalize_v2(df, MAJOR_OXIDES, '_primary_wt', '_primary_wt', 100.)

    SUM = 0.
    for i in range(0, len(MAJOR_OXIDES_DRY), 1):
        OXIDE = MAJOR_OXIDES_DRY[i]
        df[OXIDE + '_primary_mol_dry'] = df[OXIDE + '_primary_wt_dry'] / OXIDE_WEIGHT[i] * ANO_N[i]
        SUM += df[OXIDE + '_primary_mol_dry']
    O_FACTOR = 1. / SUM   

    for i in range(0, len(MAJOR_OXIDES_DRY), 1):
        OXIDE = MAJOR_OXIDES_DRY[i]        
        df[OXIDE + '_primary_mol_dry'] = df[OXIDE + '_primary_mol_dry'] * O_FACTOR / ANO_N[i] * CAT_N[i]
    
    return

class PF16:
    
    def __init__(self, df):
        self.df = df
        self.P_err = 0.24
        self.T_err = 39.

    def compute_water_correction(self):
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
        
        correction = (
            40.4*self.df['H2O_primary_wt'] - 
            2.97*self.df['H2O_primary_wt']**2. + 
            0.0761*self.df['H2O_primary_wt']**3.
            )
        return correction

    def compute_CO2_correction(self):
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
        correction = (self.df['SiO2_primary_wt_dry'] - 50.3) / -0.12804
        return correction

    def compute_temperature(self):
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
            7.85*self.df['Mg4Si2O8_dry'] +
            8545./self.df['Si4O8_dry'] -
            5.96*self.df['Al16/3O8_dry']
            )
        return temperature

    def compute_pressure1(self, T):
        """
        Calculate equlibration pressure.

        Equation (2), Plank and Forsyth (2016, G-cubed).

        Parameters
        ----------
        df : pandas dataframe
            The sample compositions.
            Should contain only one row. To use with a multi-row dataframe use
            df.apply().
            
        Returns
        -------
        pressure : pandas dataframe
            The calculated pressure(s).
        """
        pressure = (
            np.log(self.df['Si4O8_dry']) -
            4.045 +
            0.0114*self.df['Fe4Si2O8_dry'] +
            0.00052*self.df['Ca4Si2O8_dry']**2. +
            0.0024*self.df['Mg4Si2O8_dry']) / (-336.3/T - 0.0007*np.sqrt(T)
            )
        return pressure

    def compute_pressure(self, T):
        """
        Calculate equlibration pressure.

        Equation (2), Plank and Forsyth (2016, G-cubed).

        Parameters
        ----------
        df : pandas dataframe
            The sample compositions.
            Should contain only one row. To use with a multi-row dataframe use
            df.apply().
        T : float
            Temperature in oC.
            
        Returns
        -------
        pressure : pandas dataframe
            The calculated pressure(s).
        """
        compute_components_species(self.df)
        
        P = (
            np.log(self.df['Si4O8_dry']) -
            4.045 +
            0.0114*self.df['Fe4Si2O8_dry'] +
            0.00052*self.df['Ca4Si2O8_dry']**2. +
            0.0024*self.df['Mg4Si2O8_dry']) / (-336.3/(T+273.15) - 0.0007*np.sqrt(T+273.15)
            )
                                               
        return {'P': P, 'T': T, 'P_err': self.P_err, 'T_err': T * self.P_err/P}

    def compute_pressure_temperature(self):
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
            Associated errors for pressure and temperature.
        """
        compute_components_species(self.df)
        water_correction = self.compute_water_correction()
        T = self.compute_temperature() - water_correction
        P = self.compute_pressure1(T)
        if P > 2.:
            T -= self.compute_CO2_correction()
            P = self.compute_pressure1(T)
        out = {'P': P, 'T': T - 273.15, 'P_err': self.P_err, 'T_err': self.T_err}
        return out

class L09:
        
    def __init__(self, df):
        self.df = df
        self.P_err = 0.2
        self.T_err = 40.
        
    def compute_pressure1(self, T):
        """
        Calculate equlibration temperature.

        Equation (2), Lee et al (2009, G-cubed).
        
        Parameters
        ----------
        df : pandas dataframe
            The sample compositions in species mol%.
            Should contain only one row. To use with a multi-row dataframe use
            df.apply().
        T : float
            temperature in oC
            
        Returns
        -------
        pressure : pandas dataframe
            The calculated pressure(s).
        """
        
        pressure = (
            (np.log(self.df['Si4O8']) - 
                4.019 + 
                0.0165*self.df['Fe4Si2O8'] + 
                0.0005*(self.df['Ca4Si2O8']**2.)) /
            (-770*(T**(-1.)) + 0.0058*(T**0.5) - 0.003*self.df['H16O8'])
            )
        return pressure
    
    def compute_pressure(self, T):
        """
        Calculate equlibration temperature.

        Equation (2), Lee et al (2009, G-cubed).
        
        Parameters
        ----------
        df : pandas dataframe
            The sample compositions in species mol%.
            Should contain only one row. To use with a multi-row dataframe use
            df.apply().
        T : float
            temperature in oC
            
        Returns
        -------
        pressure : pandas dataframe
            The calculated pressure(s).
        """
        compute_components_species(self.df)
        P = (
            (np.log(self.df['Si4O8']) - 
                4.019 + 
                0.0165*self.df['Fe4Si2O8'] + 
                0.0005*(self.df['Ca4Si2O8']**2.)) /
            (-770*((T+273.15)**(-1.)) + 0.0058*(T**0.5) - 0.003*self.df['H16O8'])
            )  
        return {'P': P, 'T': T, 'P_err': self.P_err, 'T_err': T * self.P_err/P}
    
    def compute_temperature(self):
        """
        Calculate equlibration temperature.

        Equation (3), Lee et al (2009, G-cubed).
        
        Parameters
        ----------
        df : pandas dataframe
            The sample compositions.
            Should contain only one row. To use with a multi-row dataframe use
            df.apply().
            
        Returns
        -------
        temperature : pandas dataframe
            The calculated temperature(s).
        """
        
        temperature = (
            916.45 + 
            (13.68 * self.df['Mg4Si2O8']) + 
            (4580. / self.df['Si4O8']) -
            (0.509 * self.df['H16O8'] * self.df['Mg4Si2O8'])
            )
        return temperature
        
    def compute_pressure_temperature(self):
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
            Associated errors for pressure and temperature.
        """
        
        compute_components_species(self.df)
        T = self.compute_temperature()
        P = self.compute_pressure1((T+273.15))
        return {'P': P, 'T': T, 'P_err': self.P_err, 'T_err': self.T_err}

class TGK12_SPL:
    
    def __init__(self, df):
        self.df = df
        self.P_err = 0.15
        self.T_err = 10.78

    def compute_temperature(self, P):
        """
        Compute equilibration pressure and temperature for a given sample.
        
        Assumes spinel is stable mantle phase.
        Lines (3) and (4) of Table (5), Till et al (2012, JGR).
    
        Parameters
        ----------
        df : pandas dataframe
            Dataframe containing the sample primary composition.
            Should contain only one row. To use with a multi-row dataframe use
            df.apply().
        P : 
            Pressure in GPa.
        Returns
        -------
        out : dict
            The equilibration pressure and temperature result.
            Temperature is in degrees celcius.
            Pressure is in GPa.
            Associated errors for pressure and temperature.
        """
        
        compute_components_compound(self.df)       
        MINS = compute_minerals(self.df)
        T = compute_parameter_TGK12(self.df, MINS, P, [1212,0.,119.9,-97.33,-87.76,3.44,-4.58])  
        return {'P': P, 'T': T, 'P_err': P * self.T_err/T, 'T_err': self.T_err}         
        

    def compute_pressure_temperature(self):    
        """
        Compute equilibration pressure and temperature for a given sample.
        
        Assumes spinel is stable mantle phase.
        Lines (3) and (4) of Table (5), Till et al (2012, JGR).
    
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
            Associated errors for pressure and temperature.
        """
        
        compute_components_compound(self.df)       
        MINS = compute_minerals(self.df)
        P = compute_parameter_TGK12(self.df, MINS, 0., [-0.862,9.471,0.,-2.383,2.922,0.218,-0.146])
        T = compute_parameter_TGK12(self.df, MINS, P, [1212,0.,119.9,-97.33,-87.76,3.44,-4.58])  
        return {'P': P, 'T': T, 'P_err': self.P_err, 'T_err': self.T_err}   
    
class TGK12_PLG:
 
    def __init__(self, df):
        self.df = df
        self.P_err = 0.08
        self.T_err = 11.7
    
    def compute_temperature(self, P):
        """
        Compute equilibration pressure and temperature for a given sample.
        
        Assumes plagioclase is stable mantle phase.
        Line (2) of Table (5), Till et al (2012, JGR).
    
        Parameters
        ----------
        df : pandas dataframe
            Dataframe containing the sample primary composition.
            Should contain only one row. To use with a multi-row dataframe use
            df.apply().
        P : float
            Pressure in GPa.
        Returns
        -------
        out : dict
            The equilibration pressure and temperature result.
            Temperature is in degrees celcius.
            Pressure is in GPa.
            Associated errors for pressure and temperature.
        """
        
        compute_components_compound(self.df)
        MINS = compute_minerals(self.df)
        T = compute_parameter_TGK12(self.df, MINS, P, [1216,0.,104.4,-72.83,-194.9,24.08,-1.55])        
        return {'P': P, 'T': T, 'P_err': P * self.T_err/T, 'T_err': self.T_err}        

    def compute_pressure_temperature(self):
        """
        Compute equilibration pressure and temperature for a given sample.
        
        Assumes plagioclase is stable mantle phase.
        Line (1 and 2) of Table (5), Till et al (2012, JGR).
    
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
            Associated errors for pressure and temperature.
        """
        
        compute_components_compound(self.df)
        MINS = compute_minerals(self.df)
        P = compute_parameter_TGK12(self.df, MINS, 0., [-1.64,12.94,0.,-2.363,3.51,0.152,-0.176])
        T = compute_parameter_TGK12(self.df, MINS, P, [1216,0.,104.4,-72.83,-194.9,24.08,-1.55])        
        return {'P': P, 'T': T, 'P_err': self.P_err, 'T_err': self.T_err}  

class G13:

    def __init__(self, df):
        self.df = df
        self.P_err = 0.24
        self.T_err = 24.

    def compute_pressure(self):
        """
        Calculate equlibration pressure.

        Line (1) of Table (4), Grove et al (2013, Contrib Min Pet).

        Parameters
        ----------
        df : pandas dataframe
            The sample compositions.
            
        Returns
        -------
        pressure : float
            The calculated pressure(s).
        """
        compute_components_compound(self.df) 
        MINS = compute_minerals(self.df)
        pressure = -17.3
        pressure += 146.4 * MINS[0]
        pressure -= 38.4 * (1. - Mg_num(self.df)) 
        pressure += 19.6 * NaK_num(self.df)
        pressure += 4.81 * self.df['P2O5_primary_wt_dry'] 
        pressure *= 0.1
        return pressure

    def compute_temperature1(self, P):
        """
        Calculate equlibration temperature.

        Line (6) of Table (4), Grove et al (2013, Contrib Min Pet).
        
        Parameters
        ----------
        df : pandas dataframe
            The sample compositions.
        P: pressure in GPa
            
        Returns
        -------
        temperature : float
            The calculated temperature(s) in oC.
        """
        temperature = 1313.67
        temperature += 8.423 * ((P*10.)-1.)
        temperature -= 149.92 * (1. - Mg_num(self.df)) 
        temperature += 55.02 * NaK_num(self.df)
        temperature -= 59.69 * self.df['P2O5_primary_wt_dry'] 
        return temperature

    def compute_temperature(self, P):
        """
        Calculate equlibration temperature.

        Line (6) of Table (4), Grove et al (2013, Contrib Min Pet).
        
        Parameters
        ----------
        df : pandas dataframe
            The sample compositions.
        P: pressure in GPa
            
        Returns
        -------
        temperature : float
            The calculated temperature(s) in oC.
        """
        compute_components_compound(self.df)  
        T = 1313.67
        T += 8.423 * ((P/0.1)-1)
        T -= 149.92 * (1. - Mg_num(self.df))
        T += 55.02 * NaK_num(self.df)
        T -= 59.69 * self.df['P2O5_primary_wt_dry']
        return {'P': P, 'T': T, 'P_err': P * self.T_err/T, 'T_err': self.T_err} 
    
    def compute_pressure_temperature(self):
        """
        Compute equilibration pressure and temperature for a given sample.
        
        Assumes garnet is stable mantle phase.
        Line (1) and (6) of Table (4), Grove et al (2013, Contrib Min Pet).
    
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
            Associated errors for pressure and temperature.
        """      
        compute_components_compound(self.df)
        P = self.compute_pressure()
        T = self.compute_temperature1(P)
        return {'P': P, 'T': T, 'P_err': self.P_err, 'T_err': self.T_err} 

class BK21_GNT:

    def __init__(self, df):
        self.df = df
        self.P_err = 1.26
        self.T_err = 18.72

    def compute_temperature(self, P):
        """
        Calculate equlibration temperature.
        
        Equation 1G-T of Table (1), Brown Krein et al (2021, JGR: Solid Earth).
        
        Parameters
        ----------
        df : pandas dataframe
            Dataframe containing the sample primary composition.
            Should contain only one row. To use with a multi-row dataframe use
            df.apply().
        P : float
        Pressure in GPa
            
        Returns
        -------
        out : dict
            The equilibration pressure and temperature result.
            Temperature is in degrees celcius.
            Pressure is in GPa.
            Associated errors for pressure and temperature.
        """                  
        T = compute_parameter_BK21(self.df, 0., P/0.1, [1136.36052,0,8.73915,184.88412,-19.48245,29.00187,-23.41730,-22.48205])
        return {'P': P, 'T': T, 'P_err': P * self.T_err/T, 'T_err': self.T_err} 
    
    def compute_pressure_temperature(self):
        """
        Same as TGK12_GNT but equations have been recalibrated with additional data.

        Equations 1G-P and 1G-T of Table (1), Brown Krein et al (2021, JGR: Solid Earth).
        
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
            Associated errors for pressure and temperature.
        """       
        compute_components_compound(self.df)
        MINS = compute_minerals(self.df)
        P = 0.1 * compute_parameter_BK21(self.df, MINS[0], 0., [-77.53423,139.87872,0,46.43963,28.72057,1.52179,3.25095,18.62211])
        T = compute_parameter_BK21(self.df, 0., P/0.1, [1136.36052,0,8.73915,184.88412,-19.48245,29.00187,-23.41730,-22.48205])
        return {'P': P, 'T': T, 'P_err': self.P_err, 'T_err': self.T_err} 
    
class BK21_SPL:

    def __init__(self, df):
        self.df = df
        self.P_err = 1.53
        self.T_err = 14.28

    def compute_temperature(self, P):
        """
        Calculate equlibration temperature.
        
        Line (17) of Table (1), Brown Krein et al (2021, JGR: Solid Earth).
        
        Parameters
        ----------
        df : pandas dataframe
            Dataframe containing the sample primary composition.
            Should contain only one row. To use with a multi-row dataframe use
            df.apply().
        P : float
        Pressure in GPa
            
        Returns
        -------
        out : dict
            The equilibration pressure and temperature result.
            Temperature is in degrees celcius.
            Pressure is in GPa.
            Associated errors for pressure and temperature.
        """                  
        T = compute_parameter_BK21(self.df, 0., P/0.1, [1049.24929,0,12.71243,63.47484,-3.32516,2.65802,-12.03128,117.75713])
        return {'P': P, 'T': T, 'P_err': P * self.T_err/T, 'T_err': self.T_err} 

    def compute_pressure_temperature(self):
        """
        Same as TGK12_SPL but equations have been recalibrated with additional data.

        Line (3 and 17) of Table (1), Brown Krein et al (2021, JGR: Solid Earth).
        
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
            Associated errors for pressure and temperature.
        """       
        compute_components_compound(self.df)
        MINS = compute_minerals(self.df)
        P = 0.1 * compute_parameter_BK21(self.df, MINS[0], 0., [-29.5,84.82,0,25.,24.67,2.79,-0.138,-1.848]) 
        T = compute_parameter_BK21(self.df, 0., P/0.1, [1049.24929,0,12.71243,63.47484,-3.32516,2.65802,-12.03128,117.75713])
        return {'P': P, 'T': T, 'P_err': self.P_err, 'T_err': self.T_err} 

class BK21_PLG:
   
    def __init__(self, df):
        self.df = df
        self.P_err = 0.79
        self.T_err = 11.5
        
    def compute_temperature(self, P):
        """
        Calculate equlibration temperature.
        
        Line (6) of Table (1), Brown Krein et al (2021, JGR: Solid Earth).
        
        Parameters
        ----------
        df : pandas dataframe
            Dataframe containing the sample primary composition.
            Should contain only one row. To use with a multi-row dataframe use
            df.apply().
        P : float
        Pressure in GPa
            
        Returns
        -------
        out : dict
            The equilibration pressure and temperature result.
            Temperature is in degrees celcius.
            Pressure is in GPa.
            Associated errors for pressure and temperature.
        """                  
        T = compute_parameter_BK21(self.df, 0., P/0.1, [1074.38633,0,11.86431,65.55420,-138.22714,20.55173,5.85532,79.01883])
        return {'P': P, 'T': T, 'P_err': P * self.T_err/T, 'T_err': self.T_err} 
    
    def compute_pressure_temperature(self):
        """
        Calculate equlibration conditions.
        
        Line (1 and 7) of Table (1), Brown Krein et al (2021, JGR: Solid Earth).
        
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
            Associated errors for pressure and temperature.
        """     
        compute_components_compound(self.df)
        MINS = compute_minerals(self.df)
        P = 0.1 * compute_parameter_BK21(self.df, MINS[0], 0., [-43.58532,136.94746,0,24.54202,37.54688,2.10967,-0.77250,1.63584])
        T = compute_parameter_BK21(self.df, 0., P/0.1, [1074.38633,0,11.86431,65.55420,-138.22714,20.55173,5.85532,79.01883])
        return {'P': P, 'T': T, 'P_err': self.P_err, 'T_err': self.T_err} 

class BK21:
    
    def __init__(self, df):
        self.df = df
        self.P_err = [0.79, 1.53, 1.26]
        self.T_err = [11.5, 14.28, 18.72]
    
    def compute_pressure(self):
        """
        Calculate equlibration pressure.

        Brown Krein et al (2021, JGR: Solid Earth).

        Parameters
        ----------
        df : pandas dataframe
            The sample compositions in wt%.
            
        Returns
        -------
        pressure : array
            The calculated pressure(s) in GPa
            in plg, spl and gnt stability fields.
        """
        MINS = compute_minerals(self.df)
        pressure = [0]*3
        # Compute Pressure Plag        
        pressure[0] = 0.1 * compute_parameter_BK21(self.df, MINS[0], 0., [-43.58532,136.94746,0,24.54202,37.54688,2.10967,-0.77250,1.63584]) 
        # Compute Pressure Spl
        pressure[1] = 0.1 * compute_parameter_BK21(self.df, MINS[0], 0., [-29.5,84.82,0,25.,24.67,2.79,-0.138,-1.848]) 
        # Compute Pressure Gnt 
        pressure[2] = 0.1 * compute_parameter_BK21(self.df, MINS[0], 0., [-77.53423,139.87872,0,46.43963,28.72057,1.52179,3.25095,18.62211])
  
        return pressure

    def compute_temperature1(self, P):
        """
        Calculate equlibration temperature.

        Brown Krein et al (2021, JGR: Solid Earth).

        Parameters
        ----------
        df : pandas dataframe
            The sample compositions in wt%.
        P : float
            Pressure in GPa
                
        Returns
        -------
        temperature : array
            The calculated pressure(s) in GPa 
            in plg, spl and gnt stability fields.
        """
        temperature = [0]*3
        # Compute Temperature Plag
        temperature[0] = compute_parameter_BK21(self.df, 0., P[0]/0.1, [1074.38633,0,11.86431,65.55420,-138.22714,20.55173,5.85532,79.01883])
        # Compute Temperature Spl
        temperature[1] = compute_parameter_BK21(self.df, 0., P[1]/0.1, [1049.24929,0,12.71243,63.47484,-3.32516,2.65802,-12.03128,117.75713])           
        # Compute Temperature Gnt
        temperature[2] = compute_parameter_BK21(self.df, 0., P[2]/0.1, [1136.36052,0,8.73915,184.88412,-19.48245,29.00187,-23.41730,-22.48205])           
                  
        return temperature
    
    def compute_temperature(self, P):
        """
        Calculate equlibration temperature assuming a known pressure.

        Brown Krein et al (2021, JGR: Solid Earth).

        Parameters
        ----------
        df : pandas dataframe
            The sample compositions in wt%.
        P : float
            Pressure in GPa
            
        Returns
        -------
        temperature : array
            The calculated pressure(s) in GPa 
            in plg, spl and gnt stability fields.
        """
        compute_components_compound(self.df)
        P = self.compute_pressure()
        RMSD = self.compute_mins(P)
        T = [0]*3
        # Compute Temperature Plag
        T[0] = compute_parameter_BK21(self.df, 0., P[0]/0.1, [1074.38633,0,11.86431,65.55420,-138.22714,20.55173,5.85532,79.01883])
        # Compute Temperature Spl
        T[1] = compute_parameter_BK21(self.df, 0., P[1]/0.1, [1049.24929,0,12.71243,63.47484,-3.32516,2.65802,-12.03128,117.75713])           
        # Compute Temperature Gnt
        T[2] = compute_parameter_BK21(self.df, 0., P[2]/0.1, [1136.36052,0,8.73915,184.88412,-19.48245,29.00187,-23.41730,-22.48205])           
                  
        return {'P': P[np.argmin(RMSD)], 'T': T[np.argmin(RMSD)],
                'P_err': P[np.argmin(RMSD)] * self.T_err[np.argmin(RMSD)]/T[np.argmin(RMSD)], 
                'T_err': self.T_err[np.argmin(RMSD)]}   

    def compute_mins(self, P):
        """
        Determine root mean squared differences between predicted
        sample mineralogy and mineralogies in plg, spl and gnt fields.

        Brown Krein et al (2021, JGR: Solid Earth).

        Parameters
        ----------
        df : pandas dataframe
            The sample compositions in wt%.
        P: float
            pressure in GPa
        
        Returns
        -------
        RMSD : array
            Difference between calculated and predicted
            mineralogic make-up for plg, spl and gnt fields.
        """
        MINS = compute_minerals(self.df)
        RMSD = [0]*3
        # Calculated minerals order: OL, CPX, PLG, QTZ
        # Mineral Residuals for Plagioclase
        Mplg = [0]*4
        Mplg[0] = compute_parameter_BK21(self.df, 0., P[0]/0.1, [0.412,0,0.005,-0.162,-0.363,-0.011,-0.001,-0.126])
        Mplg[0] -= MINS[0]
        Mplg[1] = compute_parameter_BK21(self.df, 0., P[0]/0.1, [-0.256,0,0,0.032,0.339,-0.004,-0.01,0.497])
        Mplg[1] -= MINS[1]        
        Mplg[2] = compute_parameter_BK21(self.df, 0., P[0]/0.1, [0.429,0,0.012,0.132,0.503,0.001,-0.073,-0.157])        
        Mplg[2] -= MINS[2]        
        Mplg[3] = compute_parameter_BK21(self.df, 0., P[0]/0.1, [0.423,0,-0.018,0.001,-0.487,-0.001,0.023,-0.225])       
        Mplg[3] -= MINS[3]
        
        RMSD[0] = np.sqrt((Mplg[0]**2. + Mplg[1]**2. + Mplg[2]**2. + Mplg[3]**2.)/4.)
        
        # Mineral Residuals for Spinel
        Mspl = [0]*4
        Mspl[0] = compute_parameter_BK21(self.df, 0., P[1]/0.1, [0.294,0,0.009,-0.215,-0.208,-0.027,-0.003,-0.046])
        Mspl[0] -= MINS[0]        
        Mspl[1] = compute_parameter_BK21(self.df, 0., P[1]/0.1, [-0.256,0,0,0.016,0.336,-0.001,-0.006,0.525])
        Mspl[1] -= MINS[1]        
        Mspl[2] = compute_parameter_BK21(self.df, 0., P[1]/0.1, [0.826,0,-0.004,0.116,0.392,-0.012,-0.072,-0.482]) 
        Mspl[2] -= MINS[2]        
        Mspl[3] = compute_parameter_BK21(self.df, 0., P[1]/0.1, [0.148,0,-0.005,0.085,-0.53,0.023,0.02,-0.107])  
        Mspl[3] -= MINS[3]    
        
        RMSD[1] = np.sqrt((Mspl[0]**2. + Mspl[1]**2. + Mspl[2]**2. + Mspl[3]**2.)/4.)
         
        # Mineral Residuals for Garnet
        Mgnt = [0]*4
        Mgnt[0] = compute_parameter_BK21(self.df, 0., P[2]/0.1, [0.531,0,0.007,-0.314,-0.185,-0.01,-0.023,-0.107])
        Mgnt[0] -= MINS[0]          
        Mgnt[1] = compute_parameter_BK21(self.df, 0., P[2]/0.1, [-0.192,0,-0.003,0.109,0.358,-0.005,-0.011,0.413])
        Mgnt[1] -= MINS[1]          
        Mgnt[2] = compute_parameter_BK21(self.df, 0., P[2]/0.1, [0.613,0,-0.006,0.237,0.365,-0.018,-0.042,-0.243]) 
        Mgnt[2] -= MINS[2]          
        Mgnt[3] = compute_parameter_BK21(self.df, 0., P[2]/0.1, [0.048,0,0.003,-0.026,-0.547,0.019,0.010,-0.069])  
        Mgnt[3] -= MINS[3]   

        RMSD[2] = np.sqrt((Mgnt[0]**2. + Mgnt[1]**2. + Mgnt[2]**2. + Mgnt[3]**2.)/4.)
        
        return RMSD      
    
    def compute_pressure_temperature(self):
        """
        Calculate equlibration conditions.
        
        Brown Krein et al (2021, JGR: Solid Earth).
        
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
            Associated errors for pressure and temperature.
        """    
        compute_components_compound(self.df)
        P = self.compute_pressure()
        T = self.compute_temperature1(P)
        RMSD = self.compute_mins(P)
        return {'P': P[np.argmin(RMSD)], 'T': T[np.argmin(RMSD)],
                'P_err': self.P_err[np.argmin(RMSD)], 'T_err': self.T_err[np.argmin(RMSD)]} 

class HA15:   

    def __init__(self, df):
        self.df = df
        self.P_err = np.nan
        self.T_err = 31.
    
    def compute_temperature(self, P):
        """
        Calculate temperate of olivine liquidus at 1 bar
        using Beattie (1993) Cointrib Min & Pet 115:103-111 
        Equations 10 and 12.

        Convert to desired pressure using Equation from 
        Herzberg & O'Hara, (2002) G3 9:9  
        
        Parameters
        ----------
        df : pandas dataframe
            Dataframe containing the sample primary composition.
            Should contain only one row. To use with a multi-row dataframe use
            df.apply().
        P : float
            assumed pressure in GPa.
            
        Returns
        -------
        out : dict
            The equilibration pressure and temperature result.
            Temperature is in degrees celcius.
            Pressure is in GPa.
            Associated errors for pressure and temperature.
        """    
        compute_components_cation(self.df)
        DMg = compute_DMg(self.df)   
        
        #####
        Temp_Denom = 52.05/8.3143
        Temp_Denom += 2.*np.log(DMg)
        Temp_Denom += 2.*np.log(1.5*(self.df["FeO_primary_mol_dry"]+self.df["MnO_primary_mol_dry"]+self.df["MgO_primary_mol_dry"]+self.df["CaO_primary_mol_dry"]+self.df["NiO_primary_mol_dry"]))
        Temp_Denom += 2.*np.log(3.*self.df["SiO2_primary_mol_dry"])
        Temp_Denom -= 3.5*np.log(1-self.df["Al2O3_primary_mol_dry"]) + 7*np.log(1-self.df["TiO2_primary_mol_dry"])
        Temp_1bar = (113100/8.3143 + P*0.00000411/8.3143)
        Temp_1bar /= Temp_Denom
        Temp_1bar -= 273.15
        #####
        T = Temp_1bar + 54.*P - 2.*P**2.
        
        return {'P': P, 'T': T, 'P_err': P * self.T_err/T, 'T_err': self.T_err} 

class P08:
    
    def __init__(self, df):
        self.df = df
        self.P_err = 0.29
        self.T_err = 52.
    
    def compute_pressure(self, T):
        """
        Calculate temperate and pressure of melt
        equilibration by simultaneously solving
        Equation 4 of Putirka et al., (2007, Chemical Geology)
        and Equation 42 of Putirka (2008).
    
        Parameters
        ----------
        df : pandas dataframe
            Dataframe containing the sample primary composition.
            Should contain only one row. To use with a multi-row dataframe use
            df.apply().
        T : float
            Temperature in oC
            
        Returns
        -------
        out : dict
            The equilibration pressure and temperature result.
            Temperature is in degrees celcius.
            Pressure is in GPa.
            Associated errors for pressure and temperature.
        """    
        
        compute_components_cation(self.df)
        aSiO2 = (3*self.df["SiO2_primary_mol_dry"])**-2*(1-self.df["Al2O3_primary_mol_dry"])**(7/2)*(1-self.df["TiO2_primary_mol_dry"])**(7/2)
        P = (231.5+0.186*T+0.1244*T*(np.log(aSiO2))-528.*(aSiO2)**0.5+103.3*self.df["TiO2_primary_mol_dry"]+69.9*(self.df["Na2O_primary_mol_dry"]+self.df["K2O_primary_mol_dry"])+77.3*(self.df["Al2O3_primary_mol_dry"]/(self.df["Al2O3_primary_mol_dry"]+self.df["SiO2_primary_mol_dry"])))/10.
            
        return {'P': P, 'T': T, 'P_err': self.P_err, 'T_err': T * self.P_err/P} 

    
    def compute_pressure_temperature(self):
        """
        Calculate temperate and pressure of melt
        equilibration by simultaneously solving
        Equation 4 of Putirka et al., (2007, Chemical Geology)
        and Equation 42 of Putirka (2008).
    
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
            Associated errors for pressure and temperature.
        """    
        
        compute_components_cation(self.df)
        DMg = compute_DMg(self.df)   
        aSiO2 = (3*self.df["SiO2_primary_mol_dry"])**-2*(1-self.df["Al2O3_primary_mol_dry"])**(7/2)*(1-self.df["TiO2_primary_mol_dry"])**(7/2)

        P, T = sym.symbols('P, T')
        eq1 = sym.Eq(1 / ((np.log(DMg) + 2.158 - 5.115*10**-2*(self.df["Na2O_primary_wt"]+self.df["K2O_primary_wt"]) + 6.213*10**-2*self.df["H2O_primary_wt"]) / (55.09*P + 4430.)), T)
        eq2 = sym.Eq((231.5+0.186*T+0.1244*T*(np.log(aSiO2))-528.*(aSiO2)**0.5+103.3*self.df["TiO2_primary_mol_dry"]+69.9*(self.df["Na2O_primary_mol_dry"]+self.df["K2O_primary_mol_dry"])+77.3*(self.df["Al2O3_primary_mol_dry"]/(self.df["Al2O3_primary_mol_dry"]+self.df["SiO2_primary_mol_dry"])))/10., P)
        result = sym.solve([eq1,eq2],(P,T))
        T = result[T]
        P = result[P]
            
        return {'P': P, 'T': T, 'P_err': self.P_err, 'T_err': self.T_err} 

class SD20:

    def __init__(self, df):
        self.df = df
        self.P_err = 0.5
        self.T_err = 49.   
    
    def compute_temperature(self, P):
        """
        Sun & Dasgupta, (2020), EPSL.
        Temperature calculated using equation 6.
    
        Parameters
        ----------
        df : pandas dataframe
            Dataframe containing the sample primary composition.
            Should contain only one row. To use with a multi-row dataframe use
            df.apply().
        P : float
            pressure in GPa
            
        Returns
        -------
        out : dict
            The equilibration pressure and temperature result.
            Temperature is in degrees celcius.
            Pressure is in GPa.
            Associated errors for pressure and temperature.
        """          

        compute_components_per_oxygen(self.df)
        
        # Phi and Theta defined in Eq 3,
        # Omega defined in Eq 5
        phi = 2412. + 9.2 * (self.df['CO2_primary_wt']/self.df['SiO2_primary_mol_dry'])
        phi += 11.19 * (self.df['H2O_primary_wt']/self.df['SiO2_primary_mol_dry'])
        
        theta = 4.7 + 3.91 * np.log(self.df['SiO2_primary_mol_dry']) 
        theta += 1.19*((self.df['K2O_primary_mol_dry']+self.df['TiO2_primary_mol_dry'])/self.df['MgO_primary_mol_dry'])
        
        omega = 2.59 + 3.5*(self.df['CaO_primary_mol_dry'] - (2.*self.df['K2O_primary_mol_dry'])) 
        omega += 4.85*self.df['TiO2_primary_mol_dry'] 
        omega += 1.4*(self.df['MgO_primary_mol_dry']/(self.df['MgO_primary_mol_dry']+self.df['FeO_primary_mol_dry'])) 
        omega += 0.5 * self.df['MgO_primary_mol_dry'] * sym.sqrt(self.df['CO2_primary_wt']) 
        omega += 0.057 * self.df['H2O_primary_wt']      
        
        T = (10.**4. / (omega - 0.34*np.sqrt(P) - 1.26*np.log(self.df['MgO_primary_mol_dry']))) - 273.15
        
        return {'P': P, 'T': float(T), 'P_err': float(P * self.T_err/T), 'T_err': self.T_err} 
    
    def compute_pressure_temperature(self):
        """
        Sun & Dasgupta, (2020), EPSL.
        Pressure and temperature calculated using
        equations 4 and 6, repsectively.
        Pressure and temperature are functions of each other
        and must be calculated simultaneously.
    
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
            Associated errors for pressure and temperature.
        """          

        compute_components_per_oxygen(self.df)
        
        # Phi and Theta defined in Eq 3,
        # Omega defined in Eq 5
        phi = 2412. + 9.2 * (self.df['CO2_primary_wt']/self.df['SiO2_primary_mol_dry'])
        phi += 11.19 * (self.df['H2O_primary_wt']/self.df['SiO2_primary_mol_dry'])
        
        theta = 4.7 + 3.91 * np.log(self.df['SiO2_primary_mol_dry']) 
        theta += 1.19*((self.df['K2O_primary_mol_dry']+self.df['TiO2_primary_mol_dry'])/self.df['MgO_primary_mol_dry'])
        
        omega = 2.59 + 3.5*(self.df['CaO_primary_mol_dry'] - (2.*self.df['K2O_primary_mol_dry'])) 
        omega += 4.85*self.df['TiO2_primary_mol_dry'] 
        omega += 1.4*(self.df['MgO_primary_mol_dry']/(self.df['MgO_primary_mol_dry']+self.df['FeO_primary_mol_dry'])) 
        omega += 0.5 * self.df['MgO_primary_mol_dry'] * sym.sqrt(self.df['CO2_primary_wt']) 
        omega += 0.057 * self.df['H2O_primary_wt']        

        # Pressure calculated using Eq 4.
        # Temperature calculated using Eq 6.
        # Must be solved simultaneously.       
        P, T = sym.symbols('P, T')
        eq1 = sym.Eq(13. - sym.sqrt((T/45.9) * (np.log(self.df['Al2O3_primary_mol_dry']/self.df['MgO_primary_mol_dry']) - theta) - (phi/45.9) + 169.), P)
        eq2 = sym.Eq(10.**4. / (omega - 0.34*sym.sqrt(P) - 1.26*np.log(self.df['MgO_primary_mol_dry'])), T)
        result = sym.solve([eq1,eq2],(P,T))
        if len(result) == 0:
            T = np.nan
            P = np.nan
        else:
            T = result[0][1] - 273.15
            P = result[0][0]
        return {'P': float(P), 'T': float(T), 'P_err': self.P_err, 'T_err': self.T_err} 

class B93:

    def __init__(self, df):
        self.df = df
        self.T_err = 20.
    
    def compute_temperature(self, P):
        """
        Calculate temperate of olivine liquidus at P GPa
        using Equation 10 Beattie (1993, Contrib. to Min. and Pet.)
        
        Parameters
        ----------
        df : pandas dataframe
            Dataframe containing the sample primary composition.
            Should contain only one row. To use with a multi-row dataframe use
            df.apply().
        P : float
            pressure in GPa
            
        Returns
        -------
        out : dict
            The equilibration temperature result.
            Temperature is in degrees celcius.
            Assumed pressure is in GPa.
            Associated errors for pressure and temperature.
        """    
        compute_components_cation(self.df)
        DMg = compute_DMg(self.df)       
        
        # Eq 20 of Putirka 2008 Rev in Min. & Geochem. 69:1
        # DMg = (0.666-(-0.049*self.df["MnO_primary_mol"]+0.027*self.df["FeO_primary_mol"]))/(self.df["MgO_primary_mol"]+0.259*self.df["MnO_primary_mol_dry"]+0.299*self.df["FeO_primary_mol"])
        C_NM = self.df["MgO_primary_mol"] + self.df["FeO_primary_mol"] 
        C_NM += self.df["CaO_primary_mol"] + self.df["MnO_primary_mol"]
        NF = 7/2*np.log(1.-self.df["Al2O3_primary_mol"]) 
        NF += 7.*np.log(1.-self.df["TiO2_primary_mol"])

        # Eq 19 of Putirka 2008 Rev in Min. & Geochem. 69:1
        T = ((113100./8.3144) + ((P*10.**9.-10.**-5.)*((4.11*10**-6)/8.3144)))
        T /= ((52.05/8.3144) + 2.*np.log(DMg) + 2.*np.log(1.5*C_NM) + 2.*np.log(3.*self.df["SiO2_primary_mol"]) - NF)
        T -= 273.15
  
        return {'P': P, 'T': T, 'P_err': P * self.T_err/T, 'T_err': self.T_err} 

class P07_2:

    def __init__(self, df):
        self.df = df
        self.T_err = 52.
    
    def compute_temperature(self, P):
        """
        Calculate temperate of olivine liquidus at P GPa
        using Equation 2 of Putirka et al., (2007, Chemical Geology)
        
        Parameters
        Parameters
        ----------
        df : pandas dataframe
            Dataframe containing the sample primary composition.
            Should contain only one row. To use with a multi-row dataframe use
            df.apply().
        P : float
            pressure in GPa
            
        Returns
        -------
        out : dict
            The equilibration temperature result.
            Temperature is in degrees celcius.
            Assumed pressure is in GPa.
            Associated errors for pressure and temperature.
        """
        compute_components_cation(self.df)
        DMg = compute_DMg(self.df)   

        # Eq 21 of Putirka 2008 Rev in Min. & Geochem. 69:1
        T = 1 / ((np.log(DMg) + 2.158 - 5.115*10**-2*(self.df["Na2O_primary_wt"]+self.df["K2O_primary_wt"]) + 6.213*10**-2*self.df["H2O_primary_wt"]) / (55.09*P + 4430))
  
        return {'P': P, 'T': T, 'P_err': P * self.T_err/T, 'T_err': self.T_err} 

class P07_4:

    def __init__(self, df):
        self.df = df
        self.T_err = 52.
    
    def compute_temperature(self, P):
        """
        Calculate temperate of olivine liquidus at P GPa
        using Equation 4 of Putirka et al., (2007, Chemical Geology)
        
        Parameters
        ----------
        df : pandas dataframe
            Dataframe containing the sample primary composition.
            Should contain only one row. To use with a multi-row dataframe use
            df.apply().
        P : float
            pressure in GPa
            
        Returns
        -------
        out : dict
            The equilibration temperature result.
            Temperature is in degrees celcius.
            Assumed pressure is in GPa.
            Associated errors for pressure and temperature.    
        """
        
        compute_components_cation(self.df)
        DMg = compute_DMg(self.df)   

        C_NM = self.df["MgO_primary_mol"] + self.df["FeO_primary_mol"] 
        C_NM += self.df["CaO_primary_mol"] + self.df["MnO_primary_mol"]
        
        NF = 7/2*np.log(1.-self.df["Al2O3_primary_mol"]) 
        NF += 7.*np.log(1.-self.df["TiO2_primary_mol"])

        # Eq 22 of Putirka 2008 Rev in Min. & Geochem. 69:1        
        T = (15294.6+1318.8*P+2.4834*P**2)/(8.048+2.8352*np.log(DMg)+2.097*np.log(1.5*C_NM)+2.575*np.log(3*self.df["SiO2_primary_mol"])-1.41*NF+0.222*self.df["H2O_primary_wt"]+0.5*P)
            
        return {'P': P, 'T': T, 'P_err': P * self.T_err/T, 'T_err': self.T_err} 

def compute_sample_pressure_temperature(df, method="PF16", min_SiO2=0.):
    """
    Calculate temperate and pressure of melt equilibration.
    
    Parameters
    ----------
    df : pandas dataframe.
        Dataframe containing the sample primary composition.
        Should contain only one row. To use with a multi-row dataframe use
        df.apply().    
    
    method: string
        choice of thermobaric method
        PF16 = Plank & Forsyth (2016), G3
        L09 = Lee et al., (2009), G3
        TGK12_PLG = Till et al., (2012) in plagioclase stability field.
        TGK12_SPL = Till et al., (2012) in spinel stability field.
        G13 = Grove et al., (2013) in garnet stability field.
        BK21_PLG = Brown Krein et al., (2021) in plagioclase stability field.
        BK21_SPL = Brown Krein et al., (2021) in spinel stability field.
        BK21_GNT = Brown Krein et al., (2021) in garnet stability field.
        BK21 = Brown Krein et al., (2021)
        P07_P08 = Putirka et al., (2007) Eq 4. and Putirka (2008) combined.
        SD20 = Sun & Dasgupta (2020), EPSL low SiO2 thermobarometer.
            
    Returns
    -------
    PT  : dict
            The calculated temperature(s) in oC,
            pressures in GPa, and associated errors.  
    """     
    
    if df['SiO2_primary_wt'] < min_SiO2:
        out = {'P': np.nan, 'P_err': np.nan, 'T': np.nan, 'T_err': np.nan}
    elif inspect.isclass(method):
        out = method(df).compute_pressure_temperature()
    elif method == "PF16":
        out = PF16(df).compute_pressure_temperature()
    elif method == "L09":
        out = L09(df).compute_pressure_temperature()
    elif method == "TGK12_PLG":
        out = TGK12_PLG(df).compute_pressure_temperature()
    elif method == "TGK12_SPL":
        out = TGK12_SPL(df).compute_pressure_temperature()
    elif method == "G13":
        out = G13(df).compute_pressure_temperature()
    elif method == "BK21_PLG":
        out = BK21_PLG(df).compute_pressure_temperature()
    elif method == "BK21_SPL":
        out = BK21_SPL(df).compute_pressure_temperature()
    elif method == "BK21_GNT":
        out = BK21_GNT(df).compute_pressure_temperature()
    elif method == "BK21":
        out = BK21(df).compute_pressure_temperature()        
    elif method == "P08":
        out = P08(df).compute_pressure_temperature()        
    elif method == "SD20":
        out = SD20(df).compute_pressure_temperature()        
    else:
        out = thermobar(df).compute_pressure_temperature()
        
    return out


def compute_sample_temperature(df, method="HA15", P=1., min_SiO2=0.):
    """
    Calculate temperate of melt equilibration.
    
    Parameters
    ----------
    df : pandas dataframe.
        Dataframe containing the sample primary composition.
        Should contain only one row. To use with a multi-row dataframe use
        df.apply().  
                                    
    method: sring
        choice of thermobaric method
        HA15 = Herzberg & Asimow (2015)
        P07_2 = Putirka et al., (2007) equation 2.
        P07_4 = Putirka et al., (2007) equation 4.
        TGK12_PLG = Till et al., (2012) in plagioclase stability field.
        TGK12_SPL = Till et al., (2012) in spinel stability field.        
        G13 = Grove et al., (2013) in garnet stability field.  
        BK21_PLG = Brown Krein et al., (2021) in plagioclase stability field.
        BK21_SPL = Brown Krein et al., (2021) in spinel stability field.
        BK21_GNT = Brown Krein et al., (2021) in garnet stability field.
        BK21 = Brown Krein et al., (2021).
        SD20 = Sun & Dasgupta (2020), EPSL low SiO2 thermobarometer.
        
    P  : float
        pressure in GPa  
        
    Returns
    -------
    PT  : pandas dataframe
        The calculated temperature(s) in oC,
        assumed pressures in GPa and associated errors.
    """     
    
    if df['SiO2_primary_wt'] < min_SiO2:
        out = {'P': np.nan, 'P_err': np.nan, 'T': np.nan, 'T_err': np.nan}
    elif inspect.isclass(method):
        out = method(df).compute_temperature(P)
    elif method == "HA15":
        out = HA15(df).compute_temperature(P) 
    elif method == "P07_2":
        out = P07_2(df).compute_temperature(P)          
    elif method == "P07_4":
        out = P07_4(df).compute_temperature(P)
    elif method == "B93":
        out = B93(df).compute_temperature(P)        
    elif method == "TGK12_PLG":
        out = TGK12_PLG(df).compute_temperature(P)        
    elif method == "TGK12_SPL":
        out = TGK12_SPL(df).compute_temperature(P)  
    elif method == "G13":
        out = G13(df).compute_temperature(P)      
    elif method == "BK21_PLG":
        out = BK21_PLG(df).compute_temperature(P)        
    elif method == "BK21_SPL":
        out = BK21_SPL(df).compute_temperature(P)  
    elif method == "BK21_GNT":
        out = BK21_GNT(df).compute_temperature(P)  
    elif method == "BK21":
        out = BK21(df).compute_temperature(P)        
    elif method == "SD20":
        out = SD20(df).compute_temperature(P)            
    else:
        out = thermobar(df).compute_temperature(P)
        
    return out

def compute_sample_pressure(df, method="PF16", T=1300., min_SiO2=0.):
    """
    Calculate temperate of melt equilibration
    
    Parameters
    ----------
    df : pandas dataframe.
        Dataframe containing the sample primary composition.
        Should contain only one row. To use with a multi-row dataframe use
        df.apply().  
                                    
    method: sring
        choice of thermobaric method
        PF16 = Plank & Forsyth (2016), G3
        L09 = Lee et al., (2009), G3
        P08 = Putirka (2008) Eq 42.
        
    T  : float
        Temperature in oC
        
    Returns
    -------
    PT  : pandas dataframe
        The calculated temperature(s) in oC,
        assumed pressures in GPa and associated errors.
    """     
    
    if df['SiO2_primary_wt'] < min_SiO2:
        out = {'P': np.nan, 'P_err': np.nan, 'T': np.nan, 'T_err': np.nan}
    elif inspect.isclass(method):
        out = method(df).compute_pressure(T)
    elif method == "PF16":
        out = PF16(df).compute_pressure(T) 
    elif method == "L09":
        out = L09(df).compute_pressure(T)          
    elif method == "P08":
        out = P08(df).compute_pressure(T)        
    else:
        out = thermobar(df).compute_temperature(T)
        
    return out