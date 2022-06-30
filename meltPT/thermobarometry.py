import warnings

import numpy as np
import pandas as pd
import sympy as sym

def Mg_num(df):
    """
    Calculate Mg#
        
    Parameters
    ----------
    df : pandas dataframe
    The sample compositions in anhydrous mol%.
                
    Returns
    -------
    Mg#
    """   
    return df['MgO_primary_mol_dry'] / (df['MgO_primary_mol_dry'] + df['FeO_primary_mol_dry'])

def NaK_num(df):
    """
    Calculate NaK#
        
    Parameters
    ----------
    df : pandas dataframe
    The sample compositions in anhydrous mol%.
                
    Returns
    -------
    NaK#
    """            
    return (
        (df['Na2O_primary_wt_dry'] + df['K2O_primary_wt_dry']) / 
        (df['Na2O_primary_wt_dry'] + df['K2O_primary_wt_dry'] + df['CaO_primary_wt_dry'])
        )

def normalize_v2(df, list_parameters, input_suffix, output_suffix, end_sum):
    """
    Normalize list of parameters to sum to 100
    
    Parameters
    ----------
    df : pandas datframe.
    list_parameters : list of columns in dataframe.
    input_suffix : original suffix for columns to use.
    output_suffix : add suffix to make a new column name.
    end_sum : amount to sum to.

    Returns
    -------
    df : pandas datframe.
    """
    SUM = 0.
    for i in range(0, len(list_parameters), 1):
        parameter = list_parameters[i]
        SUM = df[parameter + input_suffix] + SUM

    for i in range(0, len(list_parameters), 1):
        parameter = list_parameters[i]        
        df[parameter + output_suffix] = end_sum * df[parameter + input_suffix] / SUM
        
    return df

def compute_olivine(df):
    """
    Calculate mineral components using equations
    listed in Grove (1982) Contrib Mineral Petrol 80:160-182
        
    Parameters
    ----------
    df : pandas dataframe
    The sample compositions in mol% anhydrous.
                
    Returns
    -------
    OL: predicted proportion of olivine
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
    ol = feo + mgo + 5/3*p2o5
    ol += 0.5*(al2o3 - k2o - na2o)
    ol -= cao + tio2 + cr2o3
    ol *= 0.5
    ol /= Total*2
    ####
    qz = sio2
    qz -= 0.5*(feo+mgo)
    qz -= 1.5*cao
    qz -= 0.25*al2o3   
    qz -= 2.75*(k2o+na2o)
    qz += 0.5*(cr2o3+tio2)
    qz += 2.5*p2o5
    qz /= Total    
    ####
    cpx = cao
    cpx -= 0.5*al2o3
    cpx += 0.5*(k2o+na2o)
    cpx -= 5/3*p2o5 
    cpx /= Total    
    ####    
    plg = 0.5*(al2o3+na2o-k2o)
    plg /= Total*4
    ####    
    opx = k2o/Total*4    
    #### 
    spl = (cr2o3+tio2)/Total*1.5
    ####   
    ap = p2o5/3/Total*6    
    ####     
    Min_Tot = ol + qz + plg + cpx + opx + spl + ap
    ####       
    OL = ol / Min_Tot
    OL /= (ol+qz+plg+cpx) / Min_Tot 
    return OL

def compute_components_species(df):
    """
    Calculate species proportions using equations
    listed in Appendix A of Lee et al (2009) G3
        
    Parameters
    ----------
    df : pandas dataframe
    The sample compositions in wt%.
                
    Returns
    -------
    df with additional columns of oxide mol%, species proportions and
    anhydrous species proportions.
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
                
    Returns
    -------
    df with additional columns of oxide mol%.
    """ 
        
    MAJOR_OXIDES = ['SiO2','Al2O3','FeO','MgO','CaO','Na2O','K2O','TiO2',
        'MnO','Cr2O3','P2O5', 'NiO']
    OXIDE_WEIGHT = [60.085,101.962,71.846,40.311,56.079,61.979,
        94.203,79.899,70.937,151.99,141.945,74.71]

    df['FeO_primary_wt'] = df['FeO_primary_wt'] + 0.8998 * df['Fe2O3_primary_wt']
    df['Fe2O3_primary_wt'] = 0.

    normalize_v2(df, MAJOR_OXIDES, '_primary_wt', '_primary_wt_dry', 100.)

    for i in range(0, len(MAJOR_OXIDES), 1):
        OXIDE = MAJOR_OXIDES[i]
        df[OXIDE + '_primary_mol_dry'] = df[OXIDE + '_primary_wt_dry'] / OXIDE_WEIGHT[i]
        
    return

def compute_components_cation(df):
    """
    Calculate proportions of cations in mol% assuming the sample is anhydrous.
        
    Parameters
    ----------
    df : pandas dataframe
    The sample compositions in wt%.
                
    Returns
    -------
    df with additional columns of oxide mol%.
    """
        
    MAJOR_OXIDES = ['SiO2','Al2O3','FeO','Fe2O3','MgO','CaO','Na2O','K2O',
        'TiO2','MnO','Cr2O3', 'P2O5', 'NiO']
    CAT_WEIGHT = [60.085,101.962/2.,71.846,159.96/2.,40.311,56.079,61.979/2.,
        94.203/2.,79.899,70.937,151.99/2.,141.945/2.,74.71]

    normalize_v2(df, MAJOR_OXIDES, '_primary_wt', '_primary_wt_dry', 100.) 

    for i in range(0, len(MAJOR_OXIDES), 1):
        OXIDE = MAJOR_OXIDES[i]
        df[OXIDE + '_primary_mol_dry'] = df[OXIDE + '_primary_wt_dry'] / CAT_WEIGHT[i]
    
    normalize_v2(df, MAJOR_OXIDES, '_primary_mol_dry', '_primary_mol_dry', 100.)     
        
    return

class PF16:
    
    def __init__(self, df):
        self.df = df

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

    def compute_pressure(self, T):
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
            np.log(self.df['Si4O8_dry']) -
            4.045 +
            0.0114*self.df['Fe4Si2O8_dry'] +
            0.00052*self.df['Ca4Si2O8_dry']**2. +
            0.0024*self.df['Mg4Si2O8_dry']) / (-336.3/T - 0.0007*np.sqrt(T)
            )
        return pressure

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
        """
        compute_components_species(self.df)
        water_correction = self.compute_water_correction()
        T = self.compute_temperature() - water_correction
        P = self.compute_pressure(T)
        if P > 2.:
            T -= self.compute_CO2_correction()
            P = self.compute_pressure(T)
        out = {'P': P, 'T': T - 273.15}
        return out

class L09:
        
    def __init__(self, df):
        self.df = df
        
    def compute_pressure(self, T):
        """
        Calculate equlibration temperature.

        Equation (2), Lee et al (2009, G-cubed).
        
        Parameters
        ----------
        df : pandas dataframe
            The sample compositions in species mol%.
        T : temperature in oC
            
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
    
    def compute_temperature(self):
        """
        Calculate equlibration temperature.

        Equation (3), Lee et al (2009, G-cubed).
        
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
            916.45 + 
            (13.68 * self.df['Mg4Si2O8']) + 
            (4580. / self.df['Si4O8']) -
            (0.509 * self.df['H16O8'] * self.df['Mg4Si2O8'])
            )
        return temperature
        
    def compute_pressure_temperature(self):
        compute_components_species(self.df)
        T = self.compute_temperature()
        P = self.compute_pressure(T)
        return {'P': P, 'T': T}

class TGK12_SPL:
    
    def __init__(self, df):
        self.df = df
    
    def compute_pressure(self):
        """
        Calculate equlibration pressure assuming
        spinel is stable.

        Line (8) of Table (5), Till et al (2012, JGR).
        
        Parameters
        ----------
        df : pandas dataframe
            The sample compositions.
            
        Returns
        -------
        temperature : pandas dataframe
            The calculated pressure(s).
        """
        pressure = -0.862
        pressure += 9.471 * compute_olivine(self.df)
        pressure -= 2.383 * (1. - Mg_num(self.df)) 
        pressure += 2.922 * NaK_num(self.df)
        pressure += 0.218 * self.df['TiO2_primary_wt_dry'] 
        pressure -= 0.146 * self.df['K2O_primary_wt_dry']
        return pressure

    def compute_temperature(self, P):
        """
        Calculate equlibration temperature assuming
        spinel is stable.

        Line (7) of Table (5), Till et al (2012, JGR).
        
        Parameters
        ----------
        df : pandas dataframe
            The sample compositions.
        P : equilibration pressure
            
        Returns
        -------
        temperature : pandas dataframe
            The calculated temperature(s).
        """
        temperature = 1212.
        temperature += 119.9 * P
        temperature -= 97.33 * (1. - Mg_num(self.df)) 
        temperature -= 87.76 * NaK_num(self.df)
        temperature += 3.44 * self.df['TiO2_primary_wt_dry'] 
        temperature -= 4.58 * self.df['K2O_primary_wt_dry']
        return temperature
    
    def compute_pressure_temperature(self):
        compute_components_compound(self.df)       
        P = self.compute_pressure()
        T = self.compute_temperature(P)
        return {'P': P, 'T': T}   
    
class TGK12_PLG:
    """
    Same as TGK12_SPL but assuming plagioclase is stable.

    Line (1 and 2) of Table (5), Till et al (2012, JGR).
        
    Parameters
    ----------
    df : pandas dataframe
        The sample compositions.
    P : equilibration pressure
            
    Returns
    -------
    PT : pandas dataframe
       The calculated pressure(s) and temperature(s).
    """    
    def __init__(self, df):
        self.df = df
    
    def compute_pressure(self):
        pressure = - 1.64
        pressure += 12.94 * compute_olivine(self.df)
        pressure -= 2.363 * (1. - Mg_num(self.df)) 
        pressure += 3.51 * NaK_num(self.df)
        pressure += 0.152 * self.df['TiO2_primary_wt_dry'] 
        pressure -= 0.176 * self.df['K2O_primary_wt_dry']
        return pressure

    def compute_temperature(self, P):
        temperature = 1216.
        temperature += 104.4 * P
        temperature -= 72.83 * (1. - Mg_num(self.df)) 
        temperature -= 194.9 * NaK_num(self.df)
        temperature += 24.08 * self.df['TiO2_primary_wt_dry'] 
        temperature -= 1.55 * self.df['K2O_primary_wt_dry']
        return temperature
    
    def compute_pressure_temperature(self):
        compute_components_compound(self.df)
        P = self.compute_pressure()
        T = self.compute_temperature(P)
        return {'P': P, 'T': T}  

class G13:
    """
    Same as TGK12_SPL but assuming garnet is stable.

    Line (1 and 6) of Table (4), Grove et al (2013, Contrib Min Pet).
        
    Parameters
    ----------
    df : pandas dataframe
        The sample compositions.
    P : equilibration pressure
            
    Returns
    -------
    PT : pandas dataframe
       The calculated pressure(s) and temperature(s).
    """ 
    def __init__(self, df):
        self.df = df   

    def compute_pressure(self):
        pressure = - 1.73
        pressure += 14.64 * compute_olivine(self.df)
        pressure -= 3.84 * (1. - Mg_num(self.df)) 
        pressure += 1.96 * NaK_num(self.df)
        pressure += 0.481 * self.df['P2O5_primary_wt_dry'] 
        return pressure

    def compute_temperature(self, P):
        temperature = 1313.67
        temperature += 8.423 * ((P/0.1)-1)
        temperature -= 149.92 * (1. - Mg_num(self.df)) 
        temperature += 55.02 * NaK_num(self.df)
        temperature -= 59.69 * self.df['P2O5_primary_wt_dry'] 
        return temperature
    
    def compute_pressure_temperature(self):
        compute_components_compound(self.df)
        P = self.compute_pressure()
        T = self.compute_temperature(P)
        return {'P': P, 'T': T} 

class K21_GNT:
    """
    Same as TGK12_GNT but equations have been recalibrated with additional data.

    Line (5 and 27) of Table (1), Krein et al (2021, JGR: Solid Earth).
        
    Parameters
    ----------
    df : pandas dataframe
        The sample compositions.
    P : equilibration pressure
            
    Returns
    -------
    PT : pandas dataframe
       The calculated pressure(s) and temperature(s).
    """ 
    def __init__(self, df):
        self.df = df
        
    def compute_pressure(self):
        pressure = - 77.53
        pressure += 139.9 * compute_olivine(self.df)
        pressure += 46.44 * Mg_num(self.df)
        pressure += 28.72 * NaK_num(self.df)
        pressure += 1.52 * self.df['TiO2_primary_wt_dry'] 
        pressure -= 3.25 * self.df['K2O_primary_wt_dry']         
        pressure += 18.62 * (self.df['CaO_primary_wt_dry']/self.df['Al2O3_primary_wt_dry']) 
        pressure *= 0.1              
        return pressure

    def compute_temperature(self, P):
        temperature = 1136.
        temperature += 8.739 * (P/0.1)
        temperature += 184.9 * Mg_num(self.df)
        temperature -= 19.48 * NaK_num(self.df)
        temperature += 29 * self.df['TiO2_primary_wt_dry'] 
        temperature -= 23.42 * self.df['K2O_primary_wt_dry']        
        temperature -= 22.48 * (self.df['CaO_primary_wt_dry'] / self.df['Al2O3_primary_wt_dry'])    
        return temperature
    
    def compute_pressure_temperature(self):
        compute_components_compound(self.df)
        P = self.compute_pressure()
        T = self.compute_temperature(P)
        return {'P': P, 'T': T} 
    
class K21_SPL:
    """
    Same as TGK12_SPL but equations have been recalibrated with additional data.

    Line (3 and 17) of Table (1), Krein et al (2021, JGR: Solid Earth).
        
    Parameters
    ----------
    df : pandas dataframe
        The sample compositions.
    P : equilibration pressure
            
    Returns
    -------
    PT : pandas dataframe
       The calculated pressure(s) and temperature(s).
    """ 
    def __init__(self, df):
        self.df = df

    def compute_pressure(self):
        pressure = - 29.5
        pressure += 84.82 * compute_olivine(self.df)
        pressure += 25 * Mg_num(self.df)
        pressure += 24.67 * NaK_num(self.df)
        pressure += 2.79 * self.df['TiO2_primary_wt_dry'] 
        pressure -= 0.138 * self.df['K2O_primary_wt_dry']         
        pressure -= 1.848 * (self.df['CaO_primary_wt_dry']/self.df['Al2O3_primary_wt_dry']) 
        pressure *= 0.1              
        return pressure

    def compute_temperature(self, P):
        temperature = 1049.
        temperature += 12.71 * (P/0.1)
        temperature += 63.47 * Mg_num(self.df)
        temperature -= 3.325 * NaK_num(self.df)
        temperature += 2.658 * self.df['TiO2_primary_wt_dry'] 
        temperature -= 12.03 * self.df['K2O_primary_wt_dry']        
        temperature += 117.8 * (self.df['CaO_primary_wt_dry'] / self.df['Al2O3_primary_wt_dry'])    
        return temperature
    
    def compute_pressure_temperature(self):
        compute_components_compound(self.df)
        P = self.compute_pressure()
        T = self.compute_temperature(P)
        return {'P': P, 'T': T} 

class K21_PLG:
    """
    Same as TGK12_PLG but equations have been recalibrated with additional data.

    Line (1 and 7) of Table (1), Krein et al (2021, JGR: Solid Earth).
        
    Parameters
    ----------
    df : pandas dataframe
        The sample compositions.
    P : equilibration pressure
            
    Returns
    -------
    PT : pandas dataframe
       The calculated pressure(s) and temperature(s).
    """    
    def __init__(self, df):
        self.df = df
    
    def compute_pressure(self):
        pressure = - 43.59
        pressure += 136.9 * compute_olivine(self.df)
        pressure += 24.54 * Mg_num(self.df)
        pressure += 37.55 * NaK_num(self.df)
        pressure += 2.11 * self.df['TiO2_primary_wt_dry'] 
        pressure -= 0.773 * self.df['K2O_primary_wt_dry']         
        pressure += 1.636 * (self.df['CaO_primary_wt_dry']/self.df['Al2O3_primary_wt_dry']) 
        pressure *= 0.1              
        return pressure

    def compute_temperature(self, P):
        temperature = 1074.
        temperature += 11.86 * (P/0.1)
        temperature += 65.55 * Mg_num(self.df)
        temperature -= 138.2 * NaK_num(self.df)
        temperature += 20.55 * self.df['TiO2_primary_wt_dry'] 
        temperature += 5.855 * self.df['K2O_primary_wt_dry']        
        temperature += 79.02 * (self.df['CaO_primary_wt_dry'] / self.df['Al2O3_primary_wt_dry'])    
        return temperature
    
    def compute_pressure_temperature(self):
        compute_components_compound(self.df)
        P = self.compute_pressure()
        T = self.compute_temperature(P)
        return {'P': P, 'T': T} 

class HA15:
    """
    Method of Herzberg and O'Hara (2008, )
    
    Calculate temperate of olivine liquidus at 1 bar
    using Beattie (1993) Cointrib Min & Pet 115:103-111 
    Equations 10 and 12.
    
    Convert to desired pressure using Equation from 
    Herzberg & O'Hara, (2008) G3 9:9
    
    Parameters
    ----------
    df : pandas dataframe
        The sample compositions in mol fraction.
            
    Returns
    -------
    T: temperature 
    """     

    def __init__(self, df):
        self.df = df 
    
    def compute_temperature(self, P):
        compute_components_cation(self.df)
        Sum_A = (self.df["FeO_primary_mol"]*0.279)
        Sum_A += (self.df["MnO_primary_mol"]*0.259)
        Sum_A += (self.df["MgO_primary_mol"]*1.)
        Sum_A += (self.df["CaO_primary_mol"]*0.0056)
        Sum_A += (self.df["NiO_primary_mol"]*3.346)
        #####
        Sum_B = (self.df["FeO_primary_mol"]*0.031)
        Sum_B -= (self.df["MnO_primary_mol"]*0.049)
        Sum_B += (self.df["CaO_primary_mol"]*0.0135)
        Sum_B += (self.df["NiO_primary_mol"]*-3.665)
        Sum_B += 0.0001*(self.df["TiO2_primary_mol"]+self.df["Al2O3_primary_mol"]+self.df["Cr2O3_primary_mol"]+self.df["Fe2O3_primary_mol"]+self.df["Na2O_primary_mol"]+self.df["K2O_primary_mol"]+self.df["P2O5_primary_mol"])
        #####
        DMGO2 = (2./3. - Sum_B) / Sum_A
        #####
        Temp_Denom = 52.05/8.3143
        Temp_Denom += 2.*np.log(DMGO2)
        Temp_Denom += 2.*np.log(1.5*(self.df["FeO_primary_mol"]+self.df["MnO_primary_mol"]+self.df["MgO_primary_mol"]+self.df["CaO_primary_mol"]+self.df["NiO_primary_mol"]))
        Temp_Denom += 2.*np.log(3.*self.df["SiO2_primary_mol"])
        Temp_Denom -= 3.5*np.log(1-self.df["Al2O3_primary_mol"]) + 7*np.log(1-self.df["TiO2_primary_mol"])
        Temp_1bar = (113100/8.3143 + P*0.00000411/8.3143)
        Temp_1bar /= Temp_Denom
        Temp_1bar -= 273.15
        #####
        T = Temp_1bar + 54.*P - 2.*P**2.
        
        return {'P': P, 'T': T} 

class P07:

    def __init__(self, df):
        self.df = df
    
    def compute_temperature(self, P):
        compute_components_cation(self.df)
        C_NM = self.df["MgO_primary_mol"] + self.df["FeO_primary_mol"] 
        C_NM += self.df["CaO_primary_mol"] + self.df["MnO_primary_mol"]
        NF = 7/2*np.log(1.-self.df["Al2O3_primary_mol"]) 
        NF += 7.*np.log(1.-self.dfc["TiO2_primary_mol"])
        ln_DMg, ln_DFe, Temp_Pk = sym.symbols('ln_DMg, ln_DFe, Temp_Pk')
        eq1 = sym.Eq((self.df["MgO_primary_mol"] / ln_DMg) + (self.df["FeO_primary_mol"] / ln_DFe), 0.667)
        eq2 = sym.Eq(-2.158 + (55.09*(P/Temp_Pk)) - (6.213*10**-2 * self.df["H2O_primary_wt"]) + (4430. / Temp_Pk) + (5.115*10**-2 * (self.df["Na2O_primary_wt"]+self.df["K2O_primary_wt"])),ln_DMg)
        eq3 = sym.Eq(-3.3 + (47.57*(P/Temp_Pk)) - (5.192*10**-2 * self.df["H2O_primary_wt"]) + (3344. / Temp_Pk) + (5.595*10**-2 * (self.df["Na2O_primary_wt"]+self.df["K2O_primary_wt"])) + (1.633*10**-2 * self.df["SiO2_primary_wt"]), ln_DFe)
        result = sym.solve([eq1,eq2,eq3],(ln_DMg, ln_DFe, Temp_Pk))

        if result[1][2] < 1700.:
            T = result[1][2]
        elif result[0][2] < 1700.:
            T = result[0][2] 
            
        return {'P': P, 'T': T} 


def compute_sample_pressure_temperature(df, method="PF16"):
    
    if method == "PF16":
        out = PF16(df).compute_pressure_temperature()
    elif method == "L09":
        out = L09(df).compute_pressure_temperature()
    elif method == "TGK12_PLG":
        out = TGK12_PLG(df).compute_pressure_temperature()
    elif method == "TGK12_SPL":
        out = TGK12_SPL(df).compute_pressure_temperature()
    elif method == "G13":
        out = G13(df).compute_pressure_temperature()
    elif method == "K21_PLG":
        out = K21_PLG(df).compute_pressure_temperature()
    elif method == "K21_SPL":
        out = K21_SPL(df).compute_pressure_temperature()
    elif method == "K21_GNT":
        out = K21_GNT(df).compute_pressure_temperature()
    else:
        out = thermobar(df).compute_pressure_temperature()
        
    return out


def compute_sample_temperature(df, P, method="HA15"):
    if method == "HA15":
        out = HA15(df).compute_temperature()     
    elif method == "P07":
        out = P07(df).compute_temperature()            
    else:
        out = thermobar(df).compute_pressure_temperature()
        
    return out