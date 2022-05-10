import warnings

import numpy as np
import pandas as pd

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
        
        pressure = (
            (np.log(self.df['Si4O8']) - 
                4.019 + 
                0.0165*self.df['Fe4Si2O8'] + 
                0.0005*(self.df['Ca4Si2O8']**2.)) /
            (-770*(T**(-1.)) + 0.0058*(T**0.5) - 0.003*self.df['H16O8'])
            )
        return pressure
    
    def compute_temperature(self):
        
        temperature = (
            916.45 + 
            (13.68 * self.df['Mg4Si2O8']) + 
            (4580. / self.df['Si4O8']) -
            (0.509 * self.df['H16O8'] * self.df['Mg4Si2O8'])
            )
        return temperature
        
    def compute_pressure_temperature(self):
        T = self.compute_temperature()
        P = self.compute_pressure(T)
        return {'P': P, 'T': T}

def compute_sample_pressure_temperature(df, method="PF16"):
    
    if method == "PF16":
        out = PF16(df).compute_pressure_temperature()
    elif method == "L09":
        out = L09(df).compute_pressure_temperature()
    else:
        out = thermobar(df).compute_pressure_temperature()
        
    return out