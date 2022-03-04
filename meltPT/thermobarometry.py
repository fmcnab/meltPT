import warnings

import numpy as np
import pandas as pd

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
