import warnings

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
import shapely.geometry as shp

import pyMelt as m


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
        bounds:  min --> intersection of solidus with surface
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