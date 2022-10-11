
"""
===============
The Suite class
===============

Process suites of samples.

"""

from .parse import *
from .backtrack_compositions import *
from .thermobarometry import *
from .fit_melting_paths import *

MAJOR_OXIDES = [
    'SiO2','Al2O3','FeO','Fe2O3','MgO','CaO','Na2O','K2O','TiO2','MnO',
    'Cr2O3', 'P2O5', 'NiO', 'CoO', 'H2O', 'CO2']

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
    Ce_to_H2O : float
        Ratio of Ce to H2O in mantle source.
    src_FeIII_totFe : float
        Ratio of Fe3+ to total Fe in the mantle source.
    min_MgO : float, optional
        Minimum amount of MgO in sample to be accepted.
    read_as_primary : bool
        If true, data from input_csv is assumed to be primary
        and backtracking is avoided.
    param_co2 : bool
        If true, CO2 is calculated from SiO2 concentration.
    read_PT : bool
        If true, pressures and temperatures are not calculated.
        
    Attributes
    ----------
    data : pandas dataframe
        The raw data from the provided csv.
    primary : pandas dataframe or NoneType
        If backtrack_compositions has been run, or read_as_primary is set 
        to true, will contain the estimated primary compositions in 
        various forms.
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
    """

    def __init__(self, input_csv, Ce_to_H2O=0., src_FeIII_totFe=0., min_MgO=0., read_as_primary=False, param_co2=False, read_PT=False):
        self.data = parse_csv(
            input_csv,
            Ce_to_H2O=Ce_to_H2O,
            src_FeIII_totFe=src_FeIII_totFe,
            min_MgO=min_MgO,
            param_co2=param_co2)
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
        
        if read_as_primary:
            self.primary = self.data.filter(items=MAJOR_OXIDES, axis=1)
            primary_labels = {x: x + "_primary_wt" for x in MAJOR_OXIDES}
            self.primary = self.primary.rename(primary_labels, axis=1)
        else:
            self.primary = None
            
        if read_PT:
            self.PT = self.data.filter(items=['P', 'T'], axis=1)
            self.data = self.data.drop(labels=['P', 'T'], axis=1)
        else:
            self.PT = None
            self.data = self.data.drop(labels=['P', 'T'], axis=1, errors='ignore')

    def backtrack_compositions(self, target_Fo=0.9, Kd=False, dm=0.0005, verbose=False, max_olivine_addition=0.3):
        """
        Backtrack compositions for entire suite.
        
        Applies backtrack_sample_composition to the "data" property. Result is
        saved in the "primary" property.
        
        Parameters
        ----------
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
            Maximum fraction of olivine to add before backtracking is
            abandoned.
        """
        self.primary = self.data.apply(
            backtrack_sample_composition,
            axis=1,
            result_type="expand",
            args=(target_Fo,Kd,dm,verbose)
            )

    def compute_pressure_temperature(self, method="PF16", min_SiO2=0.):
        """
        Compute equilibration pressures and temperatures for entire suite.
        
        Applies "compute_sample_pressure_temperature" to the "primary"
        property. Result is saved in the "PT" property.
        
        Parameters
        ----------
        method : str
            Code corresponding to desired thermobarometer.
            Current options are:
              - P08: Putirka et al. (2008, Revs. in Min. and Geo.)
              - L09: Lee et al. (2009, EPSL)
              - TGK12_PLG: Till et al. (2012, JGR: Solid Earth), plagioclase
              - TGK12_SPL: Till et al. (2012, JGR: Solid Earth), spinel
              - PF16: Plank and Forsyth (2016, G-cubed)
              - SD20: Sun and Dasgupta (2020, EPSL)
              - BK21: Brown Krien et al. (2021, JGR: Solid Earth), stable phase
              - BK21_PLG: Brown Krien et al. (2021, JGR: Solid Earth), plagioclase
              - BK21_SPL: Brown Krien et al. (2021, JGR: Solid Earth), spinel
              - BK21_GNT: Brown Krien et al. (2021, JGR: Solid Earth), garnet
        min_SiO2 : float
            Threshold SiO2 content below which samples will be ignored.
        """
        self.PT = self.primary.apply(
            compute_sample_pressure_temperature, 
            axis=1, 
            result_type="expand",
            args=(method,min_SiO2,))
    
    def compute_temperature(self, method="PF16", P=1., min_SiO2=0.):
        """
        Compute equilibration temperatures for entire suite at given pressure.
        
        Applies "compute_sample_temperature" to the "primary"
        property. Result is saved in the "PT" property.
        
        Parameters
        ----------
        method : str
            Code corresponding to desired thermometer.
            Current options are:
              - P08: Putirka et al. (2008, Revs. in Min. and Geo.)
              - L09: Lee et al. (2009, EPSL)
              - TGK12_PLG: Till et al. (2012, JGR: Solid Earth), plagioclase
              - TGK12_SPL: Till et al. (2012, JGR: Solid Earth), spinel
              - PF16: Plank and Forsyth (2016, G-cubed)
              - SD20: Sun and Dasgupta (2020, EPSL)
              - BK21: Brown Krien et al. (2021, JGR: Solid Earth), stable phase
              - BK21_PLG: Brown Krien et al. (2021, JGR: Solid Earth), plagioclase
              - BK21_SPL: Brown Krien et al. (2021, JGR: Solid Earth), spinel
              - BK21_GNT: Brown Krien et al. (2021, JGR: Solid Earth), garnet
              - B93: Beattie (1993, Contrib. to Min. and Pet.)
              - P07_2: Putirka et al. (2007, Chem. Geol.), Equation 2
              - P07_4: Putirka et al. (2007, Chem. Geol.), Equation 4
              - HA15: Herzberg and Asimow (2015, G-cubed)
        P : float
            Pressure at which temperature whould be calculated.
        min_SiO2 : float
            Threshold SiO2 content below which samples will be ignored.
        """
        self.PT = self.primary.apply(
            compute_sample_temperature, 
            axis=1, 
            result_type="expand",
            args=(method,P,min_SiO2,)
            )    

    def compute_pressure(self, method="PF16", T=1300., min_SiO2=0.):
        """
        Compute equilibration pressures for entire suite at given temperature.
        
        Applies "compute_sample_pressure" to the "primary"
        property. Result is saved in the "PT" property.
        
        Parameters
        ----------
        method : str
            Code corresponding to desired barometer.
            Current options are:
              - P08: Putirka et al. (2008, Revs. in Min. and Geo.)
              - L09: Lee et al. (2009, EPSL)
              - TGK12_PLG: Till et al. (2012, JGR: Solid Earth), plagioclase
              - TGK12_SPL: Till et al. (2012, JGR: Solid Earth), spinel
              - PF16: Plank and Forsyth (2016, G-cubed)
              - SD20: Sun and Dasgupta (2020, EPSL)
              - BK21: Brown Krien et al. (2021, JGR: Solid Earth), stable phase
              - BK21_PLG: Brown Krien et al. (2021, JGR: Solid Earth), plagioclase
              - BK21_SPL: Brown Krien et al. (2021, JGR: Solid Earth), spinel
              - BK21_GNT: Brown Krien et al. (2021, JGR: Solid Earth), garnet
        T : float
            Temperature at which pressure should be calculated.
        min_SiO2 : float
            Threshold SiO2 content below which samples will be ignored.
        """
        self.PT = self.primary.apply(
            compute_sample_pressure, 
            axis=1, 
            result_type="expand",
            args=(method,T,min_SiO2,)
            )  
        
    def check_samples_for_fitting(self, mantle, filters=(None,), args=((None,))):
        """
        Determine which samples should be fitted.
        
        Samples below the solidus will be dropped. Also applies additional
        filters if provided by the user. Samples which fail any filter will
        be rejected.
        
        Result is saved in "PT_to_fit" property. Same as "PT" but rejected
        samples are assigned nan pressure and temperature.
        
        Parameters
        ----------
        mantle : instance of pyMelt.mantle_class.mantle
            The mantle object to be used to calculate the solidus.
        filters : tuple of functions, optional
            Set of functions to filter samples before fitting.
            Filters take the form of a function that reads a single-row
            dataframe and returns either true or force.
        args : tuple of tuples, optional, same length as filters
            Extra arguments to be passed to the filter functions.
        """
        def combine(df):
            return {'fit': df.to_numpy().all()}

        def above_solidus(df, mantle):
            return {'fit': max([l.TSolidus(df['P']) for l in mantle.lithologies]) - df['T'] < df['T_err']}

        def isnotnan(df):
            return {'fit': ~np.isnan(np.array([df['P'], df['T']])).any()}

        combined = pd.concat([self.data, self.PT], axis=1)
        to_fit = combined.apply(above_solidus, axis=1, result_type="expand", args=(mantle,))
        to_fit = pd.concat([to_fit, combined.apply(isnotnan, axis=1, result_type="expand")], axis=1)
        if filters[0]:
            for f,a in zip(filters,args):
                to_fit = pd.concat([to_fit, combined.apply(f, axis=1, args=a, result_type="expand")], axis=1)
        to_fit = to_fit.apply(combine, axis=1, result_type="expand")
        self.PT['Fit_Tp'] = to_fit.copy()
    
    def find_individual_melt_fractions(self, mantle, path, filters=(None,), filter_args=(None,)):
        """
        Find best-fitting melt fractions for each sample relative to given melt
        path.
        
        Applies "find_sample_melt_fraction" to "PT_to_fit" dataframe. Results
        saved in "individual_melt_fractions" property.
        
        Parameters
        ----------
        mantle : instance of pyMelt.mantle_class.mantle
            The mantle object to be used to calculate the solidus.
        path : instance of pyMelt.meltingcolumn_classes.meltingColumn
            The melting path to be used.
        filters : tuple of functions, optional
            Set of functions to filter samples before fitting.
            Filters take the form of a function that reads a single-row
            dataframe and returns either true or force.
        args : tuple of tuples, optional, same length as filters
            Extra arguments to be passed to the filter functions.
        """
        self.check_samples_for_fitting(mantle, filters, filter_args)
        self.individual_melt_fractions = self.PT.apply(
            find_sample_melt_fraction, 
            axis=1, 
            result_type="expand", 
            args=(path,))

    def find_individual_potential_temperatures(self, mantle, filters=(None,), filter_args=(None,)):
        """
        Find best-fitting melting paths for each sample.
        
        Applies "find_sample_potential_temperature" to "PT_to_fit" dataframe.
        Result is saved in "individual_potential_temperatures" property.

        Parameters
        ----------
        mantle : instance of pyMelt.mantle_class.mantle
            The mantle object to be used to calculate the solidus and melting
            paths.
        filters : tuple of functions, optional
            Set of functions to filter samples before fitting.
            Filters take the form of a function that reads a single-row
            dataframe and returns either true or force.
        args : tuple of tuples, optional, same length as filters
            Extra arguments to be passed to the filter functions.
        """
        self.check_samples_for_fitting(mantle, filters, filter_args)
        self.individual_potential_temperatures = self.PT.apply(
            find_sample_potential_temperature, 
            axis=1, 
            result_type="expand", 
            args=(mantle,))

    def find_suite_potential_temperature(self, mantle, find_bounds=False, bounds_threshold=(2./3.), filters=(None,), filter_args=(None,)):
        """
        Find best-fitting potential temperature for suite.
        
        Uses scipy's minimize_scalar to minimize mean distance from suite
        pressure-temperature estimates to melting path. Options used are:
            method: "bounded"
            bounds:  min --> intersection of solidus with surface
                     max --> 1600 oC
            bracket: min --> intersection of solidus with surface
                     max --> 1600 oC
        
        Parameters
        ----------
        mantle : instance of pyMelt.mantle_class.mantle
            The mantle object to be used to calculate the solidus and melting
            paths.
        find_bounds : bool, optional
            If True, uses find_bounds to place upper and lower bounds on
            best-fitting potential temperature.
        bounds_threshold : float
            The threshold fraction of points to be lie between the best-fitting
            and bounding temperature melting paths.
        filters : tuple of functions, optional
            Set of functions to filter samples before fitting.
            Filters take the form of a function that reads a single-row
            dataframe and returns either true or force.
        args : tuple of tuples, optional, same length as filters
            Extra arguments to be passed to the filter functions.         
        
        Results (saved as class properties)
        -------
        potential_temperature : float
            The best-fitting potential temperature.
        path : instance of pyMelt.meltingcolumn_classes.meltingColumn
            The best-fitting melting path.
        suite_melt_fractions : pandas dataframe
            Nearest melt fractions, pressures and temperatures along best-
            fitting melting path for each sample.
        upper_potential_temperature : float, if find_bounds is True
            The upper bouding potential temperature.
        lower_potential_temperature : float, if find_bounds is True
            The lower bounding potential temperature.
        upper_path : instance of pyMelt.meltingcolumn_classes.meltingColumn, 
                        if finds_bounds is True
            The upper bounding melting path.
        lower_path : instance of pyMelt.meltingcolumn_classes.meltingColumn, 
                        if finds_bounds is True
            The lower bounding melting path.
        """
        self.check_samples_for_fitting(mantle, filters, filter_args)
        max_Tp = find_max_potential_temperature(mantle)
        Tp_fit = minimize_scalar(
            compute_suite_potential_temperature_misfit, 
            bracket=(min([lith.TSolidus(0.) for lith in mantle.lithologies]),max_Tp), 
            bounds=(min([lith.TSolidus(0.) for lith in mantle.lithologies]),max_Tp),
            args=(self.PT, mantle),
            method="bounded")
        self.path = mantle.adiabaticMelt(Tp_fit.x, Pstart=max(mantle.solidusIntersection(Tp_fit.x))+0.01, dP=-0.01)
        self.suite_melt_fractions = self.PT.apply(find_sample_melt_fraction, axis=1, result_type="expand", args=(self.path,))
        self.potential_temperature = Tp_fit.x
        
        if find_bounds:
            upper_points = []
            lower_points = []
            for i in range(len(self.PT)):
                if self.PT['T'].iloc[i] - self.suite_melt_fractions['T'].iloc[i] > 0:
                    upper_points.append(shp.Point(self.PT['T'].iloc[i], self.PT['P'].iloc[i]))
                elif self.PT['T'].iloc[i] - self.suite_melt_fractions['T'].iloc[i] < 0.:
                    lower_points.append(shp.Point(self.PT['T'].iloc[i], self.PT['P'].iloc[i]))

            self.upper_potential_temperature, self.upper_path = find_bounding_potential_temperature(upper_points, self.potential_temperature, mantle, threshold=bounds_threshold)
            self.lower_potential_temperature, self.lower_path = find_bounding_potential_temperature(lower_points, self.potential_temperature, mantle, lower=True, threshold=bounds_threshold)

    def write_to_csv(self, outfile, write_primary=True, write_PT=True, 
                        write_suite_Tp=False, write_individual_Tp=False):
        """
        Write results to csv.
        
        Combine desired outputs, give appropriate column names, and save as
        csv to a specified location.
        
        Parameters
        ----------
        outfile : str
            Path to location where csv should be saved.
        write_primary : bool
            Whether hydrous wt% concentrations should be saved.
        write_PT : bool
            Whether equilibration pressure/temperature estimates should be
            saved.
        write_suite_Tp : bool
            Whether results of fitting melting paths through suite of pressure/
            temperature estimates should be saved.
        write_individual_Tp : bool
            Whether results of fitting melting paths through individual
            pressure/temperature points should be saved.
        """
        output_df = self.data.copy()
        if write_primary and self.primary is not None:
            output_df = pd.concat([output_df, self.primary], axis=1)
        if write_PT and self.PT is not None:
            output_df = pd.concat([output_df, self.PT], axis=1)
        if write_suite_Tp and self.suite_melt_fractions is not None:
            rename_dict = {
                'F': 'F_suite_fit',
                'P': 'P_suite_fit',
                'T': 'T_suite_fit',
                'misfit': 'misfit_suite_fit'
            }
            suite_out = self.suite_melt_fractions.rename(columns=rename_dict)
            suite_out['Tp_suite_fit'] = self.potential_temperature * self.PT['Fit_Tp']
            suite_out['Tp_max_suite_fit'] = self.upper_potential_temperature * self.PT['Fit_Tp']
            suite_out['Tp_min_suite_fit'] = self.lower_potential_temperature * self.PT['Fit_Tp']
            suite_out = suite_out.replace(0, np.nan)
            output_df = pd.concat([output_df, suite_out], axis=1)
        if write_individual_Tp and self.individual_potential_temperatures is not None:
            rename_dict = {
                'F': 'F_ind_fit',
                'Tp': 'Tp_ind_fit'
            }
            ind_out = self.individual_potential_temperatures.rename(columns=rename_dict)
            ind_out = ind_out.drop(["P", "T", "path", "misfit"], axis=1)
            output_df = pd.concat([output_df, ind_out], axis=1)
        output_df.to_csv(outfile, index=False)