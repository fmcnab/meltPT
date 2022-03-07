from .parse import *
from .backtrack_compositions import *
from .thermobarometry import *
from .fit_melting_paths import *

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
        Backtrack compositions for entire suite.
    compute_pressure_temperature :
        Compute equilibration pressures and temperatures for entire suite.
    check_samples_for_fitting :
        Determine which samples should be fitted.
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

    def compute_pressure_temperature(self):
        """
        Compute equilibration pressures and temperatures for entire suite.
        
        Applies "compute_sample_pressure_temperature" to the "primary"
        property. Result is saved in the "PT" property.
        """
        self.PT = self.primary.apply(
            compute_sample_pressure_temperature, 
            axis=1, 
            result_type="expand")
        
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
        self.individual_melt_fractions = self.PT_to_fit.apply(
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
        self.individual_potential_temperatures = self.PT_to_fit.apply(
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

        self.upper_potential_temperature, self.upper_path = find_bounding_potential_temperature(upper_points, self.potential_temperature, mantle, threshold=bounds_threshold)
        self.lower_potential_temperature, self.lower_path = find_bounding_potential_temperature(lower_points, self.potential_temperature, mantle, lower=True, threshold=bounds_threshold)

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