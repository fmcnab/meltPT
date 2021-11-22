from MeltPT import *

# ---- north sea data
infile='North_Sea.csv'
outfile='North_Sea_PT.csv'
df = parse_csv(infile, src_FeIII_totFe=0.2, min_SiO2=40., min_MgO=8.5)
df = backtrack_compositions(df)
df = compute_pressure_temperature(df)


# ---- attempted melt fraction optimization

import pyMelt as m
from scipy.optimize import minimize, minimize_scalar

def melt_fraction_to_pressure_temperature(F, path):
    P = np.interp(F, path.F_total, path.P)
    T = np.interp(F, path.F_total, path.T)
    return P, T

def melt_fraction_misfit(F, df, path, full_output=True):
    model_P, model_T = melt_fraction_to_pressure_temperature(F, path)
    misfit = np.sqrt( (abs(df['P']-model_P)/0.24)**2. + (abs(df['T']-model_T)/39.)**2. )
    return misfit

def find_melt_fraction(df, path, full_output=True):
    if max([l.TSolidus(df['P']) for l in path.mantle.lithologies]) - df['T'] > 39.:
        return {'F': np.nan, 'misfit': np.nan}
    else:
        fit = minimize_scalar(melt_fraction_misfit, bounds=(0.,max(path.F_total)), bracket=(0.,max(path.F_total)), args=(df, path,False), method="bounded")
        if full_output:
            return {'F': fit.x, 'misfit': fit.fun}
        else:
            return fit.fun

def sample_potential_temperature_misfit(Tp, df, mantle):
    path = mantle.AdiabaticMelt_1D(Tp, Pstart=max(mantle.solidus_intersection(Tp))+0.01, steps=101)
    fit = minimize_scalar(melt_fraction_misfit, bounds=(0.,max(path.F_total)), bracket=(0.,max(path.F_total)), args=(df, path), method="bounded")
    return fit.fun

def find_sample_potential_temperature(df, mantle):
    Tp_fit = minimize_scalar(sample_potential_temperature_misfit, bracket=(lz.TSolidus(0.),1600.), bounds=(lz.TSolidus(0.),1600.), args=(df,mantle), method="bounded")
    path = mantle.AdiabaticMelt_1D(Tp_fit.x, Pstart=max(mantle.solidus_intersection(Tp_fit.x))+0.01, steps=101)
    F_fit = find_melt_fraction(df, path)
    return {'Tp': Tp_fit.x, 'F': F_fit['F'], 'path': path}

def potential_temperature_misfit(Tp, df, mantle):
    path = mantle.AdiabaticMelt_1D(Tp, Pstart=max(mantle.solidus_intersection(Tp))+0.01, steps=101)
    misfits = df.apply(find_melt_fraction, axis=1, result_type="expand", args=(path,True))['misfit']
    return np.nanmean(misfits)

def find_potential_temperature(df, mantle):
    Tp_fit = minimize_scalar(potential_temperature_misfit, bracket=(lz.TSolidus(0.),1600.), bounds=(lz.TSolidus(0.),1600.), args=(df, mantle), method="bounded")
    path = mantle.AdiabaticMelt_1D(Tp_fit.x, Pstart=max(mantle.solidus_intersection(Tp_fit.x))+0.01, steps=101)
    misfits = df.apply(find_melt_fraction, axis=1, result_type="expand", args=(path,True))
    return {'Tp': Tp_fit.x, 'F': misfits['F'], 'path': path}

# def find_potential_temperature(df, mantle, full_output=True):
#

def pandas_to_csv(df, outfile):
    # remove unwanted columns
    df = df.drop(['SiO2_primary_wt_dry','Al2O3_primary_wt_dry','FeO_primary_wt_dry','Fe2O3_primary_wt_dry','MgO_primary_wt_dry','CaO_primary_wt_dry','Na2O_primary_wt_dry','K2O_primary_wt_dry','TiO2_primary_wt_dry','MnO_primary_wt_dry','Cr2O3_primary_wt_dry','SiO2_primary_mol','Al2O3_primary_mol','FeO_primary_mol','Fe2O3_primary_mol','MgO_primary_mol','CaO_primary_mol','Na2O_primary_mol','K2O_primary_mol','TiO2_primary_mol','MnO_primary_mol','Cr2O3_primary_mol','H2O_primary_mol','Si4O8','Al16/3O8','Fe4Si2O8','Fe16/3O8','Mg4Si2O8','Ca4Si2O8','Na2Al2Si2O8','K2Al2Si2O8','Ti4O8','Mn4Si2O8','Cr16/3O8'], axis=1)
    # replace zeros with blanks
    df = df.replace(0, np.nan)
    # save to same csv as the input
    df.to_csv(outfile, index=False)
    return print("Script Finished")

def above_solidus(df, mantle):
    return {'fit': max([l.TSolidus(df['P']) for l in mantle.lithologies]) - df['T'] < 39.}

def combine(df):
    return {'fit': df.to_numpy().all()}

def samples_to_fit(df, mantle, filters=(None,), args=((None,))):
    to_fit = df.apply(above_solidus, axis=1, result_type="expand", args=(mantle,))
    if filters[0]:
        for f,a in zip(filters,args):
            to_fit = pd.concat([to_fit, df.apply(f, axis=1, args=a)], axis=1)
        to_fit = to_fit.apply(combine, axis=1, result_type="expand")
    return to_fit

to_fit = samples_to_fit(
    df,
    mantle,
    (lambda df,__: df['SiO2'] > 45., lambda df,__: df['MgO'] > 10.),
    ((None,), (None,))
    )


# ---- set up mantle domain
# lz = m.lithologies.matthews.klb1()
lz = m.lithologies.katz.lherzolite()
# lz = m.hydrous_lithology(lz, H2O=0.01)
mantle = m.Mantle([lz], [1], ['Lz'])
P_sol = np.arange(0., 6., 0.1)
T_sol = [lz.TSolidus(P) for P in P_sol]

# ---- example for single sample, fit F for specific Tp
path = mantle.AdiabaticMelt_1D(1330., Pstart=6., steps=101)
melt_fraction_fit = find_melt_fraction(df.iloc[0], path)
P_fit, T_fit = melt_fraction_to_pressure_temperature(melt_fraction_fit['F'], path)
plt.plot(path.T, path.P)
plt.plot(T_sol, P_sol, "--")
plt.plot([T_fit, df.iloc[0]['T']], [P_fit,df.iloc[0]['P']])
plt.scatter(df.iloc[0]['T'], df.iloc[0]['P'], marker="*")
plt.scatter(T_fit, P_fit)
plt.gca().invert_yaxis()
plt.show()


# ---- example for single sample, fit Tp
Tp_fit = find_sample_potential_temperature(df.iloc[0], mantle)
plt.plot(Tp_fit['path'].T, Tp_fit['path'].P)
plt.plot(T_sol, P_sol, "--")
plt.scatter(df.iloc[0]['T'], df.iloc[0]['P'], marker="*")
plt.gca().invert_yaxis()
plt.show()


# ---- example with all samples
Tp_full_fit = find_potential_temperature(df, mantle)
plt.plot(Tp_full_fit['path'].T, Tp_full_fit['path'].P)
plt.plot(T_sol, P_sol, "--")
plt.scatter(df['T'], df['P'])
plt.gca().invert_yaxis()
plt.show()

# # ---- print df to csv
# pandas_to_csv(df, outfile)
