from meltPT import *

# ---- North Sea data
s = Suite("North_Sea.csv", src_FeIII_totFe=0.2)
s.backtrack_compositions()
s.compute_pressure_temperature()

# ---- Set up mantle
lz = m.lithologies.katz.lherzolite()
mantle = m.Mantle([lz], [1], ['Lz'])
P_sol = np.arange(0., 6., 0.1)
T_sol = [lz.TSolidus(P) for P in P_sol]
path = mantle.AdiabaticMelt_1D(1330., Pstart=6., steps=101)

# ---- Fit individual sample
melt_fraction_fit = find_sample_melt_fraction(s.PT.iloc[0], path)
plt.plot(T_sol, P_sol, "k")
plt.plot(path.T, path.P, "--")
plt.plot(
    [melt_fraction_fit['T_path_ind'], s.PT.iloc[0]['T']], 
    [melt_fraction_fit['P_path_ind'], s.PT.iloc[0]['P']])
plt.scatter(s.PT.iloc[0]['T'], s.PT.iloc[0]['P'], marker="*")
plt.scatter(melt_fraction_fit['T_path_ind'], melt_fraction_fit['P_path_ind'])
plt.xlabel(r"Temperature [$^\circ$C]")
plt.ylabel("Pressure [GPa]")
plt.gca().invert_yaxis()
plt.show()

# ---- Apply to suite
s.find_individual_melt_fractions(mantle, path)
plt.plot(T_sol, P_sol, "k")
plt.plot(path.T, path.P, "--")
for i in range(len(s.PT)):
    plt.plot(
        [s.individual_melt_fractions.iloc[i]['T_path_ind'], s.PT.iloc[i]['T']], 
        [s.individual_melt_fractions.iloc[i]['P_path_ind'], s.PT.iloc[i]['P']])
    plt.scatter(s.PT.iloc[i]['T'], s.PT.iloc[i]['P'], marker="*")
    plt.scatter(s.individual_melt_fractions.iloc[i]['T_path_ind'], s.individual_melt_fractions.iloc[i]['P_path_ind'])
plt.xlabel(r"Temperature [$^\circ$C]")
plt.ylabel("Pressure [GPa]")
plt.gca().invert_yaxis()
plt.show()

# ---- Fit Tp to single sample
Tp_fit = find_sample_potential_temperature(s.PT.iloc[0], mantle)
plt.plot(T_sol, P_sol, "k")
plt.plot(Tp_fit['path'].T, Tp_fit['path'].P, "--")
plt.scatter(s.PT.iloc[0]['T'], s.PT.iloc[0]['P'], marker="*")
plt.xlabel(r"Temperature [$^\circ$C]")
plt.ylabel("Pressure [GPa]")
plt.gca().invert_yaxis()
plt.show()

# ---- Fit Tp to each sample
s.find_individual_potential_temperatures(mantle)
plt.plot(T_sol, P_sol, "k")
for i in range(len(s.PT)):
    if type(s.individual_potential_temperatures.iloc[i]['path']) is not float:
        plt.plot(s.individual_potential_temperatures.iloc[i]['path'].T, s.individual_potential_temperatures.iloc[i]['path'].P, "--")
plt.scatter(s.PT['T'], s.PT['P'], marker="*")
plt.xlabel(r"Temperature [$^\circ$C]")
plt.ylabel("Pressure [GPa]")
plt.gca().invert_yaxis()
plt.show()

# ---- Fit Tp to suite
s.find_suite_potential_temperature(mantle, find_bounds=True)
plt.plot(T_sol, P_sol, "k")
lab=r"$T_p$ = $%i^{+%i}_{-%i}$ $^\circ$C" % (
    s.potential_temperature, 
    s.upper_potential_temperature - s.potential_temperature,
    s.potential_temperature - s.lower_potential_temperature)
plt.plot(s.path.T, s.path.P, "--", label=lab)
plt.plot(s.upper_path.T, s.upper_path.P, ":")
plt.plot(s.lower_path.T, s.lower_path.P, ":")
plt.scatter(s.PT['T'], s.PT['P'], marker="*")
plt.xlabel(r"Temperature [$^\circ$C]")
plt.ylabel("Pressure [GPa]")
plt.legend()
plt.gca().invert_yaxis()
plt.show()

# ---- Write CSV
s.write_to_csv("North_Sea_output.csv")




