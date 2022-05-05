from meltPT import *
import matplotlib.pyplot as plt


# ---- CAVP data
s = Suite("Central_Anatolia.csv", src_FeIII_totFe=0.2, min_MgO=8.5, min_SiO2=40.)
s.backtrack_compositions()
s.compute_pressure_temperature()


# ---- Set up mantle
lz = m.lithologies.katz.lherzolite()
mantle = m.mantle([lz], [1], ['Lz'])
max_P = 5.
P_sol = np.arange(0., max_P, 0.1)
T_sol = [lz.TSolidus(P) for P in P_sol]


# ---- Filter
K_weight = 39.0983
O_weight = 15.999
pc_to_ppm = 10000.
K_Nb_filter = lambda df: df['K2O']/df['Nb'] * (2.*K_weight) / (2.*K_weight + O_weight) * pc_to_ppm < 500.
where = np.where(K_Nb_filter(s.data))[0]
where_not = np.where(K_Nb_filter(s.data) == False)[0]


# ---- Find T
s.find_suite_potential_temperature(mantle, find_bounds=True, filters=(K_Nb_filter,), filter_args=(None,))


# ---- Plot
lab=r"$T_p$ = $%i^{+%i}_{-%i}$ $^\circ$C" % (
    s.potential_temperature, 
    s.upper_potential_temperature - s.potential_temperature,
    s.potential_temperature - s.lower_potential_temperature)
plt.plot(s.path.T, s.path.P, "--", label=lab)
plt.plot(s.upper_path.T, s.upper_path.P, ":")
plt.plot(s.lower_path.T, s.lower_path.P, ":")
plt.plot(T_sol, P_sol, "k")
plt.scatter(s.PT['T'].iloc[where], s.PT['P'].iloc[where], marker="o")
plt.scatter(s.PT['T'].iloc[where_not], s.PT['P'].iloc[where_not], marker="^")
plt.xlabel(r"Temperature [$^\circ$C]")
plt.ylabel("Pressure [GPa]")
plt.legend()
plt.gca().invert_yaxis()
plt.show()