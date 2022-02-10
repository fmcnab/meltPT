from meltPT import *

# ---- check against plank & forsyth supplementary
s = Suite("PF16_S7.csv", src_FeIII_totFe=0.19)
s.backtrack_compositions(Kd=0.3, verbose=True)
s.compute_pressure_temperature()
print()
print("Result from PF16 supplementary 7: P = 2.09 GPa, T = 1347 oC.")
print("Our result:                       P = %.2f GPa, T = %i oC." % 
    (s.PT['P'], s.PT['T']))


# ---- Set up mantle
lz = m.lithologies.katz.lherzolite()
mantle = m.Mantle([lz], [1], ['Lz'])
P_sol = np.arange(0., 3., 0.1)
T_sol = [lz.TSolidus(P) for P in P_sol]
path = mantle.AdiabaticMelt_1D(1330., Pstart=6., steps=101)


# ---- Perform fit to melt path
s.find_individual_potential_temperatures(mantle)
print()
print("Best-fitting melting model: Tp = %i oC, F = %.2f %%." % 
    (s.individual_potential_temperatures['Tp'], 
    s.individual_potential_temperatures['F_path_ind']*100.))
    

# ---- Plot
melt_label = "Melting path, $T_p = %i" % s.individual_potential_temperatures.iloc[0]['Tp']
plt.plot(T_sol, P_sol, "k", label="Solidus")
plt.plot(
    s.individual_potential_temperatures.iloc[0]['path'].T, 
    s.individual_potential_temperatures.iloc[0]['path'].P, 
    "--", label=melt_label)
plt.scatter(s.PT['T'], s.PT['P'], marker="*", label="Sample")
plt.xlabel(r"Temperature [$^\circ$C]")
plt.ylabel("Pressure [GPa]")
plt.legend()
plt.gca().invert_yaxis()
plt.show()
