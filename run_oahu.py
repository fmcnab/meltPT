from meltPT import *
import matplotlib.pyplot as plt

# ---- Oahu data
s = Suite("Oahu.csv", src_FeIII_totFe=0.15, min_MgO=8.5, min_SiO2=40.)
s.backtrack_compositions()
s.compute_pressure_temperature()


# ---- Set up mantle
lz = m.lithologies.katz.lherzolite()
# lz.parameters['Mcpx'] = 0.15 # shorttle
# lz.CP = 1187. # shorttle
# lz.alphas = 30. # shorttle
# lz.DeltaS = 407. # shorttle
mantle = m.mantle([lz], [1], ['Lz'])
# max_P = -lz.parameters['A2'] / (2.*lz.parameters['A3'])
max_P = 7.
P_sol = np.arange(0., max_P, 0.1)
T_sol = [lz.TSolidus(P) for P in P_sol]

# ---- Fit Tp to suite
s.find_suite_potential_temperature(mantle, find_bounds=True)

# ---- Conductive geotherm stuff
density = 3300.
g = 9.81
z_lab = 60.e3
P_lab = z_lab * density * g / 1.e9
T_lab = mantle.adiabat(P_lab, s.potential_temperature)

# ---- Plot
plt.plot([0.,T_lab], [0.,P_lab])
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
plt.xlim((T_sol[0]),T_sol[-1])
plt.legend()
plt.gca().invert_yaxis()
plt.show()

# ---- Plot Pressure vs. La/Sm

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2,sharey="row")

ax1.scatter( s.data['Sm']/s.data['Yb'], s.PT['P'] )
ax1.set_ylabel("P [Gpa]")
ax1.set_xlabel("Sm/Yb")
ax1.set_box_aspect(1)
ax1.invert_yaxis()

inds = np.where(s.data['Nb'] != 0)[0]
ax2.scatter( s.data['Nb'].iloc[inds]/s.data['Zr'].iloc[inds], s.PT['P'].iloc[inds] )
ax2.set_xlabel("Nb/Zr")
ax2.set_box_aspect(1)

ax3.scatter( s.data['La']/s.data['Sm'], s.suite_melt_fractions['F'] )
ax3.set_ylabel("F")
ax3.set_xlabel("La/Sm")
ax3.set_box_aspect(1)

inds = np.where(s.data['Nb'] != 0)[0]
ax4.scatter( s.data['Nb'].iloc[inds]/s.data['Zr'].iloc[inds], s.suite_melt_fractions['F'].iloc[inds] )
ax4.set_xlabel("Nb/Zr")
ax4.set_box_aspect(1)

plt.tight_layout()
plt.show()

