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
z_lab = 50.e3
P_lab = z_lab * density * g / 1.e9
T_lab = mantle.adiabat(P_lab, s.potential_temperature)
z = s.path.P / density / g * 1.e9
k_m = 4.5
k_c = 2.
z_moho = 20.e3
C_1_m = T_lab / (z_lab + z_moho*(k_m/k_c - 1.))
C_1_c = C_1_m * k_m / k_c
T = C_1_c * z
T[z > z_moho] = C_1_m * (z[z > z_moho] - z_lab) + T_lab

# ---- Make complete geotherm
combi_T = np.column_stack((s.path.T, T)).min(axis=1)

# ---- Xenoliths, Guest et al. 2020
Tx = np.array([1041., 1058., 1096., 1033., 1056.5, 984.])
Tx_err = np.array([31., 8., 19., 16., 8.5, 10.])
Px = np.array([1.46, 1.685, 1.685, 1.555, 1.555, 0.745])
Px_err = np.array([0.15, 0.115, 0.145, 0.185, 0.115, 0.035])

# ---- Plot
plt.plot(combi_T, s.path.P)
plt.plot(T_sol, P_sol, "k")
lab=r"$T_p$ = $%i^{+%i}_{-%i}$ $^\circ$C" % (
    s.potential_temperature, 
    s.upper_potential_temperature - s.potential_temperature,
    s.potential_temperature - s.lower_potential_temperature)
plt.plot(s.path.T, s.path.P, "--", label=lab)
plt.plot(s.upper_path.T, s.upper_path.P, ":")
plt.plot(s.lower_path.T, s.lower_path.P, ":")
plt.scatter(s.PT['T'], s.PT['P'], marker="*")
plt.errorbar(Tx, Px, xerr=Tx_err, yerr=Px_err, fmt="o")
plt.xlabel(r"Temperature [$^\circ$C]")
plt.ylabel("Pressure [GPa]")
plt.xlim((900.),T_sol[-1])
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

# # ---- Fit Tp to suite
# s.find_individual_potential_temperatures(mantle)
# hist = plt.hist(s.individual_potential_temperatures.Tp, bins=15)
# def gauss(x, mean, sd):
#     # return (1./(sd*np.sqrt(2.*np.pi))) * np.exp(-((x-mean)/(np.sqrt(2)*sd))**2.)
#     return np.exp(-((x-mean)/(np.sqrt(2)*sd))**2.)
# ts = np.linspace(hist[1].min(), hist[1].max(), 100)
# plt.plot(ts, hist[0].max()*gauss(ts, s.potential_temperature, s.upper_potential_temperature-s.potential_temperature))
# plt.show()
# 
# fig, (ax1, ax2) = plt.subplots(1,2)
# 
# ax1.scatter( s.data['La']/s.data['Sm'], s.suite_melt_fractions['F'] )
# ax1.set_ylabel("F (Suite)")
# ax1.set_xlabel("La/Sm")
# ax1.set_box_aspect(1)
# 
# ax2.scatter( s.data['La']/s.data['Sm'], s.individual_potential_temperatures['F'] )
# ax2.set_ylabel("F (Individual)")
# ax2.set_xlabel("La/Sm")
# ax2.set_box_aspect(1)
# 
# plt.show()
