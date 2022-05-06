from meltPT import *
import matplotlib.pyplot as plt

# ---- Import All Hawaii Data

#df = pd.read_csv("Hawaii.csv", sep=',')
#df1 = df.loc[df['Province']=="Lanai"]
#df = df.loc[df['Province']=="Kahoolawe"]
#df = df.append(df1)
#df.to_csv("province.csv", sep=',')

# ---- Province data
s = Suite("Hawaii.csv", src_FeIII_totFe=0.15, min_MgO=8.5, min_SiO2=40.)
s.backtrack_compositions()
s.compute_pressure_temperature()


# ---- Set up mantle
lz = m.lithologies.katz.lherzolite()
# lz.parameters['Mcpx'] = 0.15 # shorttle
# lz.CP = 1187. # shorttle
# lz.alphas = 30. # shorttle
# lz.DeltaS = 407. # shorttle
mantle = m.mantle([lz], [1], ['Lz'])
max_P = -lz.parameters['A2'] / (2.*lz.parameters['A3'])
P_sol = np.arange(0., max_P, 0.1)
T_sol = [lz.TSolidus(P) for P in P_sol]


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

s.write_to_csv("Hawaii_output.csv")

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
