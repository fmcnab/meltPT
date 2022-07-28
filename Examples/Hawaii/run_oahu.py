from meltPT import *
import matplotlib.pyplot as plt

"""
Script to make Hawaii Figure panels for Oahu
"""

############################
# ---- Isolate Oahu Data
# ---- Split into Shield, Post-Shield and Rejuvenated Phases
# ---- Backtrack and Estimate PT
############################

# ---- Shield Data
df = pd.read_csv("Hawaii.csv", sep=',')
df1 = df.loc[(df['Stage']=="Shield") & (df['Province']=="Oahu") & (df['Sm']>0) & (df['Yb']>0)]
df1.to_csv("province.csv", sep=',')
s = Suite("province.csv", src_FeIII_totFe=0.15, Ce_to_H2O=144., min_MgO=8.5, min_SiO2=40.)
s.backtrack_compositions()
s.compute_pressure_temperature()

# ---- Post-Shield Data
df = pd.read_csv("Hawaii.csv", sep=',')
df1 = df.loc[(df['Stage']=="Post-Shield") & (df['Province']=="Oahu") & (df['Sm']>0) & (df['Yb']>0)]
df1.to_csv("province.csv", sep=',')
s2 = Suite("province.csv", src_FeIII_totFe=0.15, Ce_to_H2O=136., min_MgO=8.5, min_SiO2=40.)
s2.backtrack_compositions()
s2.compute_pressure_temperature()

# ---- Rejuvenated Data
df = pd.read_csv("Hawaii.csv", sep=',')
df1 = df.loc[(df['Stage']=="Rejuvenated") & (df['Province']=="Oahu") & (df['Sm']>0) & (df['Yb']>0)]
df1.to_csv("province.csv", sep=',')
s3 = Suite("province.csv", src_FeIII_totFe=0.17, Ce_to_H2O=211., min_MgO=8.5, min_SiO2=40.)
s3.backtrack_compositions()
s3.compute_pressure_temperature()

# ---- Combine Phases and Estimate PT
df_hawaii = pd.concat([s.primary,s2.primary,s3.primary])

hawaii_dict = {'SiO2' : df_hawaii['SiO2_primary_wt'],
               'TiO2' : df_hawaii['TiO2_primary_wt'],
               'Al2O3' : df_hawaii['Al2O3_primary_wt'],
               'FeO' : df_hawaii['FeO_primary_wt'],
               'Fe2O3' : df_hawaii['Fe2O3_primary_wt'],
               'MnO' : df_hawaii['MnO_primary_wt'],
               'MgO' : df_hawaii['MgO_primary_wt'],
               'CaO' : df_hawaii['CaO_primary_wt'],
               'Na2O' : df_hawaii['Na2O_primary_wt'],
               'K2O' : df_hawaii['K2O_primary_wt'],
               'P2O5' : df_hawaii['P2O5_primary_wt'],
               'H2O' : df_hawaii['H2O_primary_wt']}
df_hawaii = pd.DataFrame.from_dict(hawaii_dict)
df_hawaii.to_csv("Oahu.csv", sep=',')

s4 = Suite("Oahu.csv", read_as_primary=True)
s4.compute_pressure_temperature()

############################
# ---- Set up mantle
############################
lz = m.lithologies.katz.lherzolite()
mantle = m.mantle([lz], [1], ['Lz'])
max_P = -lz.parameters['A2'] / (2.*lz.parameters['A3'])
P_sol = np.arange(0., max_P, 0.1)
T_sol = [lz.TSolidus(P) for P in P_sol]

############################
# ---- Fit Tp to Oahu Data
############################
s4.find_suite_potential_temperature(mantle, find_bounds=True)

############################
# ---- Plot Oahu Xenoliths
############################
# # ---- Xenoliths, Guest et al. 2020
# Tx = np.array([1041., 1058., 1096., 1033., 1056.5, 984.])
# Tx_err = np.array([31., 8., 19., 16., 8.5, 10.])
# Px = np.array([1.46, 1.685, 1.685, 1.555, 1.555, 0.745])
# Px_err = np.array([0.15, 0.115, 0.145, 0.185, 0.115, 0.035])

############################
# ---- Plot Hawaii Figure For Oahu
############################
fig, ((ax1, ax2),(ax3, ax4)) = plt.subplots(2,2, sharey=True)

# ---- Plot Pressure vs. Temperature
ax1.plot(T_sol, P_sol, "k")
ax1.plot(s4.path.T, s4.path.P, "-", color="k", zorder=1)
ax1.plot(s4.upper_path.T, s4.upper_path.P, "--", color="k", zorder=1)
ax1.plot(s4.lower_path.T, s4.lower_path.P, "--", color="k", zorder=1)
ax1.text(0.95, 0.95, "$T_p$ = $%i^{+%i}_{-%i}$ $^\circ$C" % (
    s4.potential_temperature, 
    s4.upper_potential_temperature - s4.potential_temperature,
    s4.potential_temperature - s4.lower_potential_temperature), verticalalignment='top', 
         horizontalalignment='right', transform=ax1.transAxes, fontsize=12)
ax1.text(0.05, 0.05, 'a', verticalalignment='bottom', 
         horizontalalignment='left', transform=ax1.transAxes, fontsize=12)
ax1.scatter(s.PT['T'], s.PT['P'], marker="^", facecolors="orange", edgecolor="k", zorder=2)
ax1.scatter(s3.PT['T'], s3.PT['P'], marker="o", facecolors="dodgerblue", edgecolor="k", zorder=2)
ax1.set_xlabel("Temperature [$^\circ$C]")
ax1.xaxis.set_tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)
ax1.xaxis.set_label_position('top')
ax1.set_ylabel("Pressure [GPa]")
ax1.set_xlim((1325.),(1725.))
ax1.set_ylim((0.5),(5.5))
ax1.invert_yaxis()

# ---- Plot Pressure vs. La/Sm
ax2.scatter( s.data['Sm']/s.data['Yb'], s.PT['P'], marker="^", facecolors="orange", edgecolor="k")
ax2.scatter( s3.data['Sm']/s3.data['Yb'], s3.PT['P'], marker="o", facecolors="dodgerblue", edgecolor="k")
ax2.set_xlabel("Sm/Yb")
ax2.xaxis.set_label_position('top')
ax2.xaxis.set_tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)
ax2.text(0.05, 0.05, 'b', verticalalignment='bottom', 
         horizontalalignment='left', transform=ax2.transAxes, fontsize=12)
ax2.set_ylim((0.5),(5.5))
ax2.invert_yaxis()

# ---- Plot Pressure vs. Nb/Zr
ax3.scatter( s.data.loc[s.data['Nb']>0., 'Nb']/s.data.loc[s.data['Nb']>0., 'Zr'], s.PT.loc[s.data['Nb']>0., 'P'], marker="^", facecolors="orange", edgecolor="k")
ax3.scatter( s3.data.loc[s3.data['Nb']>0., 'Nb']/s3.data.loc[s3.data['Nb']>0., 'Zr'], s3.PT.loc[s3.data['Nb']>0., 'P'], marker="o", facecolors="dodgerblue", edgecolor="k")
ax3.set_xlabel("Nb/Zr")
ax3.set_ylabel("Pressure [GPa]")
ax3.text(0.05, 0.05, 'c', verticalalignment='bottom', 
         horizontalalignment='left', transform=ax3.transAxes, fontsize=12)
ax3.set_ylim((0.5),(5.5))
ax3.invert_yaxis()

# ---- Plot Pressure vs. Nd Isotopes
ax4.scatter( s.data.loc[s.data['143Nd/144Nd']>0., '143Nd/144Nd'], s.PT.loc[s.data['143Nd/144Nd']>0., 'P'], marker="^", facecolors="orange", edgecolor="k")
ax4.scatter( s3.data.loc[s3.data['143Nd/144Nd']>0., '143Nd/144Nd'], s3.PT.loc[s3.data['143Nd/144Nd']>0., 'P'], marker="o", facecolors="dodgerblue", edgecolor="k")
ax4.set_xlabel("$^{143}$Nd/$^{144}$Nd")
ax4.text(0.05, 0.05, 'd', verticalalignment='bottom', 
         horizontalalignment='left', transform=ax4.transAxes, fontsize=12)
ax4.set_ylim((0.5),(5.5))
ax4.invert_yaxis()

# ---- Finish Plot
plt.tight_layout()
plt.show()