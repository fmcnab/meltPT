from meltPT import *
import matplotlib.pyplot as plt

"""
Script to make Figure panels for Harrat Rahat
"""

df = pd.read_csv("Harrat_Rahat.csv", sep=',')
df1 = df.loc[(df['Province']=="Harrat_Rahat") & (df['Ce']>0.) & (df['Sm']>0.) & (df['Yb']>0.)]

############################
# ---- Split into Shield, Post-Shield and Rejuvenated Phases
# ---- Backtrack and Estimate PT
############################
# Calculate water based on H2O/Ce values
df1['H2O'] = 200. * df1['Ce'] / 10000
# Assign src_FeIII_totFe values
df1['src_FeIII_totFe'] = 0.18
#df1.loc[(df1['Unique_Id']=="Medinah"), 'src_FeIII_totFe'] = 0.24
#df1.loc[(df1['Unique_Id']=="Shawahit"), 'src_FeIII_totFe'] = 0.15
#df1.loc[(df1['Unique_Id']=="Hammah"), 'src_FeIII_totFe'] = 0.24

df1.to_csv("province.csv", sep=',')
s = Suite("province.csv", min_MgO=8.5)
s.backtrack_compositions()
s.compute_pressure_temperature(min_SiO2=40.)

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
s.find_suite_potential_temperature(mantle, find_bounds=True)

############################
# ---- Plot Hawaii Figure For Oahu
############################
fig, (ax1, ax2) = plt.subplots(1,2, sharey=True)

# ---- Plot Pressure vs. Temperature
ax1.plot(T_sol, P_sol, "k")
ax1.plot(s.path.T, s.path.P, "-", color="k", zorder=1)
ax1.plot(s.upper_path.T, s.upper_path.P, "--", color="k", zorder=1)
ax1.plot(s.lower_path.T, s.lower_path.P, "--", color="k", zorder=1)
ax1.text(0.95, 0.95, "$T_p$ = $%i^{+%i}_{-%i}$ $^\circ$C" % (
    s.potential_temperature, 
    s.upper_potential_temperature - s.potential_temperature,
    s.potential_temperature - s.lower_potential_temperature), verticalalignment='top', 
          horizontalalignment='right', transform=ax1.transAxes, fontsize=12)
ax1.text(0.05, 0.05, 'a', verticalalignment='bottom', 
          horizontalalignment='left', transform=ax1.transAxes, fontsize=12)
ax1.scatter(s.PT['T'][s.data['Unique_Id']=="Shawahit"], s.PT['P'][s.data['Unique_Id']=="Shawahit"],
            marker="^", facecolors="orange", edgecolor="k", zorder=2)
ax1.scatter(s.PT['T'][s.data['Unique_Id']=="Medinah"], s.PT['P'][s.data['Unique_Id']=="Medinah"],
            marker="s", facecolors="deeppink", edgecolor="k", zorder=2)


#ax1.scatter(s.PT['T'], s.PT['P'],
#            marker="^", facecolors="orange", edgecolor="k", zorder=2)
ax1.set_xlabel("Temperature [$^\circ$C]")
ax1.xaxis.set_tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)
ax1.xaxis.set_label_position('top')
ax1.set_ylabel("Pressure [GPa]")
ax1.set_xlim((1325.),(1725.))
ax1.set_ylim((0.5),(5.5))
ax1.invert_yaxis()

# ---- Plot Pressure vs. Sm/Yb
ax2.scatter(s.data['Sm']/s.data['Yb'], s.PT['P'], 
            marker="^", facecolors="orange", edgecolor="k")
ax2.set_xlabel("Sm/Yb")
ax2.xaxis.set_label_position('top')
ax2.xaxis.set_tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)
ax2.text(0.05, 0.05, 'b', verticalalignment='bottom', 
          horizontalalignment='left', transform=ax2.transAxes, fontsize=12)
ax2.set_ylim((0.5),(5.5))
ax2.invert_yaxis()

# ---- Finish Plot
plt.tight_layout()
plt.show()

# ---- Save to csv
s.write_to_csv("Harrat_Rahat_output.csv")
