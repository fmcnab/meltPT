#------------------------------------------------------------------------------
# ---- Figure 6
#
# Recreate elements of Figure 6 from McNab & Ball (2023).
# Comment / uncomment lines 21/22 to plot Loa/Kea trends.
# This will only work with v1.1.0!
#------------------------------------------------------------------------------

# ---- Importing
from meltPT import *
import matplotlib.pyplot as plt

# ---- Isolate Island Data
# Read in compilation dataset of Hawaiian samples from McNab & Ball, (2022)
df = pd.read_csv("Dataset_S1.csv", sep=',')

# We want take only samples from our Hawaii database that correspond to
# the Mauna Kea Trend: Moloka`i, Maui and the Hawai`ian volcanoes Kohala, 
# Mauna Kea and K\={\i}lauea.
# Or the Mauna Loa Trend: L\=ana`i, Kaho`olawe, the Hawai`ian volcanoes of 
# Hual\=alai and Mauna Loa, as well as Penguin Bank, M\=ahukona and L\=oi`hi
df1 = df.loc[(df['Province']=="Molokai") | (df['Province']=="Maui") | (df['Volcano']=="Kohala") | (df['Volcano']=="Mauna_Kea") | (df['Volcano']=="Kilauea")]
#df1 = df.loc[(df['Province']=="Lanai") | (df['Province']=="Kahoolawe") | (df['Volcano']=="Hualalai") | (df['Volcano']=="Mauna_Loa") | (df['Province']=="Penguin_Bank")]


# ---- Split into Shield, Post-Shield and Rejuvenated Phases
# Here, assume that different phases of volcanism have different source 
# compositions. Unlike in Tutorials 1 and 2, here different samples require 
# different H2O/Ce and src_FeIII_totFe values. So here we add H2O and 
# src_FeIII_totFe columns to the input pandas dataframe.

# Calculate water based on different H2O/Ce values in main text
df1.loc[(df1['Stage']=="Shield"), 'H2O'] = 144. * df1.loc[(df1['Stage']=="Shield"), 'Ce'] / 10000
df1.loc[(df1['Stage']=="Post-Shield"), 'H2O'] = 136. * df1.loc[(df1['Stage']=="Post-Shield"), 'Ce'] / 10000
df1.loc[(df1['Stage']=="Rejuvenated"), 'H2O'] = 211. * df1.loc[(df1['Stage']=="Rejuvenated"), 'Ce'] / 10000

# Assign src_FeIII_totFe values in main text
df1.loc[(df1['Stage']=="Shield"), 'src_FeIII_totFe'] = 0.15
df1.loc[(df1['Stage']=="Post-Shield"), 'src_FeIII_totFe'] = 0.15
df1.loc[(df1['Stage']=="Rejuvenated"), 'src_FeIII_totFe'] = 0.17

# Only include samples that have Ce
df = df1.loc[(df1['Ce']>0)]

# Remove samples not assigned to any phase.
df = df.loc[~df['Stage'].isnull()]

# Save to a csv
df.to_csv("province.csv", sep=',')

# ---- Reading data and initialising the Suite object
s = Suite("province.csv", min_MgO=8.5)

# ---- Backtrack and Estimate pressure and temperature
# See Tutorial 1 for comprehensive explanation
b = BacktrackOlivineFractionation()
s.backtrack_compositions(backtracker=b)
s.compute_pressure_temperature(min_SiO2=40.)

# ---- Fit Tp to Oahu Data
# See Tutorial 1 for comprehensive explanation

# Set up mantle lithology
lz = m.lithologies.katz.lherzolite()
mantle = m.mantle([lz], [1], ['Lz'])
max_P = -lz.parameters['A2'] / (2.*lz.parameters['A3'])
P_sol = np.arange(0., max_P, 0.1)
T_sol = [lz.TSolidus(P) for P in P_sol]

# calculate Tp
s.find_suite_potential_temperature(mantle, find_bounds=True)

# ---- Plot Figure For Oahu
fig, (ax1,ax2) = plt.subplots(1,2)

# Plot Pressure vs. Temperature figure in top left panel

# Plot solidus
ax1.plot(T_sol, P_sol, "k")

# Plot best fitting melt path
ax1.plot(s.path.T, s.path.P, "-", color="k", zorder=1)

# Plot bounding melt paths
ax1.plot(s.upper_path.T, s.upper_path.P, "--", color="k", zorder=1)
ax1.plot(s.lower_path.T, s.lower_path.P, "--", color="k", zorder=1)

# Plot data

ax1.scatter(s.PT['T'][s.data['Stage']=="Post-Shield"], s.PT['P'][s.data['Stage']=="Post-Shield"],
            marker="s", facecolors="deeppink", edgecolor="k", zorder=2)
ax1.scatter(s.PT['T'][s.data['Stage']=="Rejuvenated"], s.PT['P'][s.data['Stage']=="Rejuvenated"],
            marker="o", facecolors="dodgerblue", edgecolor="k", zorder=2)
ax1.scatter(s.PT['T'][s.data['Stage']=="Shield"], s.PT['P'][s.data['Stage']=="Shield"],
            marker="^", facecolors="orange", edgecolor="k", zorder=2)

# Organise axes
ax1.text(0.95, 0.95, "$T_p$ = $%i^{+%i}_{-%i}$ $^\circ$C" % (
   s.potential_temperature, 
   s.upper_potential_temperature - s.potential_temperature,
   s.potential_temperature - s.lower_potential_temperature), verticalalignment='top', 
         horizontalalignment='right', transform=ax1.transAxes, fontsize=12)
ax1.text(0.05, 0.05, 'a', verticalalignment='bottom', 
          horizontalalignment='left', transform=ax1.transAxes, fontsize=12)
ax1.set_xlabel("Temperature [$^\circ$C]")
ax1.set_ylabel("Pressure [GPa]")
ax1.set_xlim((1325.),(1725.))
ax1.set_ylim((0.5),(5.5))
ax1.invert_yaxis()
ax1.set_box_aspect(1)

ax2.scatter(
    s.data.loc[(s.data['143Nd/144Nd']>0.) & (s.data['Nb']>0.) & (s.data['Stage']=="Post-Shield"), '143Nd/144Nd'],
    s.data.loc[(s.data['143Nd/144Nd']>0.) & (s.data['Nb']>0.) & (s.data['Stage']=="Post-Shield"), 'Nb']/
    s.data.loc[(s.data['143Nd/144Nd']>0.) & (s.data['Nb']>0.) & (s.data['Stage']=="Post-Shield"), 'Zr'],
    marker="s", facecolors="deeppink", edgecolor="k", zorder=2
)
ax2.scatter(
    s.data.loc[(s.data['143Nd/144Nd']>0.) & (s.data['Nb']>0.) & (s.data['Stage']=="Rejuvenated"), '143Nd/144Nd'],
    s.data.loc[(s.data['143Nd/144Nd']>0.) & (s.data['Nb']>0.) & (s.data['Stage']=="Rejuvenated"), 'Nb']/
    s.data.loc[(s.data['143Nd/144Nd']>0.) & (s.data['Nb']>0.) & (s.data['Stage']=="Rejuvenated"), 'Zr'],
    marker="o", facecolors="dodgerblue", edgecolor="k", zorder=2
)
ax2.scatter(
    s.data.loc[(s.data['143Nd/144Nd']>0.) & (s.data['Nb']>0.) & (s.data['Stage']=="Shield"), '143Nd/144Nd'],
    s.data.loc[(s.data['143Nd/144Nd']>0.) & (s.data['Nb']>0.) & (s.data['Stage']=="Shield"), 'Nb']/
    s.data.loc[(s.data['143Nd/144Nd']>0.) & (s.data['Nb']>0.) & (s.data['Stage']=="Shield"), 'Zr'],
    marker="^", facecolors="orange", edgecolor="k", zorder=2
)

ax2.set_xlabel("$^{143}$Nd / $^{144}$Nd")
ax2.set_ylabel("Nb / Zr")
ax2.invert_yaxis()
ax2.set_box_aspect(1)

plt.show()