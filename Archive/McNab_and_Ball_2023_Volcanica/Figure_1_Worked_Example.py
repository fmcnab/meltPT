from meltPT import *
import matplotlib.pyplot as plt
import sys

# ---- INTRODUCTION
# Recreate elements of Figure 1 from McNab & Ball (2023).
# See also Tutorial 1 for more details and a version that works with the
# current version of meltPT.
# This will only work with v1.1.0!

# ---- Load the data, do the backtracking, compute pressures and temperatures.
# ---- First we will use a fixed value for Kd.
s1 = Suite("UT09DV04.csv", src_FeIII_totFe=0.17)
b1 = BacktrackOlivineFractionation(Kd=0.3)
s1.backtrack_compositions(b1)
s1.compute_pressure_temperature()


# ---- Get the intermediary compositions using the "return all" flag.
__, comp1 = b1.backtrack_sample_composition(s1.data.iloc[0], return_all=True)

# ---- Repeat with variable Kd
s2 = Suite("UT09DV04.csv", src_FeIII_totFe=0.17)
b2 = BacktrackOlivineFractionation()
s2.backtrack_compositions(b2)
s2.compute_pressure_temperature()
__, comp2 = b2.backtrack_sample_composition(s2.data.iloc[0], return_all=True)


# ---- Compute amounts of olivine added
count1 = np.arange(0, len(comp1), 1)        
ol_added1 = count1*b1.dm / (1.+(count1*b1.dm)) * 100.
count2 = np.arange(0, len(comp2), 1)        
ol_added2 = count2*b2.dm / (1.+(count2*b2.dm)) * 100.


# ---- Compute normalised concentrations through time
phases = ['MgO', 'SiO2', 'FeO', 'Al2O3', 'CaO', 'Na2O', 'H2O']
concs1 = {}
concs2 = {}
for p in phases:
    concs1[p] = []
    concs2[p] = []
    for i,c in enumerate(comp1):
        concs1[p].append( c[p]/comp1[0][p] )
    for i,c in enumerate(comp2):
        concs2[p].append( c[p]/comp2[0][p] )


# ---- Make backtracking plots.
fig, axs = plt.subplots(1,2,sharex=True)
axs[0].plot(ol_added1, [c['Fo'] for c in comp1])
axs[0].plot(ol_added2, [c['Fo'] for c in comp2], "--")
axs[0].set_xlabel("Olivine added [%]")
axs[0].set_ylabel("Fo#")
axs[0].set_box_aspect(1)
labels = ["MgO", "SiO2", "FeO", "", "", "", "Other"]
colors = ["r", "g", "b", "grey", "grey", "grey", "grey"]
for i,p in enumerate(phases):
    axs[1].plot(ol_added1, concs1[p], color=colors[i], label=labels[i])
    axs[1].plot(ol_added2, concs2[p], "--", color=colors[i])
axs[1].set_xlabel("Olivine added [%]")
axs[1].set_ylabel(r"$C$ / $C_0$")
axs[1].set_box_aspect(1)
plt.legend()
plt.show()

    
# ---- Compute potential temperatures
lz = m.lithologies.katz.lherzolite()
mantle = m.mantle([lz], [1], ['Lz'])
s1.find_individual_potential_temperatures(mantle)
s2.find_individual_potential_temperatures(mantle)


# ---- Plot pressure/temperature
P_sol = np.arange(1., 4., 0.1)
T_sol = lz.TSolidus(P_sol)
plt.plot(T_sol, P_sol, "grey")
plt.plot(
    s1.individual_potential_temperatures['path'].iloc[0].T,
    s1.individual_potential_temperatures['path'].iloc[0].P)
plt.plot(
    s2.individual_potential_temperatures['path'].iloc[0].T,
    s2.individual_potential_temperatures['path'].iloc[0].P)
lab = r"$K_d = 0.3$"
plt.scatter(s1.PT['T'], s1.PT['P'], label=lab)
lab = r"Variable $K_d$"
plt.scatter(s2.PT['T'], s2.PT['P'], label=lab)
plt.legend()
plt.ylim(1.,4.)
plt.xlim(1300.,1500.)
plt.xlabel(r"Temperature ($^{\circ}$C)")
plt.ylabel("Pressure (GPa)")
plt.gca().invert_yaxis()
plt.show()

# ---- Try out some different schemes

# set up
s = Suite("UT09DV04.csv", src_FeIII_totFe=0.17)
b = BacktrackOlivineFractionation(Kd=0.3)
s.backtrack_compositions(b)

# start the plot
P_sol = np.arange(1., 4.5, 0.1)
T_sol = lz.TSolidus(P_sol)
plt.plot(T_sol, P_sol, "grey")

# First let's define a list of the schemes we would like to use. We will take
# a subset for now.
schemes =  ["PF16", "L09", "P08", "BK21", "SD20"]
schemesT = ["B93", "P07_2", "P07_4", "HA15"]

# Loop over the schemes
for scheme in schemes:
    s.compute_pressure_temperature(method=scheme)
    plt.errorbar(
        s.PT['T'], s.PT['P'], 
        xerr=s.PT['T_err'], yerr=s.PT['P_err'],
        label=scheme,
        marker="o")
for scheme in schemesT:
    s.compute_temperature(method=scheme, P=s1.PT['P'])
    plt.errorbar(
        s.PT['T'], s1.PT['P'], 
        xerr=s.PT['T_err'],
        label=scheme,
        marker="^")
        
# finish up
plt.ylim(1.,4.5)
plt.xlim(1300.,1550.)
plt.xlabel(r"Temperature ($^{\circ}$C)")
plt.ylabel("Pressure (GPa)")
plt.legend()
plt.gca().invert_yaxis()
plt.show()

