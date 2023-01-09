#------------------------------------------------------------------------------
# ---- Figure S2-S5
#
# Recreate elements of Figure S2-S5 from McNab & Ball (2023).
# Comment/uncomment lines 123-125/232-235 to include/exclude thermobarometer error.
# This will only work with v1.1.0!
#------------------------------------------------------------------------------


from meltPT import *
import pyMelt as m
import matplotlib.pyplot as plt
import random
import sys

class Parameter:
    def __init__(self, val, plus, minus):
        self.val = val
        self.plus = plus
        self.minus = minus
    
    def draw_normal(self):
        return random.gauss(self.val, (self.plus+self.minus)/2.)
        
    def range(self, n=10):
        return np.linspace(self.val-(3.*self.minus), self.val+(3.*self.plus), n)
        
    def gaussian(self):
        sd = (self.plus+self.minus)/2.
        return (
            (1./(sd*np.sqrt(2.*np.pi))) * 
            np.exp(-0.5*(((self.range(n=100)-self.val)/sd)**2.))
            )

    def fechner(self):
        x = self.range(n=100)
        upper = np.where(x >= self.val)
        f = (2.*self.minus)/(self.plus + self.minus)/self.minus
        f *= 1./(np.sqrt(2.*np.pi))*np.exp(-0.5*(((x-self.val)/self.minus)**2.))
        f[upper] = (2.*self.plus)/(self.plus + self.minus)/self.plus
        f[upper] *= 1./(np.sqrt(2.*np.pi))*np.exp(-0.5*(((x[upper]-self.val)/self.plus)**2.))
        return f

# ---- Values
shield_src_FeIII_totFe = Parameter(0.15, 0.02, 0.02)
shield_Ce_to_H2O = Parameter(144, 56, 56)
rej_src_FeIII_totFe = Parameter(0.17, 0.04, 0.04)
rej_Ce_to_H2O = Parameter(211, 29, 29)
src_Fo = Parameter(0.9, 0.005, 0.005)

# ---- Isolate Island Data
# Read in compilation dataset of Hawaiian samples from McNab & Ball, (2022)
df = pd.read_csv("Dataset_S1.csv", sep=',')

# We want take only samples from our Hawaii database that correspond to
# the island of Oahu. The other islands are shown below for users to try
# Island options: Niihau / Kaula / Kauai / Oahu / Molokai / Maui
#                 Kahoolawe / Lanai / Hawaii
Island = "Oahu"
df = df.loc[(df['Province']==Island)]
df = df.loc[(df['Ce']>0)]
df = df.loc[~df['Stage'].isnull()]

df.loc[(df['Stage']=="Rejuvenated"), 'src_FeIII_totFe'] = rej_src_FeIII_totFe.val
df.loc[(df['Stage']=="Shield"), 'src_FeIII_totFe'] = shield_src_FeIII_totFe.val

df.loc[(df['Stage']=="Rejuvenated"), 'H2O'] = df.loc[(df['Stage']=="Rejuvenated"), 'Ce'] * rej_Ce_to_H2O.val / 10000.
df.loc[(df['Stage']=="Shield"), 'H2O'] = df.loc[(df['Stage']=="Shield"), 'Ce'] * shield_Ce_to_H2O.val / 10000.

df.to_csv("./all.csv", sep=',')



# ---- Standard way
s = Suite(
    "./all.csv",
    min_MgO=8.5
    )

b = BacktrackOlivineFractionation()
s.backtrack_compositions(backtracker=b)
s.compute_pressure_temperature(min_SiO2=40.)
lz = m.lithologies.katz.lherzolite()
mantle = m.mantle([lz], [1], ['Lz'])
s.find_suite_potential_temperature(mantle, find_bounds=True)


# ---- MC
ss = []

for i in range(10):
    
    print(i)
    
    si = Suite("./all.csv", min_MgO=8.5, src_FeIII_totFe=0.15, src_Fo=0.9, Ce_to_H2O=144.)

    for j in si.data.index:
        
        if si.data.at[j,'Stage'] == "Shield":
            si.data.at[j,'src_FeIII_totFe'] = shield_src_FeIII_totFe.draw_normal()
            ratio = shield_Ce_to_H2O.draw_normal()
            while ratio < 0.:
                ratio = shield_Ce_to_H2O.draw_normal()
        elif si.data.at[j,'Stage'] == "Rejuvenated":
            si.data.at[j,'src_FeIII_totFe'] = rej_src_FeIII_totFe.draw_normal()
            ratio = rej_Ce_to_H2O.draw_normal()
            while ratio < 0.:
                ratio = rej_Ce_to_H2O.draw_normal()

        
        si.data.at[j,'FeO_tot'] = si.data.at[j,'FeO'] + (si.data.at[j,'Fe2O3'] * (74.84 * 2.) / 159.69)
        si.data.at[j,'FeO'] = si.data.at[j,'FeO_tot'] * (1. - si.data.at[j,'src_FeIII_totFe'])
        si.data.at[j,'Fe2O3'] = si.data.at[j, 'FeO_tot'] * si.data.at[j,'src_FeIII_totFe'] * 159.69 / 71.84 / 2.      

        si.data.at[j,'H2O'] = si.data.at[j,'Ce'] * ratio / 10000.
        
        si.data.at[j, 'src_Fo'] = src_Fo.draw_normal()
    
    bi = BacktrackOlivineFractionation()
    si.backtrack_compositions(backtracker=b)
    si.compute_pressure_temperature(min_SiO2=40.)
    
    # for j in s.data.index:    
    #     si.PT.at[j,'P'] = Parameter(si.PT.at[j,'P'], si.PT.at[j,'P_err'], si.PT.at[j,'P_err']).draw_normal()
    #     si.PT.at[j,'T'] = Parameter(si.PT.at[j,'T'], si.PT.at[j,'T_err'], si.PT.at[j,'T_err']).draw_normal()

    si.find_suite_potential_temperature(mantle, find_bounds=True)

    ss.append(si)

# ---- Plot

fig, axs = plt.subplots(2,3)

# Iteration 0: source Fo#
axs[0,0].hist(ss[0].data.loc[(ss[0].data['Stage']=="Shield"), 'src_Fo']*100.)
axs[0,0].hist(ss[0].data.loc[(ss[0].data['Stage']=="Rejuvenated"), 'src_Fo']*100.)
axs[0,0].set_ylabel("Count")
axs[0,0].set_xlabel("Fo#")

# Iteration 0: source Fe3+/totFe
axs[0,1].hist(ss[0].data.loc[(ss[0].data['Stage']=="Shield"), 'src_FeIII_totFe'])
axs[0,1].hist(ss[0].data.loc[(ss[0].data['Stage']=="Rejuvenated"), 'src_FeIII_totFe'])
axs[0,1].set_xlabel(r"Fe$^{3+}$/$\Sigma$Fe")

# Iteration 0: H2O/Ce
axs[0,2].hist(
    ss[0].data.loc[(ss[0].data['Stage']=="Shield"), 'H2O']/
    ss[0].data.loc[(ss[0].data['Stage']=="Shield"), 'Ce']*10000.)
axs[0,2].hist(
    ss[0].data.loc[(ss[0].data['Stage']=="Rejuvenated"), 'H2O']/
    ss[0].data.loc[(ss[0].data['Stage']=="Rejuvenated"), 'Ce']*10000.)
axs[0,2].set_xlabel(r"H$_2$O / Ce")

# Iteration 0: PT
max_P = -lz.parameters['A2'] / (2.*lz.parameters['A3'])
P_sol = np.arange(0., max_P, 0.1)
T_sol = [lz.TSolidus(P) for P in P_sol]
axs[1,0].plot(T_sol, P_sol, "k")
axs[1,0].plot(ss[0].path.T, ss[0].path.P, "-", color="k", zorder=1)
axs[1,0].scatter(ss[0].PT['T'][ss[0].data['Stage']=="Shield"], ss[0].PT['P'][ss[0].data['Stage']=="Shield"],
            marker="^", facecolors="orange", edgecolor="k", zorder=2)
axs[1,0].scatter(ss[0].PT['T'][ss[0].data['Stage']=="Post-Shield"], ss[0].PT['P'][ss[0].data['Stage']=="Post-Shield"],
            marker="s", facecolors="deeppink", edgecolor="k", zorder=2)
axs[1,0].scatter(ss[0].PT['T'][ss[0].data['Stage']=="Rejuvenated"], ss[0].PT['P'][ss[0].data['Stage']=="Rejuvenated"],
            marker="o", facecolors="dodgerblue", edgecolor="k", zorder=2)
axs[1,0].set_xlabel("Temperature [$^\circ$C]")
axs[1,0].set_ylabel("Pressure [GPa]")
axs[1,0].set_xlim((1325.),(1725.))
axs[1,0].set_ylim((0.5),(5.5))
axs[1,0].invert_yaxis()

# Reference: PT
axs[1,1].plot(T_sol, P_sol, "k")
axs[1,1].plot(s.path.T, s.path.P, "-", color="k", zorder=1)
axs[1,1].scatter(s.PT['T'][s.data['Stage']=="Shield"], s.PT['P'][s.data['Stage']=="Shield"],
            marker="^", facecolors="orange", edgecolor="k", zorder=2)
axs[1,1].scatter(s.PT['T'][s.data['Stage']=="Post-Shield"], s.PT['P'][s.data['Stage']=="Post-Shield"],
            marker="s", facecolors="deeppink", edgecolor="k", zorder=2)
axs[1,1].scatter(s.PT['T'][s.data['Stage']=="Rejuvenated"], s.PT['P'][s.data['Stage']=="Rejuvenated"],
            marker="o", facecolors="dodgerblue", edgecolor="k", zorder=2)
axs[1,1].set_xlabel("Temperature [$^\circ$C]")
axs[1,1].set_xlim((1325.),(1725.))
axs[1,1].set_ylim((0.5),(5.5))
axs[1,1].invert_yaxis()

# Probabilities
Tp = Parameter(
    s.potential_temperature,
    s.upper_potential_temperature-s.potential_temperature,
    s.potential_temperature-s.lower_potential_temperature)
axs[1,2].plot(Tp.range(100), Tp.fechner(), linewidth=2)
for si in ss:
    Tp = Parameter(
        si.potential_temperature,
        si.upper_potential_temperature-si.potential_temperature,
        si.potential_temperature-si.lower_potential_temperature)
    axs[1,2].plot(Tp.range(100), Tp.fechner(), linewidth=0.5)
axs[1,2].set_xlabel(r"$T_p$ [$^\circ$C]")
axs[1,2].set_ylabel("Probability")

# Show
plt.show()



# ---- MC2
ss2 = []

for i in range(10):
    
    print(i)
    
    si = Suite("./all.csv", min_MgO=8.5, src_FeIII_totFe=0.15, src_Fo=0.9, Ce_to_H2O=144.)

    si.data.loc[(si.data['Stage']=="Rejuvenated"), 'src_FeIII_totFe'] = rej_src_FeIII_totFe.draw_normal()
    si.data.loc[(si.data['Stage']=="Shield"), 'src_FeIII_totFe'] = shield_src_FeIII_totFe.draw_normal()
    si.data['FeO_tot'] = si.data['FeO'] + (si.data['Fe2O3'] * (74.84 * 2.) / 159.69)
    si.data['FeO'] = si.data['FeO_tot'] * (1. - si.data['src_FeIII_totFe'])
    si.data['Fe2O3'] = si.data['FeO_tot'] * si.data['src_FeIII_totFe'] * 159.69 / 71.84 / 2.      

    si.data.loc[(si.data['Stage']=="Rejuvenated"), 'H2O'] = si.data.loc[(si.data['Stage']=="Rejuvenated"), 'Ce'] * rej_Ce_to_H2O.draw_normal() / 10000.
    si.data.loc[(si.data['Stage']=="Shield"), 'H2O'] = si.data.loc[(si.data['Stage']=="Shield"), 'Ce'] * shield_Ce_to_H2O.draw_normal() / 10000.
    
    si.data.loc[(si.data['Stage']=="Rejuvenated"), 'src_Fo'] = src_Fo.draw_normal()
    si.data.loc[(si.data['Stage']=="Shield"), 'src_Fo'] = src_Fo.draw_normal()
    
    bi = BacktrackOlivineFractionation()
    si.backtrack_compositions(backtracker=b)
    si.compute_pressure_temperature(min_SiO2=40.)
    
    # si.PT.loc[(si.data['Stage']=="Rejuvenated"),'P'] += Value(0., 0.24, 0.24).draw_normal()
    # si.PT.loc[(si.data['Stage']=="Rejuvenated"),'T'] += Value(0., 39., 39.).draw_normal()
    # si.PT.loc[(si.data['Stage']=="Shield"),'P'] += Value(0., 0.24, 0.24).draw_normal()
    # si.PT.loc[(si.data['Stage']=="Shield"),'T'] += Value(0., 39., 39.).draw_normal()

    si.find_suite_potential_temperature(mantle, find_bounds=True)

    ss2.append(si)

# ---- Plot

fig, axs = plt.subplots(2,3)

# Source Fo#
axs[0,0].hist([si.data.loc[(si.data['Stage']=="Shield"), 'src_Fo'].iloc[0]*100. for si in ss2])
axs[0,0].hist([si.data.loc[(si.data['Stage']=="Rejuvenated"), 'src_Fo'].iloc[0]*100. for si in ss2])
axs[0,0].set_ylabel("Count")
axs[0,0].set_xlabel("Fo#")

# Source Fe3+/totFe
axs[0,1].hist([si.data.loc[(si.data['Stage']=="Shield"), 'src_FeIII_totFe'].iloc[0] for si in ss2])
axs[0,1].hist([si.data.loc[(si.data['Stage']=="Rejuvenated"), 'src_FeIII_totFe'].iloc[0] for si in ss2])
axs[0,1].set_xlabel(r"Fe$^{3+}$/$\Sigma$Fe")

# H2O/Ce
axs[0,2].hist(
    [si.data.loc[(si.data['Stage']=="Shield"), 'H2O'].iloc[0]/
     si.data.loc[(si.data['Stage']=="Shield"), 'Ce'].iloc[0]*10000. for si in ss2])
axs[0,2].hist(
    [si.data.loc[(si.data['Stage']=="Rejuvenated"), 'H2O'].iloc[0]/
     si.data.loc[(si.data['Stage']=="Rejuvenated"), 'Ce'].iloc[0]*10000. for si in ss2])
axs[0,2].set_xlabel(r"H$_2$O / Ce")

# Iteration 0: PT
max_P = -lz.parameters['A2'] / (2.*lz.parameters['A3'])
P_sol = np.arange(0., max_P, 0.1)
T_sol = [lz.TSolidus(P) for P in P_sol]
axs[1,0].plot(T_sol, P_sol, "k")
axs[1,0].plot(ss2[0].path.T, ss2[0].path.P, "-", color="k", zorder=1)
axs[1,0].scatter(ss2[0].PT['T'][ss2[0].data['Stage']=="Shield"], ss2[0].PT['P'][ss2[0].data['Stage']=="Shield"],
            marker="^", facecolors="orange", edgecolor="k", zorder=2)
axs[1,0].scatter(ss2[0].PT['T'][ss2[0].data['Stage']=="Post-Shield"], ss2[0].PT['P'][ss2[0].data['Stage']=="Post-Shield"],
            marker="s", facecolors="deeppink", edgecolor="k", zorder=2)
axs[1,0].scatter(ss2[0].PT['T'][ss2[0].data['Stage']=="Rejuvenated"], ss2[0].PT['P'][ss2[0].data['Stage']=="Rejuvenated"],
            marker="o", facecolors="dodgerblue", edgecolor="k", zorder=2)
axs[1,0].set_xlabel("Temperature [$^\circ$C]")
axs[1,0].set_ylabel("Pressure [GPa]")
axs[1,0].set_xlim((1325.),(1725.))
axs[1,0].set_ylim((0.5),(5.5))
axs[1,0].invert_yaxis()

# Reference: PT
axs[1,1].plot(T_sol, P_sol, "k")
axs[1,1].plot(s.path.T, s.path.P, "-", color="k", zorder=1)
axs[1,1].scatter(s.PT['T'][s.data['Stage']=="Shield"], s.PT['P'][s.data['Stage']=="Shield"],
            marker="^", facecolors="orange", edgecolor="k", zorder=2)
axs[1,1].scatter(s.PT['T'][s.data['Stage']=="Post-Shield"], s.PT['P'][s.data['Stage']=="Post-Shield"],
            marker="s", facecolors="deeppink", edgecolor="k", zorder=2)
axs[1,1].scatter(s.PT['T'][s.data['Stage']=="Rejuvenated"], s.PT['P'][s.data['Stage']=="Rejuvenated"],
            marker="o", facecolors="dodgerblue", edgecolor="k", zorder=2)
axs[1,1].set_xlabel("Temperature [$^\circ$C]")
axs[1,1].set_xlim((1325.),(1725.))
axs[1,1].set_ylim((0.5),(5.5))
axs[1,1].invert_yaxis()

# Probabilities
Tp = Parameter(
    s.potential_temperature,
    s.upper_potential_temperature-s.potential_temperature,
    s.potential_temperature-s.lower_potential_temperature)
axs[1,2].plot(Tp.range(100), Tp.fechner(), linewidth=2)
for si in ss2:
    Tp = Parameter(
        si.potential_temperature,
        si.upper_potential_temperature-si.potential_temperature,
        si.potential_temperature-si.lower_potential_temperature)
    axs[1,2].plot(Tp.range(100), Tp.fechner(), linewidth=0.5)
axs[1,2].set_xlabel(r"$T_p$ [$^\circ$C]")
axs[1,2].set_ylabel("Probability")

# Show
plt.show()