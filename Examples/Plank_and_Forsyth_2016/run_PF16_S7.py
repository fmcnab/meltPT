from meltPT import *
import matplotlib.pyplot as plt
import sys

# ---- check against plank & forsyth supplementary
s = Suite("PF16_S7.csv", src_FeIII_totFe=0.19)
s.backtrack_compositions(Kd=0.3, verbose=True)
s.compute_pressure_temperature()
print()
print("Result from PF16 supplementary 7: P = 2.09 GPa, T = 1347 oC.")
print("Our result:                       P = %.2f GPa, T = %i oC." % 
    (s.PT['P'], s.PT['T']))

# ---- check reading in primaries
s2 = Suite("PF16_S7_Primary.csv", read_as_primary=True)
s2.compute_pressure_temperature()
print("Our result (2):                   P = %.2f GPa, T = %i oC." % 
    (s2.PT['P'], s2.PT['T']))

# ---- Set up mantle
lz = m.lithologies.katz.lherzolite()
mantle = m.mantle([lz], [1], ['Lz'])
P_sol = np.arange(0., 3., 0.1)
T_sol = [lz.TSolidus(P) for P in P_sol]


# ---- Perform fit to melt path
s.find_individual_potential_temperatures(mantle)
print()
print("Best-fitting melting model: Tp = %i oC, F = %.2f %%." % 
    (s.individual_potential_temperatures['Tp'], 
    s.individual_potential_temperatures['F']*100.))


# ---- Plot
melt_label = "Melting path, $T_p = %i ^{\circ}$C" % s.individual_potential_temperatures.iloc[0]['Tp']
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




__, comp = backtrack_sample_composition(s.data.iloc[0], return_all=True)

fig, axs = plt.subplots(2,1,sharex=True)

phases = ['MgO', 'SiO2', 'FeO', 'Al2O3', 'CaO', 'Na2O', 'H2O']
# phases = ['MgO']
count = np.arange(0, len(comp), 1)
ol_added = count*0.0005 / (1.+(count*0.0005))
for p in phases:
    conc = []
    Fo = []
    for i,c in enumerate(comp):
        conc.append( c[p]/comp[0][p] )
        Fo.append(compute_forsterite_number(c))
    axs[0].plot(ol_added, conc, label=p)
axs[0].legend()
axs[0].set_ylabel("Concentration / Initial concentration")

axs[1].plot(ol_added, Fo)
axs[1].set_ylabel("Fo")
axs[1].set_xlabel("Fraction olivine added")
plt.show()



sys.exit()


# ---- Errors

from copy import deepcopy

def water_ox_sweeps(infile):
    
    suite = Suite(infile, src_FeIII_totFe=0.19)

    # define ranges: H2O
    min_H2O = suite.data['H2O'].iloc[0] / 2.
    max_H2O = suite.data['H2O'].iloc[0] * 1.5
    arr_H2O = np.linspace(min_H2O, max_H2O, 11)
    
    # define ranges: Fe3+
    min_Fe = 0.1
    max_Fe = 0.28
    arr_Fe = np.linspace(min_Fe, max_Fe, 11)
    
    # define ranges, Fo
    min_Fo = 0.88
    max_Fo = 0.92
    arr_Fo = np.linspace(min_Fo, max_Fo, 11)
    
    # outputs
    T = np.zeros((len(arr_H2O),len(arr_Fe),len(arr_Fo)))
    P = np.zeros((len(arr_H2O),len(arr_Fe),len(arr_Fo)))
    for i,H2O in enumerate(arr_H2O):
        for j,Fe in enumerate(arr_Fe):
            for k,Fo in enumerate(arr_Fo):
                s = Suite(infile, src_FeIII_totFe=Fe)
                s.data['H2O'] = H2O
                s.backtrack_compositions(Kd=0.3, target_Fo=Fo)
                s.compute_pressure_temperature()
                T[i,j,k] = s.PT['T']
                P[i,j,k] = s.PT['P']
                
    return arr_H2O, arr_Fe, arr_Fo, P, T
        
water, ox, Fo, P, T = water_ox_sweeps("PF16_S7.csv")

fig, axs  = plt.subplots(2,3,sharey="row",sharex="col")

axs[0,0].errorbar(s.data['H2O'].iloc[0], s.PT['P'], yerr=0.24, fmt="o")
axs[0,0].scatter(water, P[:,5,5], marker="+")
axs[0,0].set_ylabel("Pressure [GPa]")

axs[0,1].errorbar(0.19, s.PT['P'], yerr=0.24, fmt="o")
axs[0,1].scatter(ox, P[5,:,5], marker="+")

axs[0,2].errorbar(0.9, s.PT['P'], yerr=0.24, fmt="o")
axs[0,2].scatter(Fo, P[5,5,:], marker="+")

axs[1,0].errorbar(s.data['H2O'].iloc[0], s.PT['T'], yerr=39., fmt="o")
axs[1,0].scatter(water, T[:,5,5], marker="+")
axs[1,0].set_ylabel(r"Temperature [$^{\circ}$C]")
axs[1,0].set_xlabel(r"H$_2$O [wt%]")

axs[1,1].errorbar(0.19, s.PT['T'], yerr=39., fmt="o")
axs[1,1].scatter(ox, T[5,:,5], marker="+")
axs[1,1].set_xlabel(r"Fe$^{3+}$ / $\Sigma$Fe")

axs[1,2].errorbar(0.9, s.PT['T'], yerr=39., fmt="o")
axs[1,2].scatter(Fo, T[5,5,:], marker="+")
axs[1,2].set_xlabel("Fo")

for axss in axs:
    for ax in axss:
        ax.set_aspect(1./ax.get_data_ratio())

plt.show()

def plot_sweeps(grid, label):

    fig, axs = plt.subplots(2,2)

    cpt_00 = axs[0,0].contourf(ox, water, grid[:,:,5])
    # ax1.contour(ox, water, P, [2.], linestyles="dashed")
    cbar = fig.colorbar(cpt_00, ax=axs[0,0])
    cbar.set_label(label)
    axs[0,0].set_ylabel(r"H$_2$O [wt%]")
    axs[0,0].set_xlabel(r"Fe$^{3+}$ / $\Sigma$Fe")

    cpt_01 = axs[0,1].contourf(Fo, water, grid[:,5,:])
    cbar = fig.colorbar(cpt_01, ax=axs[0,1])
    cbar.set_label(label)

    cpt_11 = axs[1,1].contourf(Fo, ox, P[5,:,:])
    cbar = fig.colorbar(cpt_11, ax=axs[1,1])
    cbar.set_label(label)
    axs[1,1].set_ylabel(r"Fe$^{3+}$ / $\Sigma$Fe")
    axs[1,1].set_xlabel(r"Fo")

    for axss in axs:
        for ax in axss:
            ax.set_aspect(1./ax.get_data_ratio())
    axs[1,0].axis("off")
    fig.tight_layout()
    plt.show()

plot_sweeps(P, "Pressure [GPa]")
plot_sweeps(T, "Temperature [$^{\circ}$C]")

plt.plot(T_sol, P_sol, "k", label="Solidus")
plt.errorbar(s.PT['T'], s.PT['P'], xerr=39., yerr=0.24, fmt="o")
for t,p in zip(T,P):
    plt.plot(t, p, "+")
plt.gca().invert_yaxis()
plt.show()








# 
# import random 
# 
# def monte_carlo(infile):
# 
#     suite = Suite(infile, src_FeIII_totFe=0.19)
# 
#     # define ranges: H2O
#     min_H2O = suite.data['H2O'].iloc[0] - (0.5 * suite.data['H2O'].iloc[0])
#     max_H2O = suite.data['H2O'].iloc[0] + (0.5 * suite.data['H2O'].iloc[0])
# 
#     # define ranges: Fe3+
#     min_Fe = 0.1
#     max_Fe = 0.28
# 
#     # outputs
#     n = 200
#     H2O = np.zeros(n)
#     Fe = np.zeros(n)
#     T = np.zeros(n)
#     P = np.zeros(n)
#     for i in range(n):
#         H2O = random.uniform(min_H2O, max_H2O)
#         Fe = random.uniform(min_Fe, max_Fe)
#         s = Suite(infile, src_FeIII_totFe=Fe)
#         s.data['H2O'] = H2O
#         s.backtrack_compositions(Kd=0.3)
#         s.compute_pressure_temperature()
#         T[i] = s.PT['T']
#         P[i] = s.PT['P']
# 
#     return H2O, Fe, T, P
# 
# water, ox, T, P = monte_carlo("PF16_S7.csv")
# 
# plt.plot(T_sol, P_sol, "k", label="Solidus")
# plt.errorbar(s.PT['T'], s.PT['P'], xerr=39., yerr=0.24, fmt="o")
# plt.plot(T, P, "+")
# plt.gca().invert_yaxis()
# plt.show()
# 
# def gauss(x, mean, sd):
#     return (1./(sd*np.sqrt(2.*np.pi))) * np.exp(-((x-mean)/(np.sqrt(2)*sd))**2.)
# 
# P_rng = np.linspace(0.,4.,1000)
# pdf = np.zeros(len(P_rng))
# for i,p in enumerate(P):
#     if not np.isnan(p):
#         pdf += gauss(P_rng,p,0.24)
# plt.plot(P_rng, pdf/pdf.max())
# orig = gauss(P_rng, s.PT['P'].iloc[0], 0.24)
# plt.plot(P_rng, orig/orig.max(), "--")
# plt.show()
# 
# T_rng = np.linspace(1000.,1900.,1000)
# pdf = np.zeros(len(T_rng))
# for i,t in enumerate(T):
#     if not np.isnan(t):
#         pdf += gauss(T_rng,t,39.)
# plt.plot(T_rng, pdf/pdf.max())
# orig = gauss(T_rng, s.PT['T'].iloc[0], 39.)
# plt.plot(T_rng, orig/orig.max(), "--")
# plt.show()
