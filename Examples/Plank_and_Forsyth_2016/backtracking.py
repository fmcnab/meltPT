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


# ---- Brute force backtracking
__, comp = backtrack_sample_composition(s.data.iloc[0], Kd=0.3, verbose=True, return_all=True)

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
        Fo.append(compute_forsterite_number(c, Kd=0.3))
    axs[0].plot(ol_added, conc, label=p)
axs[0].legend()
axs[0].set_ylabel("Concentration / Initial concentration")

axs[1].plot(ol_added, Fo)
axs[1].set_ylabel("Fo")
axs[1].set_xlabel("Fraction olivine added")
plt.show()

# ---- Save

with open("backtrack_concentrations.dat", "wb") as f:
    cols = ["blue", "0/125/0", "red", "darkgrey"]
    for i,phase in enumerate(['MgO', 'FeO', 'SiO2', 'Al2O3']):
        hdr = b"> -W1p,%s\n" % (cols[i].encode("ascii"))
        f.write(hdr)
        arr = np.column_stack(( ol_added, [c[phase]/comp[0][phase] for c in comp] ))
        np.savetxt(f, arr)

with open("labels.dat", "wb") as f:
    label = ["MgO", "FeO", "SiO@-2@-", "Other"]
    for i,phase in enumerate(['MgO', 'FeO', 'SiO2', 'Al2O3']):
        hdr = b"%f %f 8p,Helvetica,%s %s\n" % (
            ol_added[-1],
            comp[-1][phase]/comp[0][phase],
            cols[i].encode("ascii"),
            label[i].encode("ascii")
            )
        f.write(hdr)

with open("forsterite.dat", "wb") as f:
    arr = np.column_stack(( ol_added, Fo ))
    np.savetxt(f, arr)
    
    
    

# ---- Brute force backtracking
__, comp = backtrack_sample_composition(s.data.iloc[0], verbose=True, return_all=True)

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

# ---- Save

with open("backtrack_concentrations_varKd.dat", "wb") as f:
    cols = ["blue", "0/125/0", "red", "darkgrey"]
    for i,phase in enumerate(['MgO', 'FeO', 'SiO2', 'Al2O3']):
        hdr = b"> -W1p,%s,-\n" % (cols[i].encode("ascii"))
        f.write(hdr)
        arr = np.column_stack(( ol_added, [c[phase]/comp[0][phase] for c in comp] ))
        np.savetxt(f, arr)

with open("labels_varKd.dat", "wb") as f:
    label = ["MgO", "FeO", "SiO@-2@-", "Other"]
    for i,phase in enumerate(['MgO', 'FeO', 'SiO2', 'Al2O3']):
        hdr = b"%f %f 8p,Helvetica,%s %s\n" % (
            ol_added[-1],
            comp[-1][phase]/comp[0][phase],
            cols[i].encode("ascii"),
            label[i].encode("ascii")
            )
        f.write(hdr)

with open("forsterite_varKd.dat", "wb") as f:
    arr = np.column_stack(( ol_added, Fo ))
    np.savetxt(f, arr)


# ---- melt path etc.

lz = m.lithologies.katz.lherzolite()
mantle = m.mantle([lz], [1], ['Lz'])
P_sol = np.arange(0., 4.1, 0.1)
T_sol = [lz.TSolidus(P) for P in P_sol]

s.find_individual_potential_temperatures(mantle)
print()
print("Best-fitting melting model: Tp = %i oC, F = %.2f %%." % 
    (s.individual_potential_temperatures['Tp'], 
    s.individual_potential_temperatures['F']*100.))

with open("temperature_pressure.dat", "wb") as f:
    arr = np.column_stack(( s.PT['T'], s.PT['P'] ))
    np.savetxt(f, arr)
    
with open("solidus.dat", "wb") as f:
    arr = np.column_stack(( T_sol, P_sol ))
    np.savetxt(f, arr)
    
with open("melt_path.dat", "wb") as f:
    sub_sol_P = np.arange(P_sol.max(), mantle.solidusIntersection(s.individual_potential_temperatures.iloc[0]['Tp']), -0.01)
    arr = np.column_stack((
        mantle.adiabat(sub_sol_P, s.individual_potential_temperatures.iloc[0]['Tp']),
        sub_sol_P
        ))
    np.savetxt(f, arr)
    arr = np.column_stack((
        s.individual_potential_temperatures.iloc[0]['path'].T,
        s.individual_potential_temperatures.iloc[0]['path'].P
        ))
    np.savetxt(f, arr)
    
with open("adiabat.dat", "wb") as f:
    supra_sol_P = np.arange(mantle.solidusIntersection(s.individual_potential_temperatures.iloc[0]['Tp']), 0., -0.01)
    arr = np.column_stack((
        mantle.adiabat(supra_sol_P, s.individual_potential_temperatures.iloc[0]['Tp']),
        supra_sol_P
        ))
    np.savetxt(f, arr)
    
    
    
# # ---- other geotherm
# 
# # def calc_geotherm(Tp, z_l, mantle):
# 
# from scipy.interpolate import interp1d, splrep, splev
# 
# def lith_geotherm(P, P_l, mantle, Tp):
#     T_l = mantle.adiabat(P_l, Tp)
#     return (T_l / P_l) * P
# 
# P_max = 4.
# Tp = 1330.
# P_l = 1.5
# P = np.arange(0., P_max, 0.001)
# init_lithT = lith_geotherm(P, P_l, mantle, Tp)
# init_asthT = mantle.adiabat(P, Tp)
# lithT = init_lithT[ np.where(init_lithT < mantle.lithologies[0].TSolidus(P)) ]
# lithP = P[ np.where(init_lithT < mantle.lithologies[0].TSolidus(P)) ]
# asthT = init_asthT[ np.where(init_asthT < mantle.lithologies[0].TSolidus(P)) ]
# asthP = P[ np.where(init_asthT < mantle.lithologies[0].TSolidus(P)) ]
# 
# geotherm = interp1d(
#     np.hstack(( lithP, asthP )),
#     np.hstack(( lithT, asthT )),
#     kind="quadratic"
#     )
# 
# # spline = splrep(
# #     np.hstack(( lithP, asthP )),
# #     np.hstack(( lithT, asthT )),
# #     k=2,
# #     s=1000000.
# #     )
# # geotherm = splev(P, spline)
# 
# T = np.arange(1., asthT.max()-1., 1.)
# plt.plot(T_sol, P_sol)
# plt.plot(lithT, lithP, "--")
# plt.plot(asthT, asthP, "--")
# plt.plot(geotherm(P), P, ":")
# plt.gca().invert_yaxis()
# plt.show()