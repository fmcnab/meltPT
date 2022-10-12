from meltPT import *
import matplotlib.pyplot as plt
import sys

# ---- check against plank & forsyth supplementary
# s = Suite("PF16_S7.csv", src_FeIII_totFe=0.19)
s = Suite("UT09DV04.csv", src_FeIII_totFe=0.17)
s.backtrack_compositions(Kd=0.3)
s.compute_pressure_temperature()
print()
print("Result from PF16 supplementary 7: P = 2.09 GPa, T = 1347 oC.")
print("Our result:                       P = %.2f GPa, T = %i oC." % 
    (s.PT['P'], s.PT['T']))


# ---- Brute force backtracking
__, comp = backtrack_sample_composition(s.data.iloc[0], Kd=0.3, verbose=True, return_all=True)

# fig, axs = plt.subplots(2,1,sharex=True)

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
#     axs[0].plot(ol_added, conc, label=p)
# axs[0].legend()
# axs[0].set_ylabel("Concentration / Initial concentration")
# 
# axs[1].plot(ol_added, Fo)
# axs[1].set_ylabel("Fo")
# axs[1].set_xlabel("Fraction olivine added")
# plt.show()

# ---- Save

out_dir = "../../../meltPT_working/"

with open(out_dir + "backtrack_concentrations.dat", "wb") as f:
    cols = ["blue", "0/125/0", "red", "darkgrey"]
    for i,phase in enumerate(['MgO', 'FeO', 'SiO2', 'Al2O3']):
        hdr = b"> -W0.7p,%s\n" % (cols[i].encode("ascii"))
        f.write(hdr)
        arr = np.column_stack(( ol_added*100., [c[phase]/comp[0][phase] for c in comp] ))
        np.savetxt(f, arr)

with open(out_dir + "labels.dat", "wb") as f:
    label = ["MgO", "FeO", "SiO@-2@-", "Other"]
    for i,phase in enumerate(['MgO', 'FeO', 'SiO2', 'Al2O3']):
        hdr = b"%f %f 7p,Helvetica,%s %s\n" % (
            ol_added[-1]*100.,
            comp[-1][phase]/comp[0][phase],
            cols[i].encode("ascii"),
            label[i].encode("ascii")
            )
        f.write(hdr)

with open(out_dir + "forsterite.dat", "wb") as f:
    arr = np.column_stack(( ol_added*100., np.array(Fo)*100. ))
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

with open(out_dir + "temperature_pressure.dat", "wb") as f:
    arr = np.column_stack(( s.PT['T'], s.PT['P'] ))
    np.savetxt(f, arr)

with open(out_dir + "solidus.dat", "wb") as f:
    arr = np.column_stack(( T_sol, P_sol ))
    np.savetxt(f, arr)

with open(out_dir + "melt_path.dat", "wb") as f:
    hdr = b"> -L%i\n" % (s.individual_potential_temperatures.iloc[0]['Tp'])
    f.write(hdr)
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
    
with open(out_dir + "adiabat.dat", "wb") as f:
    hdr = b"> -L%i\n" % (s.individual_potential_temperatures.iloc[0]['Tp'])
    f.write(hdr)
    supra_sol_P = np.arange(mantle.solidusIntersection(s.individual_potential_temperatures.iloc[0]['Tp']), 0., -0.01)
    arr = np.column_stack((
        mantle.adiabat(supra_sol_P, s.individual_potential_temperatures.iloc[0]['Tp']),
        supra_sol_P
        ))
    np.savetxt(f, arr)
    


























# ---- check against plank & forsyth supplementary
# s = Suite("PF16_S7.csv", src_FeIII_totFe=0.19)
s = Suite("UT09DV04.csv", src_FeIII_totFe=0.17)
s.backtrack_compositions()
s.compute_pressure_temperature()
print()
print("Result from PF16 supplementary 7: P = 2.09 GPa, T = 1347 oC.")
print("Our result:                       P = %.2f GPa, T = %i oC." % 
    (s.PT['P'], s.PT['T']))



# ---- Brute force backtracking
__, comp = backtrack_sample_composition(s.data.iloc[0], verbose=True, return_all=True)

# fig, axs = plt.subplots(2,1,sharex=True)

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
#     axs[0].plot(ol_added, conc, label=p)
# axs[0].legend()
# axs[0].set_ylabel("Concentration / Initial concentration")
# 
# axs[1].plot(ol_added, Fo)
# axs[1].set_ylabel("Fo")
# axs[1].set_xlabel("Fraction olivine added")
# plt.show()

# ---- Save

with open(out_dir + "backtrack_concentrations_varKd.dat", "wb") as f:
    cols = ["blue", "0/125/0", "red", "darkgrey"]
    for i,phase in enumerate(['MgO', 'FeO', 'SiO2', 'Al2O3']):
        hdr = b"> -W0.7p,%s,-\n" % (cols[i].encode("ascii"))
        f.write(hdr)
        arr = np.column_stack(( ol_added*100., [c[phase]/comp[0][phase] for c in comp] ))
        np.savetxt(f, arr)

with open(out_dir + "labels_varKd.dat", "wb") as f:
    label = ["MgO", "FeO", "SiO@-2@-", "Other"]
    for i,phase in enumerate(['MgO', 'FeO', 'SiO2', 'Al2O3']):
        hdr = b"%f %f 7p,Helvetica,%s %s\n" % (
            ol_added[-1]*100.,
            comp[-1][phase]/comp[0][phase],
            cols[i].encode("ascii"),
            label[i].encode("ascii")
            )
        f.write(hdr)

with open(out_dir + "forsterite_varKd.dat", "wb") as f:
    arr = np.column_stack(( ol_added*100., np.array(Fo)*100. ))
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

with open(out_dir + "temperature_pressure_varKd.dat", "wb") as f:
    arr = np.column_stack(( s.PT['T'], s.PT['P'] ))
    np.savetxt(f, arr)

# with open(out_dir + "solidus.dat", "wb") as f:
#     arr = np.column_stack(( T_sol, P_sol ))
#     np.savetxt(f, arr)

with open(out_dir + "melt_path_varKd.dat", "wb") as f:
    hdr = b"> -L%i\n" % (s.individual_potential_temperatures.iloc[0]['Tp'])
    f.write(hdr)
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
    
with open(out_dir + "adiabat_varKd.dat", "wb") as f:
    hdr = b"> -L%i\n" % (s.individual_potential_temperatures.iloc[0]['Tp'])
    f.write(hdr)
    supra_sol_P = np.arange(mantle.solidusIntersection(s.individual_potential_temperatures.iloc[0]['Tp']), 0., -0.01)
    arr = np.column_stack((
        mantle.adiabat(supra_sol_P, s.individual_potential_temperatures.iloc[0]['Tp']),
        supra_sol_P
        ))
    np.savetxt(f, arr)
    
    
    
