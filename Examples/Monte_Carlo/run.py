from meltPT import *
import pyMelt as m
import matplotlib.pyplot as plt
import random

class Value:
    def __init__(self, val, plus, minus):
        self.val = val
        self.plus = plus
        self.minus = minus
    
    def draw_normal(self):
        return random.gauss(self.val, (self.plus+self.minus)/2.)
        
    def range(self, n=10):
        return np.linspace(self.val-(2.*self.minus), self.val+(2.*self.plus), n)
        
    def gaussian(self):
        sd = (self.plus+self.minus)/2.
        return (
            (1./(sd*np.sqrt(2.*np.pi))) * 
            np.exp(-0.5*(((self.range(n=100)-self.val)/sd)**2.))
            )

def gaussian(x, val, sd):
    return (1./(sd*np.sqrt(2.*np.pi)))*np.exp(-0.5*(((x-val)/sd)**2.))

# Path to input file
infile = "../../meltPT/Examples/Tutorials/Data/PF16_UT09DV04.csv"

# Properties
src_FeIII_totFe = Value(0.17, 0.01, 0.01)
src_Fo = Value(0.9, 0.01, 0.01)
Kd = Value(0.3, 0.01, 0.01)
H2O = Value(1.376785, 0.1, 0.1)


# Mantle
lz = m.lithologies.katz.lherzolite()
mantle = m.mantle([lz], [1], ['Lz'])
P = np.arange(1., 3.5, 0.01)

# Basic run
s = Suite(infile, src_FeIII_totFe=src_FeIII_totFe.val)
s.backtrack_compositions(Kd=0.3, target_Fo=src_Fo.val)
s.compute_pressure_temperature(method="PF16")
s.find_individual_potential_temperatures(mantle)

# Isolated: Fo
suites_Fo = []
for Fo in src_Fo.range():
    si = Suite(infile, src_FeIII_totFe=src_FeIII_totFe.val)
    si.backtrack_compositions(Kd=Kd.val, target_Fo=Fo)
    si.compute_pressure_temperature(method="PF16")
    si.find_individual_potential_temperatures(mantle)
    suites_Fo.append(si)
    
# Isolated: Fe
suites_Fe = []
for Fe in src_FeIII_totFe.range():
    si = Suite(infile, src_FeIII_totFe=Fe)
    si.backtrack_compositions(Kd=Kd.val, target_Fo=src_Fo.val)
    si.compute_pressure_temperature(method="PF16")
    si.find_individual_potential_temperatures(mantle)
    suites_Fe.append(si)
    
# Isolated: Kd
suites_Kd = []
for kd in Kd.range():
    si = Suite(infile, src_FeIII_totFe=src_FeIII_totFe.val)
    si.backtrack_compositions(Kd=kd, target_Fo=src_Fo.val)
    si.compute_pressure_temperature(method="PF16")
    si.find_individual_potential_temperatures(mantle)
    suites_Kd.append(si)
    
# Isolated: Water
suites_H2O = []
for h2o in H2O.range():
    si = Suite(infile, src_FeIII_totFe=src_FeIII_totFe.val)
    si.data.at[0, 'H2O'] = h2o
    si.backtrack_compositions(Kd=Kd.val, target_Fo=src_Fo.val)
    si.compute_pressure_temperature(method="PF16")
    si.find_individual_potential_temperatures(mantle)
    suites_H2O.append(si)


# Plot
fig, axs = plt.subplots(2, 4, sharey="row")
axs[0,0].plot(lz.TSolidus(P), P, c="0.75")
axs[0,0].scatter(
    [s.PT['T'].iloc[0] for s in suites_Fo],
    [s.PT['P'].iloc[0] for s in suites_Fo],
    c=src_Fo.range(),
    s=10)
axs[0,0].errorbar(s.PT['T'], s.PT['P'], xerr=s.PT['T_err'], yerr=s.PT['P_err'], marker="o")
# axs[0,0].invert_yaxis()
col = axs[1,0].scatter(
    src_Fo.range(), 
    [s.individual_potential_temperatures['Tp'].iloc[0] for s in suites_Fo],
    c=src_Fo.range(),
    s=10)
axs[1,0].scatter(
    src_Fo.val, 
    s.individual_potential_temperatures['Tp'].iloc[0])
axs[0,0].set_xlabel(r"$T$ [$^\circ$C]")
axs[0,0].set_ylabel(r"$P$ [GPa]")
axs[1,0].set_ylabel(r"$T_p$ [$^\circ$C]")
axs[1,0].set_xlabel("Fo#")
# plt.colorbar(col)

axs[0,1].plot(lz.TSolidus(P), P, c="0.75")
axs[0,1].scatter(
    [s.PT['T'].iloc[0] for s in suites_Fe],
    [s.PT['P'].iloc[0] for s in suites_Fe],
    c=src_FeIII_totFe.range(),
    s=10)
axs[0,1].errorbar(s.PT['T'], s.PT['P'], xerr=s.PT['T_err'], yerr=s.PT['P_err'], marker="o")
col = axs[1,1].scatter(
    src_FeIII_totFe.range(), 
    [s.individual_potential_temperatures['Tp'].iloc[0] for s in suites_Fe],
    c=src_FeIII_totFe.range(),
    s=10)
axs[1,1].scatter(
    src_FeIII_totFe.val, 
    s.individual_potential_temperatures['Tp'].iloc[0])
# plt.colorbar(col)
axs[0,1].set_xlabel(r"$T$ [$^\circ$C]")
axs[1,1].set_xlabel(r"Fe$^{3+}$/$\Sigma$Fe")


axs[0,2].plot(lz.TSolidus(P), P, c="0.75")
axs[0,2].scatter(
    [s.PT['T'].iloc[0] for s in suites_Kd],
    [s.PT['P'].iloc[0] for s in suites_Kd],
    c=Kd.range(),
    s=10)
axs[0,2].errorbar(s.PT['T'], s.PT['P'], xerr=s.PT['T_err'], yerr=s.PT['P_err'], marker="o")
col = axs[1,2].scatter(
    Kd.range(), 
    [s.individual_potential_temperatures['Tp'].iloc[0] for s in suites_Kd],
    c=Kd.range(),
    s=10)
axs[1,2].scatter(
    Kd.val, 
    s.individual_potential_temperatures['Tp'].iloc[0])
# plt.colorbar(col)
axs[0,2].set_xlabel(r"$T$ [$^\circ$C]")
axs[1,2].set_xlabel(r"$K_d$")


axs[0,3].plot(lz.TSolidus(P), P, c="0.75")
axs[0,3].scatter(
    [s.PT['T'].iloc[0] for s in suites_H2O],
    [s.PT['P'].iloc[0] for s in suites_H2O],
    c=H2O.range(),
    s=10)
axs[0,3].errorbar(s.PT['T'], s.PT['P'], xerr=s.PT['T_err'], yerr=s.PT['P_err'], marker="o")
axs[0,3].invert_yaxis()
col = axs[1,3].scatter(
    H2O.range(), 
    [s.individual_potential_temperatures['Tp'].iloc[0] for s in suites_H2O],
    c=H2O.range(),
    s=10)
axs[1,3].scatter(
    H2O.val, 
    s.individual_potential_temperatures['Tp'].iloc[0])
# plt.colorbar(col)
axs[0,3].set_xlabel(r"$T$ [$^\circ$C]")
axs[1,3].set_xlabel(r"H$_2$O")

for ax in axs:
    for a in ax:
        a.set_box_aspect(1)

plt.show()

# Monte Carlo
ss = []
Fos = []
Kds = []
for i in range(100):
    print(i)
    Fo = src_Fo.draw_normal()
    kd = Kd.draw_normal()
    # Fo = 0.9
    si = Suite(infile, src_FeIII_totFe=src_FeIII_totFe.draw_normal())
    si.backtrack_compositions(Kd=kd, target_Fo=src_Fo.draw_normal())
    si.compute_pressure_temperature(method="PF16")
    si.data.at[0, 'H2O'] = H2O.draw_normal()
    si.PT.at[0,'P'] = Value(s.PT.at[0,'P'], s.PT.at[0,'P_err'], s.PT.at[0,'P_err']).draw_normal()
    si.PT.at[0,'T'] = Value(s.PT.at[0,'T'], s.PT.at[0,'T_err'], s.PT.at[0,'T_err']).draw_normal()
    si.find_individual_potential_temperatures(mantle)
    ss.append(si)
    Fos.append(Fo)
    Kds.append(kd)

 
# Plot

fig, axs = plt.subplots(2, 4, sharey=True)
gs = axs[0,2].get_gridspec()
for ax in axs[0:,2:]:
    for a in ax:
        a.remove()
PT = fig.add_subplot(gs[0:, -2])
Tp = fig.add_subplot(gs[0:, -1])

h = axs[0,0].hist(Fos)
axs[0,0].plot(src_Fo.range(n=100), src_Fo.gaussian() / max(src_Fo.gaussian()) * max(h[0]))
axs[0,0].set_xlabel("Source Fo#")
axs[0,0].set_ylabel("Density")
axs[0,0].set_box_aspect(1)

h = axs[0,1].hist([s.data['src_FeIII_totFe'].iloc[0] for s in ss])
axs[0,1].plot(src_FeIII_totFe.range(n=100), src_FeIII_totFe.gaussian() / max(src_FeIII_totFe.gaussian()) * max(h[0]))
axs[0,1].set_xlabel(r"Fe$^{3+}$/$\Sigma$Fe")
axs[0,1].set_box_aspect(1)

h = axs[1,0].hist(Kds)
axs[1,0].plot(Kd.range(n=100), Kd.gaussian() / max(Kd.gaussian()) * max(h[0]))
axs[1,0].set_xlabel("$K_d$")
axs[1,0].set_ylabel("Density")
axs[1,0].set_box_aspect(1)

h = axs[1,1].hist([s.data['H2O'].iloc[0] for s in ss])
axs[1,1].plot(H2O.range(n=100), H2O.gaussian() / max(H2O.gaussian()) * max(h[0]))
axs[1,1].set_xlabel(r"H$_2$O [wt%]")
axs[1,1].set_box_aspect(1)

P = np.arange(1., 3.5, 0.01)
PT.plot(lz.TSolidus(P), P, c="0.75")
PT.scatter(
    [s.PT['T'].iloc[0] for s in ss], 
    [s.PT['P'].iloc[0] for s in ss], 
    s=0.2)
PT.errorbar(s.PT['T'], s.PT['P'], xerr=s.PT['T_err'], yerr=s.PT['P_err'], marker="o")
PT.invert_yaxis()
PT.set_xlabel(r"$T$ [$^\circ$C]")
PT.set_ylabel(r"$P$ [GPa]")

Tps = [Tp for Tp in [s.individual_potential_temperatures['Tp'].iloc[0] for s in ss] if np.isfinite(Tp)]
Tps, counts = np.unique(Tps, return_counts=True)
cumsum = np.cumsum(counts)/len(counts)
minTp = np.interp(1./3., cumsum, Tps)
maxTp = np.interp(2./3., cumsum, Tps)
Tp.hist(Tps)
lab = r"_{33-66}T_p, set"
f = Tp.fill(
    [minTp, minTp, maxTp, maxTp],
    [0., 30., 30., 0.],
    alpha=0.5,
    label=lab
)
l2 = Tp.plot([s.individual_potential_temperatures['Tp'].iloc[0]]*2, [0, 30], label=lab)
l = Tp.plot([np.median(Tps)]*2, [0, 30], "--", label=lab)
Tp.set_ylim(0,30)
Tp.set_xlabel(r"$T_p$ [$^\circ$C]")
Tp.set_ylabel(r"Count")
Tp.legend([(f[0], l[0]), (l2[0],), ], [r"$T_p$, set", r"$T_p$, individual", ])

# plt.tight_layout()
plt.show()


# fig, axs = plt.subplots(1, 4)
# 
# axs[0].hist([s.data['src_FeIII_totFe'].iloc[0] for s in ss])
# axs[0].plot(
#     np.arange(0.14, 0.2, 0.001), 
#     gaussian(
#         np.arange(0.14, 0.2, 0.001), 
#         src_FeIII_totFe.val,
#         (src_FeIII_totFe.plus+src_FeIII_totFe.minus)/2.))
# 
# axs[1].hist(Fos)
# axs[1].plot(
#     np.arange(0.87, 0.93, 0.001), 
#     gaussian(
#         np.arange(0.87, 0.93, 0.001), 
#         src_Fo.val,
#         (src_Fo.plus+src_Fo.minus)/2.))
# 
# 
# lz = m.lithologies.katz.lherzolite()
# P = np.arange(1., 3.5, 0.01)
# axs[2].plot(lz.TSolidus(P), P, c="0.75")
# axs[2].scatter(
#     [s.PT['T'].iloc[0] for s in ss], 
#     [s.PT['P'].iloc[0] for s in ss], 
#     s=0.2)
# axs[2].errorbar(s.PT['T'], s.PT['P'], xerr=s.PT['T_err'], yerr=s.PT['P_err'], marker="o")
# axs[2].invert_yaxis()
# 
# axs[3].hist([s.individual_potential_temperatures['Tp'].iloc[0] for s in ss])
# plt.show()