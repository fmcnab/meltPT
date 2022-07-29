from meltPT import *
import matplotlib.pyplot as plt
import sys
from scipy.interpolate import interp1d, splev, splrep

# ---- check against plank & forsyth supplementary
s = Suite("PF16_S7.csv", src_FeIII_totFe=0.19)
s.backtrack_compositions(Kd=0.3)
s.compute_pressure_temperature()
print()
print("Result from PF16 supplementary 7: P = 2.09 GPa, T = 1347 oC.")
print("Our result:                       P = %.2f GPa, T = %i oC." % 
    (s.PT['P'], s.PT['T']))


# ---- Set up mantle
lz = m.lithologies.katz.lherzolite()
mantle = m.mantle([lz], [1], ['Lz'])


# ---- Geotherm

class Lithospheric_Asthenospheric_Geotherm:
    
    def __init__(self, P_lab, T_p, mantle, P_max=10.):
        self.P_lab = P_lab
        self.T_lab = mantle.adiabat(P_lab, Tp)
        self.T_p = T_p
        self.mantle = mantle
        self.P_max = P_max
        self.find_conductive_solidus_intercept()
        self.construct_geotherm()
        
    def compute_conductive_T(self, P):
        return (self.T_lab/self.P_lab) * P
        
    def find_conductive_solidus_intercept(self):
        
        def f(P):
            cond = self.compute_conductive_T(P)
            sol = self.mantle.lithologies[0].TSolidus(P)
            return np.sqrt( (cond - sol)**2. )
        
        fit = minimize_scalar(f)
        self.conductive_solidus_intercept = {
            'P': fit.x,
            'T': self.compute_conductive_T(fit.x)}
        
    def construct_geotherm(self):
        
        if self.T_lab < self.mantle.lithologies[0].TSolidus(self.P_lab):
            P_cond = np.array([0., self.P_lab])
            T_cond = self.compute_conductive_T(P_cond)
            P_asth = np.linspace(self.P_lab, self.P_max, 100.)
            T_asth = self.mantle.adiabat(P_asth, self.T_p)
            P_interp = np.array([])
            T_interp = np.array([])
        else:
            P_cond = np.array([0., self.conductive_solidus_intercept['P']])[:-1]
            T_cond = self.compute_conductive_T(P_cond)
            P_asth = np.linspace(self.mantle.solidusIntersection(self.T_p)[0], self.P_max, 100.)[1:]
            T_asth = self.mantle.adiabat(P_asth, self.T_p)
            n=100
            P_cond_interp = np.linspace(self.conductive_solidus_intercept['P'], self.P_lab, n)
            T_cond_interp = self.compute_conductive_T(P_cond_interp)
            P_asth_interp = np.linspace(self.P_lab, self.mantle.solidusIntersection(self.T_p)[0], n)
            T_asth_interp = self.mantle.adiabat(P_asth_interp, self.T_p)
            P_interp = np.zeros(n)
            T_interp = np.zeros(n)
            for i in range(n):
                P_interp[i] = (1. - i/n) * P_cond_interp[i] + (i/n) * P_asth_interp[i]
                T_interp[i] = (1. - i/n) * T_cond_interp[i] + (i/n) * T_asth_interp[i]

        P_stack = np.hstack((P_cond, P_interp, P_asth))
        T_stack = np.hstack((T_cond, T_interp, T_asth))
        
        self.T = interp1d(P_stack, T_stack)

P = np.arange(0., 5., 0.01)
Pb = 1.5
Tp = 1350.
g = Lithospheric_Asthenospheric_Geotherm(Pb, Tp, mantle)

plt.plot(g.compute_conductive_T(np.array([P.min(),Pb])), [P.min(),Pb], ":")
plt.plot(mantle.adiabat(np.hstack(([Pb, P.max()])), Tp), [Pb, P.max()], ":")
plt.plot(g.T(P), P)
plt.plot(mantle.lithologies[0].TSolidus(P), P, "--")
plt.gca().invert_yaxis()
plt.show()




