from meltPT import *
import matplotlib.pyplot as plt
import sys

# ---- Set up mantle
lz = m.lithologies.katz.lherzolite()
mantle = m.mantle([lz], [1], ['Lz'])
P_sol = np.arange(0., 3., 0.1)
T_sol = [lz.TSolidus(P) for P in P_sol]

# ---- Plank and Forsyth, 2016
s = Suite("KLB1_43.csv", read_as_primary=True)
s.compute_pressure_temperature()
print()
print("Plank & Forsyth (2016):          P = %.2f GPa, T = %i oC." % 
    (s.PT['P'], s.PT['T']))
s.find_individual_potential_temperatures(mantle)
print("Best-fitting melting model:      Tp = %i oC, F = %.2f %%." % 
    (s.individual_potential_temperatures['Tp'], 
    s.individual_potential_temperatures['F']*100.))
print()

# ---- Beattie, 1993
s2 = Suite("K21_KR4003_3012.csv", read_as_primary=True)
s2.compute_temperature(method="B93", P=3.)
print("Beattie (1993):                  P = %.2f GPa, T = %i oC." % 
    (s2.PT['P'], s2.PT['T']))
s2.find_individual_potential_temperatures(mantle)
print("Best-fitting melting model:      Tp = %i oC, F = %.2f %%." % 
    (s2.individual_potential_temperatures['Tp'], 
    s2.individual_potential_temperatures['F']*100.))
print()

# ---- Lee et al 2009
s3 = Suite("K21_KR4003_3012.csv", read_as_primary=True, src_FeIII_totFe=0.15)
s3.compute_pressure_temperature(method="L09")
print("Lee et al (2009):                P = %.2f GPa, T = %i oC." % 
    (s3.PT['P'], s3.PT['T']))
s3.find_individual_potential_temperatures(mantle)
print("Best-fitting melting model:      Tp = %i oC, F = %.2f %%." % 
    (s3.individual_potential_temperatures['Tp'], 
    s3.individual_potential_temperatures['F']*100.))
print()

# ---- Till et al 2012 
s4 = Suite("K21_KR4003_3012.csv", read_as_primary=True)
s4.compute_pressure_temperature(method="TGK12_SPL")
print("Till et al (2012) SPL:           P = %.2f GPa, T = %i oC." % 
    (s4.PT['P'], round(s4.PT['T'])))
s4.find_individual_potential_temperatures(mantle)
print("Best-fitting melting model:      Tp = %i oC, F = %.2f %%." % 
    (s4.individual_potential_temperatures['Tp'], 
    s4.individual_potential_temperatures['F']*100.))
print()

s5 = Suite("K21_KR4003_3012.csv", read_as_primary=True)
s5.compute_pressure_temperature(method="TGK12_PLG")
print("Till et al (2012) PLG:           P = %.2f GPa, T = %i oC." % 
    (s5.PT['P'], round(s5.PT['T'])))
s5.find_individual_potential_temperatures(mantle)
print("Best-fitting melting model:      Tp = %i oC, F = %.2f %%." % 
    (s5.individual_potential_temperatures['Tp'], 
    s5.individual_potential_temperatures['F']*100.))
print()

# ---- Grove et al 2013
s6 = Suite("K21_KR4003_3012.csv", read_as_primary=True)
s6.compute_pressure_temperature(method="G13")
print("Grove et al (2013):              P = %.2f GPa, T = %i oC." % 
    (s6.PT['P'], s6.PT['T']))
s6.find_individual_potential_temperatures(mantle)
print("Best-fitting melting model:      Tp = %i oC, F = %.2f %%." % 
    (s6.individual_potential_temperatures['Tp'], 
    s6.individual_potential_temperatures['F']*100.))
print()

# ---- Krein et al 2021
s7 = Suite("K21_KR4003_3012.csv", read_as_primary=True)
s7.compute_pressure_temperature(method="BK21")
print("Krein et al (2021):              P = %.2f GPa, T = %i oC." % 
    (s7.PT['P'], s7.PT['T']))
s7.find_individual_potential_temperatures(mantle)
print("Best-fitting melting model:      Tp = %i oC, F = %.2f %%." % 
    (s7.individual_potential_temperatures['Tp'], 
    s7.individual_potential_temperatures['F']*100.))
print()

# ---- Herzberg & Asimow 2015
s8 = Suite("K21_KR4003_3012.csv", read_as_primary=True)
s8.compute_temperature(method="HA15", P=3.)
print("Herzberg & Asimow (2015):        P = %.2f GPa, T = %i oC." % 
    (s8.PT['P'], round(s8.PT['T'])))
s8.find_individual_potential_temperatures(mantle)
print("Best-fitting melting model:      Tp = %i oC, F = %.2f %%." % 
    (s8.individual_potential_temperatures['Tp'], 
    s8.individual_potential_temperatures['F']*100.))
print()

# ---- Putirka et al 2007
s9 = Suite("K21_KR4003_3012.csv", read_as_primary=True)
s9.compute_pressure_temperature(method="P07_P08")
print("Putirka et al (2008):            P = %.2f GPa, T = %i oC." % 
    (s9.PT['P'], s9.PT['T']))
s9.find_individual_potential_temperatures(mantle)
print("Best-fitting melting model:      Tp = %i oC, F = %.2f %%." % 
    (s9.individual_potential_temperatures['Tp'], 
    s9.individual_potential_temperatures['F']*100.))
print()

# ---- Sun & Dasgupta, 2020
s9 = Suite("K21_KR4003_3012.csv", read_as_primary=True)
s9.compute_pressure_temperature(method="SD20")
print("Sun & Dasgupta (2020):            P = %.2f GPa, T = %i oC." % 
    (s9.PT['P'], s9.PT['T']))
s9.find_individual_potential_temperatures(mantle)
print("Best-fitting melting model:      Tp = %i oC, F = %.2f %%." % 
    (s9.individual_potential_temperatures['Tp'], 
    s9.individual_potential_temperatures['F']*100.))
print()

# # ---- Putirka et al 2007
# s10 = Suite("K21_KLB2.csv", read_as_primary=True)
# s10.compute_temperature(method="P07_4", P=1.5)
# print("Putirka et al (2007) Equation 4: P = %.2f GPa, T = %i oC." % 
#     (s10.PT['P'], s10.PT['T']))
# s10.find_individual_potential_temperatures(mantle)
# print("Best-fitting melting model:      Tp = %i oC, F = %.2f %%." % 
#     (s10.individual_potential_temperatures['Tp'], 
#     s10.individual_potential_temperatures['F']*100.))
# print()

# s11 = Suite("K21_KLB2.csv", read_as_primary=True)
# s11.compute_temperature(method="P07_2", P=1.5)
# print("Putirka et al (2007) Equation 2: P = %.2f GPa, T = %i oC." % 
#     (s11.PT['P'], s11.PT['T']))
# s11.find_individual_potential_temperatures(mantle)
# print("Best-fitting melting model:      Tp = %i oC, F = %.2f %%." % 
#     (s11.individual_potential_temperatures['Tp'], 
#     s11.individual_potential_temperatures['F']*100.))
# print()

