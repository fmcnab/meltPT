from meltPT import *
import matplotlib.pyplot as plt
import sys

# ---- Plank and Forsyth, 2016
s = Suite("K21_KLB2.csv", read_as_primary=True)
s.compute_pressure_temperature()
print()
print("Plank & Forsyth (2016):          P = %.2f GPa, T = %i oC." % 
    (s.PT['P'], s.PT['T']))

# ---- Beattie, 1993
s2 = Suite("K21_KLB2.csv", read_as_primary=True)
s2.compute_temperature(method="B93", P=1.5)
print("Beattie (1993):                  P = %.2f GPa, T = %i oC." % 
    (s2.PT['P'], s2.PT['T']))

# ---- Lee et al 2009
s3 = Suite("K21_KLB2.csv", read_as_primary=True)
s3.compute_pressure_temperature(method="L09")
print("Lee et al (2009):                P = %.2f GPa, T = %i oC." % 
    (s3.PT['P'], s3.PT['T']))

# ---- Till et al 2012 
s7 = Suite("K21_KLB2.csv", read_as_primary=True)
s7.compute_pressure_temperature(method="TGK12_SPL")
print("Till et al (2012) SPL:           P = %.2f GPa, T = %i oC." % 
    (s7.PT['P'], round(s7.PT['T'])))

s7 = Suite("K21_KLB2.csv", read_as_primary=True)
s7.compute_pressure_temperature(method="TGK12_PLG")
print("Till et al (2012) PLG:           P = %.2f GPa, T = %i oC." % 
    (s7.PT['P'], round(s7.PT['T'])))

# ---- Grove et al 2013
s8 = Suite("K21_KLB2.csv", read_as_primary=True)
s8.compute_pressure_temperature(method="G13")
print("Grove et al (2013):              P = %.2f GPa, T = %i oC." % 
    (s8.PT['P'], s8.PT['T']))

# ---- Krein et al 2021
s9 = Suite("K21_KLB2.csv", read_as_primary=True)
s9.compute_pressure_temperature(method="BK21")
print("Krein et al (2021):              P = %.2f GPa, T = %i oC." % 
    (s9.PT['P'], s9.PT['T']))

# ---- Herzberg & Asimow 2015
s12 = Suite("K21_KLB2.csv", read_as_primary=True)
s12.compute_temperature(method="HA15", P=1.5)
print("Herzberg & Asimow (2015):        P = %.2f GPa, T = %i oC." % 
    (s12.PT['P'], round(s12.PT['T'])))

# ---- Putirka et al 2007
s13 = Suite("K21_KLB2.csv", read_as_primary=True)
s13.compute_pressure_temperature(method="P07_P08")
print("Putirka et al (2008):            P = %.2f GPa, T = %i oC." % 
    (s13.PT['P'], s13.PT['T']))

# ---- Putirka et al 2007
s14 = Suite("K21_KLB2.csv", read_as_primary=True)
s14.compute_temperature(method="P07_4", P=1.5)
print("Putirka et al (2007) Equation 4: P = %.2f GPa, T = %i oC." % 
    (s14.PT['P'], s14.PT['T']))

s15 = Suite("K21_KLB2.csv", read_as_primary=True)
s15.compute_temperature(method="P07_2", P=1.5)
print("Putirka et al (2007) Equation 2: P = %.2f GPa, T = %i oC." % 
    (s15.PT['P'], s15.PT['T']))