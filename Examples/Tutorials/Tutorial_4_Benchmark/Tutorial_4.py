#--------------------------
# Tutorial Four
#
# In this tutorial, we will use each thermobarometer
# in turn to estimate pressure and temperature for an
# example sample used in the original publication.
# This tutorial is a check to make sure they are coded 
# up correctly!
# Small deviations, e.g., 1-2 oC, may result from 
# rounding differences between implementations.
#--------------------------

# ---- Importing
from meltPT import *
import matplotlib.pyplot as plt
import sys

# ---- Using Thermometers
# B93, P07_2, and P07_4 can only be used as thermometers. 
# TGK12, G13, SD20, BK21 can also be used in this way.
# To calculate temperature, use the 'compute_temperature'
# function with the desired pressure in GPa, e.g., 
# for pressure = 0.8 GPa: 
# s.compute_temperature(P=0.8)

# Beattie, 1993 (B93)
# Comparing to result from Table 6 of Putirka (2008)
# Siqueros Sample
# P = 0.80 GPa, T = 1352 oC

s1 = Suite("Data/P07_Table_1.csv", read_as_primary=True)
s1.compute_temperature(method="B93", P=0.8)

print("B93 Documented Result:            P = 0.80 GPa, T = 1352 oC")
print("Our result at 0.80 GPa:           P = %.2f GPa, T = %i oC." % 
    (s1.PT['P'][0], s1.PT['T'][0]))
print()
print("-----------------------------------------------------------")

# Putirka et al, 2007 Equation 2 (P07_2)
# Comparing to result from Table 6 of Putirka (2008)
# MORB result

s2 = Suite("Data/P08_Table_6.csv", read_as_primary=True)
s2.compute_temperature(method="P07_2", P=1.56)
print("P07_2 Documented Result:          P = 1.56 GPa, T = 1374 oC")
print("Our result at 1.56 GPa:           P = %.2f GPa, T = %i oC." % 
    (s2.PT['P'][0], s2.PT['T'][0]))
print()
print("-----------------------------------------------------------")

# Putirka et al, 2007 Equation 4 (P07_4)
# Comparing to result from Table 6 of Putirka (2008)
# MORB result

s3 = Suite("Data/P07_Table_1.csv", read_as_primary=True)
s3.compute_temperature(method="P07_4", P=0.8)
print("P07_4 Documented Result:          P = 0.80 GPa, T = 1347 oC")
print("Our result at 1.56 GPa:           P = %.2f GPa, T = %i oC." % 
    (s3.PT['P'][0], s3.PT['T'][0]))
print()
print("-----------------------------------------------------------")

# ---- Using Thermobarometers
# Here, we estimate pressure and temperature in the same manner as the 
# previous tutorials.

# Lee et al, 2009
# Comparing to result for Sample s190 (The first sample in the
# spreadsheet attached to Lee et al., 2009).
# Comparing to answers from Lee is a chance to also check the 
# backtracking works correctly, unfortunately it does not. We
# think this is due to a mistake in the spreadsheet attached to
# Lee et al., (2009).

s4 = Suite("Data/L09_s190.csv", src_FeIII_totFe=0.05, Ce_to_H2O=200.)
s4.backtrack_compositions(target_Fo=0.9, dm=0.005, Kd=0.32166722315513757, verbose=False)
s4.compute_pressure_temperature(method="L09")
s4a = Suite("Data/L09_s190.csv", src_FeIII_totFe=0.05, Ce_to_H2O=200.)
s4a.backtrack_compositions(target_Fo=0.9, dm=0.005, verbose=False)
s4a.compute_pressure_temperature(method="L09")
print("L09 Documented Result             P = 2.39 GPa, T = 1503 oC.")
print("Our result (Kd=0.32167):          P = %.2f GPa, T = %i oC." %
    (s4.PT['P'], s4.PT['T']))
print("Our result (variable Kd):         P = %.2f GPa, T = %i oC." % 
    (s4a.PT['P'], s4a.PT['T']))

# Result using the primary composition outputted by the spreadsheet
s5 = Suite("Data/L09_s190_primary.csv", read_as_primary=True)
s5.compute_pressure_temperature(method="L09")
print("Our result from primary comp:     P = %.2f GPa, T = %i oC." % 
    (s5.PT['P'], s5.PT['T']))
print()

# Result backtracking using a fixed Kd value of 0.32.
s6 = Suite("Data/L09_s190.csv", src_FeIII_totFe=0.05, Ce_to_H2O=200.)
s6.backtrack_compositions(target_Fo=0.9, dm=0.005, Kd=0.32, verbose=False)
s6.compute_pressure_temperature(method="L09")
print("L09 Documented Result (Kd=0.32):  P = 2.35 GPa, T = 1498 oC.")
print("Our result:                       P = %.2f GPa, T = %i oC." % 
    (s6.PT['P'], s6.PT['T']))
print()
print("-----------------------------------------------------------")

# Till et al, 2012 
# Comparing to result for Reykjanes in Table 6 of Till et al., (2012).
s8 = Suite("Data/TGK12_T6.csv", read_as_primary=True)
s8.compute_pressure_temperature(method="TGK12_SPL")
print("TGK12 Documented Result:          P = 1.16 GPa, T = 1316 oC.")
print("Our result:                       P = %.2f GPa, T = %i oC." % 
    (s8.PT['P'][0], round(s8.PT['T'][0])))
print()
print("-----------------------------------------------------------")

# Grove et al, 2013
# Comparing to result in Table 5 for KR4003 experiment 40.07
# Note that this example onlytests the thermometer.
s9 = Suite("Data/G13_Table_5.csv", read_as_primary=True)
s9.compute_temperature(method="G13", P=4.)
print("G13 Documented Result:            P = 4.00 GPa, T = 1624 oC")
print("Our result at 4 GPa:              P = %.2f GPa, T = %i oC." % 
    (s9.PT['P'], s9.PT['T']))
print()
print("-----------------------------------------------------------")

# Herzberg & Asimow, 2015
# Comparing to result in Supplementary Iformation 1 of Herzberg &
# Asimow, (2015) for Mkea in MEGAPRIMELT3
# Again, this result only tests the thermometer.
s10 = Suite("Data/HA15_MEGAPRIMELT3.csv", read_as_primary=True)
s10.compute_temperature(method="HA15", P=0.)
print("HA15 Documented Result:           P = 0.00 GPa, T = 1431 oC")
print("Our result at 0 GPa:              P = %.2f GPa, T = %i oC." % 
    (s10.PT['P'][0], round(s10.PT['T'][0])))
print()
print("-----------------------------------------------------------")

# Plank and Forsyth, 2016
# Comparing to result in Supplementary S7 in Plank & Forsyth, (2016)
# Also provides an opportunity to test backtracking method

s11 = Suite("Data/PF16_S7.csv", src_FeIII_totFe=0.19, Ce_to_H2O=200.)
s11.backtrack_compositions(Kd=0.3, verbose=False)
s11.compute_pressure_temperature()
print("PF16 Documented Result:           P = 2.09 GPa, T = 1347 oC.")
print("Our result:                       P = %.2f GPa, T = %i oC." % 
    (s11.PT['P'], s11.PT['T']))
print()
print("-----------------------------------------------------------")

# Sun & Dasgupta, 2020 (SD20)
# Comparing to result in Supplementary 1 of Sun & Dasgupta, (2020)
# for Longhi (1995) sample.
s13 = Suite("Data/SD20_S1.csv", read_as_primary=True)
s13.compute_pressure_temperature(method="SD20")
print("SD20 Documented Result:           P = 2.90 GPa, T = 1529 oC.")
print("Our result:                       P = %.2f GPa, T = %i oC." % 
    (s13.PT['P'], s13.PT['T']))
print()
print("-----------------------------------------------------------")

# Brown Krein et al, 2021 (BK21)
# Comparing to result in Supplementary Table S1 of Sun & Dasgupta (2020)
# for KLB-1 experiment 18 (Hirose & Kushiro, 1993). This sample is in 
# the spinel stability field.
s14 = Suite("Data/K21_KLB2.csv", read_as_primary=True)
s14.compute_pressure_temperature(method="BK21")
print("BK21 Documented Result:           P = 1.51 GPa, T = 1326 oC")
print("Our result :                      P = %.2f GPa, T = %i oC." % 
    (s14.PT['P'], s14.PT['T']))
print()

# Comparing to result in Supplementary Table S1 of Sun & Dasgupta (2020)
# for KR4003 experiment 40.06 (Walter, 1998). This sample is in the 
# garnet stability field.
s15 = Suite("Data/K21_KR4003.csv", read_as_primary=True)
s15.compute_pressure_temperature(method="BK21")
print("BK21 Documented Result:           P = 3.56 GPa, T = 1584 oC")
print("Our result :                      P = %.2f GPa, T = %i oC." % 
    (s15.PT['P'], s15.PT['T']))
print()
print("-----------------------------------------------------------")
