from meltPT import *
import matplotlib.pyplot as plt
import sys

# ---- Till et al 2012 
s = Suite("TGK12_data.csv", read_as_primary=True)
s.compute_pressure_temperature()
print()
print("Spinel stability field")
print("Result from TGK12 Table 6 Reykjanes: P = 1.16 GPa, T = 1316 oC.")
print("Our result:                   P = %.2f GPa, T = %i oC." % 
    (s.PT['P'][0], s.PT['T'][0]))

# ---- Till et al 2012 
s = Suite("TGK12_data.csv", read_as_primary=True)
s.compute_pressure_temperature()
print("------------")
print("Plagioclase stability field")
# print("Result from TGK12 Table 6 Reykjanes: P = 1.16 GPa, T = 1316 oC.")
# print("Our result (2):                   P = %.2f GPa, T = %i oC." % 
#     (s.PT['P'][0], s.PT['T'][0]))




