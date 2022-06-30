from meltPT import *
import matplotlib.pyplot as plt
import sys

# ---- check reading in primaries
s2 = Suite("Example_1.csv", read_as_primary=True)
s2.compute_pressure_temperature(method="K21_GNT")
print("Result in Figure 3:               P = 0.95 GPA, T = 1298 oC")
print("Our result (2):                   P = %.2f GPa, T = %i oC." % 
    (s2.PT['P'], s2.PT['T']))