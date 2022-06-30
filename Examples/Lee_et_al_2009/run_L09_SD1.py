from meltPT import *
import matplotlib.pyplot as plt
import sys

# ---- check against plank & forsyth supplementary
s = Suite("L09_data.csv", src_FeIII_totFe=0.05)
s.backtrack_compositions(target_Fo=0.9, dm=0.005, verbose=True)
s.compute_pressure_temperature(method="L09")
print()
print("Result from L09 supplementary dataset")
print("Result for sample s190 (line 18): P = 2.39 GPa, T = 1503 oC.")
print("Our result:                       P = %.2f GPa, T = %i oC." % 
    (s.PT['P'], s.PT['T']))

# ---- check reading in primaries
s2 = Suite("L09_data_primary.csv", read_as_primary=True)
s2.compute_pressure_temperature(method="L09")
print("Our result from primary comp:     P = %.2f GPa, T = %i oC." % 
    (s2.PT['P'], s2.PT['T']))