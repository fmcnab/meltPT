from meltPT import *
import matplotlib.pyplot as plt
import sys

# ---- Plank and Forsyth, 2016
s = Suite("PF16_S7.csv", src_FeIII_totFe=0.19)
s.backtrack_compositions(Kd=0.3, verbose=False)
s.compute_pressure_temperature()
print("Plank & Forsyth (2016)")
print()
print("Result from Supplementary 7:      P = 2.09 GPa, T = 1347 oC.")
print("Our result:                       P = %.2f GPa, T = %i oC." % 
    (s.PT['P'], s.PT['T']))
# ---- check reading in primaries
s2 = Suite("PF16_S7_Primary.csv", read_as_primary=True)
s2.compute_pressure_temperature()
print("Our result from primary comp:     P = %.2f GPa, T = %i oC." % 
    (s2.PT['P'], s2.PT['T']))
print("-----------------------------------------------------------")

# ---- Lee et al 2009
s3 = Suite("L09_s190.csv", src_FeIII_totFe=0.05)
s3.backtrack_compositions(target_Fo=0.9, dm=0.005, verbose=False)
s3.compute_pressure_temperature(method="L09")
print("Lee et al (2009)")
print()
print("Result from supplementary dataset")
print("for Sample s190:                  P = 2.39 GPa, T = 1503 oC.")
print("Our result:                       P = %.2f GPa, T = %i oC." % 
    (s3.PT['P'], s3.PT['T']))

# ---- Lee et al 2009
s3 = Suite("L09_s190.csv", src_FeIII_totFe=0.05)
s3.backtrack_compositions(target_Fo=0.9, dm=0.005, Kd=0.32, verbose=False)
s3.compute_pressure_temperature(method="L09")
print("Dataset result with fixed Kd=0.32:P = 2.35 GPa, T = 1498 oC.")
print("Our result:                       P = %.2f GPa, T = %i oC." % 
    (s3.PT['P'], s3.PT['T']))

# ---- check reading in primaries
s4 = Suite("L09_s190_primary.csv", read_as_primary=True)
s4.compute_pressure_temperature(method="L09")
print("Our result from primary comp:     P = %.2f GPa, T = %i oC." % 
    (s4.PT['P'], s4.PT['T']))
print()
print("-----------------------------------------------------------")

# ---- Till et al 2012 
s5 = Suite("TGK12_T6.csv", read_as_primary=True)
s5.compute_pressure_temperature(method="TGK12_SPL")
print("Till et al (2012) - Spinel is stable")
print()
print("Result from Table 6 for Reykjanes:P = 1.16 GPa, T = 1316 oC.")
print("Our result:                       P = %.2f GPa, T = %i oC." % 
    (s5.PT['P'][0], round(s5.PT['T'][0])))
print()
print("-----------------------------------------------------------")

# ---- Grove et al 2013
s6 = Suite("K21_D-20-15.csv", read_as_primary=True)
s6.compute_pressure_temperature(method="K21_PLG")
print("Grove et al 2013")
print()
print("Result in Figure 3:               P = 0.95 GPa, T = 1298 oC")
print("Our result:                       P = %.2f GPa, T = %i oC." % 
    (s6.PT['P'], s6.PT['T']))
print()
print("-----------------------------------------------------------")

# ---- Krein et al 2021
s7 = Suite("K21_D-20-15.csv", read_as_primary=True)
s7.compute_pressure_temperature(method="K21_PLG")
print("Krein et al (2021)") 
print()
print("- Plagioclase is stable")
print()
print("Result in Figure 3 for D-20-15:   P = 0.95 GPa, T = 1298 oC")
print("Our result :                      P = %.2f GPa, T = %i oC." % 
    (s7.PT['P'], s7.PT['T']))
print()

s8 = Suite("K21_KLB2.csv", read_as_primary=True)
s8.compute_pressure_temperature(method="K21_SPL")
print("- Spinel is stable")
print()
print("Result in Table S1 for KLB-1")
print("experiment 18")
print("(Hirose & Kushiro, 1993):         P = 1.51 GPa, T = 1326 oC")
print("Our result :                      P = %.2f GPa, T = %i oC." % 
    (s8.PT['P'], s8.PT['T']))
print()

s9 = Suite("K21_KR4003.csv", read_as_primary=True)
s9.compute_pressure_temperature(method="K21_GNT")
print("- Garnet is stable")
print()
print("Result in Table S1 for KR4003")
print("experiment 40.06 (Walter 1998):   P = 3.56 GPa, T = 1584 oC")
print("Our result :                      P = %.2f GPa, T = %i oC." % 
    (s9.PT['P'], s9.PT['T']))
print()
print("-----------------------------------------------------------")