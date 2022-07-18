from meltPT import *
import matplotlib.pyplot as plt
import sys

# ---- Beattie, 1993
s1 = Suite("P07_Table_1.csv", read_as_primary=True)
s1.compute_temperature(method="B93", P=0.8)
print()
print("Beattie (1993)")
print()
print("Result in Table 6 of Putirka (2008)")
print("for Siqueiros:                    P = 0.80 GPa, T = 1352 oC")
print("Our result at 0.80 GPa:           P = %.2f GPa, T = %i oC." % 
    (s1.PT['P'][0], s1.PT['T'][0]))
print()
print("-----------------------------------------------------------")


# ---- Putirka et al, 2007
s2 = Suite("P08_Table_6.csv", read_as_primary=True)
s2.compute_temperature(method="P07_2", P=1.56)
print()
print("Putirka et al (2007) Equation 2")
print()
print("Result in Table 6 of Putirka (2008)")
print("for MORB:                         P = 1.56 GPa, T = 1374 oC")
print("Our result at 1.56 GPa:           P = %.2f GPa, T = %i oC." % 
    (s2.PT['P'][0], s2.PT['T'][0]))
print()
print("-----------------------------------------------------------")

s3 = Suite("P07_Table_1.csv", read_as_primary=True)
s3.compute_temperature(method="P07_4", P=0.8)
print("Putirka et al (2007) Equation 4")
print()
print("Result in Table 6 of Putirka (2008)")
print("for MORB:                         P = 0.80 GPa, T = 1347 oC")
print("Our result at 1.56 GPa:           P = %.2f GPa, T = %i oC." % 
    (s3.PT['P'][0], s3.PT['T'][0]))
print()
print("-----------------------------------------------------------")

# ---- Lee et al, 2009
s4 = Suite("L09_s190.csv", src_FeIII_totFe=0.05, Ce_to_H2O=200.)
s4.backtrack_compositions(target_Fo=0.9, dm=0.005, verbose=False)
s4.compute_pressure_temperature(method="L09")
print("Lee et al (2009)")
print()
print("Result from supplementary dataset")
print("for Sample s190:                  P = 2.39 GPa, T = 1503 oC.")
print("Our result:                       P = %.2f GPa, T = %i oC." % 
    (s4.PT['P'], s4.PT['T']))
# ---- check reading in primaries
s5 = Suite("L09_s190_primary.csv", read_as_primary=True)
s5.compute_pressure_temperature(method="L09")
print("Our result from primary comp:     P = %.2f GPa, T = %i oC." % 
    (s5.PT['P'], s5.PT['T']))
print()

# ---- Lee et al, 2009
s6 = Suite("L09_s190.csv", src_FeIII_totFe=0.05, Ce_to_H2O=200.)
s6.backtrack_compositions(target_Fo=0.9, dm=0.005, Kd=0.32, verbose=False)
s6.compute_pressure_temperature(method="L09")
print("Dataset result with fixed Kd=0.32:P = 2.35 GPa, T = 1498 oC.")
print("Our result:                       P = %.2f GPa, T = %i oC." % 
    (s6.PT['P'], s6.PT['T']))

# ---- check reading in primaries
s7 = Suite("L09_s190_primary_fixed_Kd.csv", read_as_primary=True)
s7.compute_pressure_temperature(method="L09")
print("Our result from primary comp:     P = %.2f GPa, T = %i oC." % 
    (s7.PT['P'], s7.PT['T']))
print()
print("-----------------------------------------------------------")

# ---- Till et al, 2012 
s8 = Suite("TGK12_T6.csv", read_as_primary=True)
s8.compute_pressure_temperature(method="TGK12_SPL")
print("Till et al (2012) - Spinel is stable")
print()
print("Result from Table 6 for Reykjanes:P = 1.16 GPa, T = 1316 oC.")
print("Our result:                       P = %.2f GPa, T = %i oC." % 
    (s8.PT['P'][0], round(s8.PT['T'][0])))
print()
print("-----------------------------------------------------------")

# ---- Grove et al, 2013
s9 = Suite("G13_Table_5.csv", read_as_primary=True)
s9.compute_temperature(method="G13", P=4.)
print("Grove et al 2013")
print()
print("Result in Table 5 for KR4003")
print("experiement 40.07 (Walter, 1998): P = 4.00 GPa, T = 1624 oC")
print("Our result at 4 GPa:              P = %.2f GPa, T = %i oC." % 
    (s9.PT['P'], s9.PT['T']))
print()
print("-----------------------------------------------------------")

# ---- Herzberg & Asimow, 2015
s10 = Suite("HA15_MEGAPRIMELT3.csv", read_as_primary=True)
s10.compute_temperature(method="HA15", P=0.)
print("Herzberg & Asimow 2015")
print()
print("Result in Supp Info 1 for")
print("Mkea in MEGAPRIMELT3:             P = 0.00 GPa, T = 1431 oC")
print("Our result at 0 GPa:              P = %.2f GPa, T = %i oC." % 
    (s10.PT['P'][0], round(s10.PT['T'][0])))
print()
print("-----------------------------------------------------------")

# ---- Plank and Forsyth, 2016
s11 = Suite("PF16_S7.csv", src_FeIII_totFe=0.19, Ce_to_H2O=200.)
s11.backtrack_compositions(Kd=0.3, verbose=False)
s11.compute_pressure_temperature()
print("Plank & Forsyth (2016)")
print()
print("Result from Supplementary 7:      P = 2.09 GPa, T = 1347 oC.")
print("Our result:                       P = %.2f GPa, T = %i oC." % 
    (s11.PT['P'], s11.PT['T']))
# ---- check reading in primaries
s12 = Suite("PF16_S7_Primary.csv", read_as_primary=True)
s12.compute_pressure_temperature()
print("Our result from primary comp:     P = %.2f GPa, T = %i oC." % 
    (s12.PT['P'], s12.PT['T']))
print()
print("-----------------------------------------------------------")

# ---- Sun & Dasgupta, 2020
s13 = Suite("SD20_S1.csv", read_as_primary=True)
s13.compute_pressure_temperature(method="SD20")
print("Sun & Dasgupta (2020)")
print()
print("Result from Supplementary 1 for")
print("Longhi (1995) no known T or P:    P = 2.90 GPa, T = 1529 oC.")
print("Our result:                       P = %.2f GPa, T = %i oC." % 
    (s13.PT['P'], s13.PT['T']))
print()
print("-----------------------------------------------------------")

# ---- Krein et al, 2021
s14 = Suite("K21_KLB2.csv", read_as_primary=True)
s14.compute_pressure_temperature(method="BK21")
print("Brown Krien et al (2021)")
print()
print("- Spinel is stable")
print("Result in Table S1 for")
print("KLB-1 experiment 18")
print("(Hirose & Kushiro, 1993):         P = 1.51 GPa, T = 1326 oC")
print("Our result :                      P = %.2f GPa, T = %i oC." % 
    (s14.PT['P'], s14.PT['T']))
print()

s15 = Suite("K21_KR4003.csv", read_as_primary=True)
s15.compute_pressure_temperature(method="BK21")
print("- Garnet is stable")
print("Result in Table S1 for KR4003")
print("experiment 40.06 (Walter 1998):   P = 3.56 GPa, T = 1584 oC")
print("Our result :                      P = %.2f GPa, T = %i oC." % 
    (s15.PT['P'], s15.PT['T']))
print()
print("-----------------------------------------------------------")
