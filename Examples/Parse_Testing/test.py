from meltPT import *
import traceback

print("Sample with no SiO2:")
print()
try:
    s = Suite("UT09DV04_noSiO2.csv", src_FeIII_totFe=0.17)
except Exception:
    traceback.print_exc()
print()
print()

print("Sample with no MgO:")
print()
try:
    s = Suite("UT09DV04_noMgO.csv", src_FeIII_totFe=0.17)
except Exception:
    traceback.print_exc()
print()
print()


print("Sample with no Fe:")
print()
try:
    s = Suite("UT09DV04_noFe.csv", src_FeIII_totFe=0.17)
except Exception:
    traceback.print_exc()
print()
print()

print("Sample with no Al2O3 (also various other things):")
print()
s = Suite("UT09DV04_noAl2O3.csv", src_FeIII_totFe=0.17)
print()
print()

print("Sample with no H2O or Ce (also various other things):")
print()
s = Suite("UT09DV04_noH2O.csv", src_FeIII_totFe=0.17)
